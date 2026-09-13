#!/usr/bin/env python3
"""Prove the default preprocess → motif path cannot inject synthetic sequence.

Contract:
  * an exact-length read emits one SEQ line of the requested width
  * a short read is dropped
  * default output contains no random ATCG padding
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {rel}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


trim_to_kmer = load_module("trim_to_kmer", "preprocessing/trim_to_kmer.py")
filter_by_length = load_module("filter_by_length", "preprocessing/filter_by_length.py")

# Deterministic 101-mer (lig147 merged width). Blocks are unique so a leftover
# 21-mer cannot appear inside the first 40-mer (the old tile+pad artefact).
LIG147_SEQ = ("A" * 40) + ("C" * 40) + ("G" * 21)
assert len(LIG147_SEQ) == 101
SHORT_SEQ = "TTGACCAA"  # 8 bp; cannot fill a 40-mer


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")


def read_seq_lines(path: Path) -> list[str]:
    return [
        line.rstrip("\n\r")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


class TestTrimToKmerDefault(unittest.TestCase):
    def test_exact_length_emits_one_requested_width_line(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            seq_out = tmp_path / "out.seq"
            write_fasta(fasta, [("exact", LIG147_SEQ)])

            stats = trim_to_kmer.trim_to_kmers(fasta, seq_out, kmer_size=40)
            lines = read_seq_lines(seq_out)

            self.assertEqual(stats["written"], 1)
            self.assertEqual(stats["dropped"], 0)
            self.assertEqual(stats["padded"], 0)
            self.assertEqual(len(lines), 1)
            self.assertEqual(len(lines[0]), 40)
            self.assertEqual(lines[0], LIG147_SEQ[:40])

    def test_short_read_is_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            seq_out = tmp_path / "out.seq"
            write_fasta(
                fasta,
                [("exact", LIG147_SEQ), ("short", SHORT_SEQ)],
            )

            stats = trim_to_kmer.trim_to_kmers(fasta, seq_out, kmer_size=40)
            lines = read_seq_lines(seq_out)

            self.assertEqual(stats["total"], 2)
            self.assertEqual(stats["written"], 1)
            self.assertEqual(stats["dropped"], 1)
            self.assertEqual(lines, [LIG147_SEQ[:40]])
            self.assertNotIn(SHORT_SEQ, "".join(lines))

    def test_default_output_has_no_random_padding(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            seq_out = tmp_path / "out.seq"
            # 101 bp would produce a 21-base leftover under the old tile+pad path.
            write_fasta(fasta, [("lig147", LIG147_SEQ), ("short", SHORT_SEQ)])

            trim_to_kmer.trim_to_kmers(fasta, seq_out, kmer_size=40)
            lines = read_seq_lines(seq_out)

            self.assertEqual(len(lines), 1, "tiling would emit 3 lines for a 101-mer")
            observed = lines[0]
            self.assertEqual(observed, LIG147_SEQ[:40])
            # Every output base must come from the source window; none invented.
            self.assertTrue(observed in LIG147_SEQ)
            leftover = LIG147_SEQ[80:]
            self.assertEqual(leftover, "G" * 21)
            self.assertFalse(any("G" in line for line in lines))
            self.assertFalse(any("C" in line for line in lines))

    def test_allow_pad_requires_seed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            seq_out = tmp_path / "out.seq"
            write_fasta(fasta, [("short", SHORT_SEQ)])
            with self.assertRaises(ValueError):
                trim_to_kmer.trim_to_kmers(
                    fasta, seq_out, kmer_size=40, allow_pad=True, seed=None
                )

    def test_cli_default_matches_library(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            seq_out = tmp_path / "out.seq"
            write_fasta(fasta, [("exact", LIG147_SEQ), ("short", SHORT_SEQ)])
            rc = trim_to_kmer.main([str(fasta), str(seq_out), "--kmer-size", "40"])
            self.assertEqual(rc, 0)
            self.assertEqual(read_seq_lines(seq_out), [LIG147_SEQ[:40]])


class TestFilterByLengthDefault(unittest.TestCase):
    def test_exact_kept_short_dropped_no_pad(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            out = tmp_path / "out.fa"
            write_fasta(
                fasta,
                [("exact", LIG147_SEQ), ("short", SHORT_SEQ), ("long", LIG147_SEQ + "TTTT")],
            )
            stats = filter_by_length.normalize_fasta_to_fixed_length(
                str(fasta),
                str(out),
                target_length=101,
                min_input_len=101,
                max_input_len=101,
            )
            from Bio import SeqIO

            with out.open() as handle:
                kept = [str(rec.seq) for rec in SeqIO.parse(handle, "fasta")]
            self.assertEqual(stats["kept"], 1)
            self.assertEqual(stats["padded"], 0)
            self.assertEqual(stats["discarded"], 2)
            self.assertEqual(kept, [LIG147_SEQ])

    def test_open_window_without_allow_pad_still_drops_shorts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            out = tmp_path / "out.fa"
            write_fasta(fasta, [("short", SHORT_SEQ), ("exact", LIG147_SEQ)])
            stats = filter_by_length.normalize_fasta_to_fixed_length(
                str(fasta),
                str(out),
                target_length=101,
                min_input_len=1,
                max_input_len=101,
            )
            self.assertEqual(stats["padded"], 0)
            self.assertEqual(stats["kept"], 1)
            self.assertEqual(stats["discarded"], 1)

    def test_allow_pad_requires_seed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            out = tmp_path / "out.fa"
            write_fasta(fasta, [("short", SHORT_SEQ)])
            with self.assertRaises(ValueError):
                filter_by_length.normalize_fasta_to_fixed_length(
                    str(fasta),
                    str(out),
                    target_length=101,
                    min_input_len=1,
                    max_input_len=101,
                    allow_pad=True,
                    seed=None,
                )

    def test_ligand_lengths_and_flank_delta_unchanged(self) -> None:
        self.assertEqual(
            filter_by_length.LIGAND_LENGTHS,
            {"lig147": 101, "lig200": 154, "ligN40": 187, "ligN70": 247},
        )
        self.assertEqual(
            filter_by_length.LIGAND_LENGTHS["ligN70"] - filter_by_length.LIGAND_LENGTHS["ligN40"],
            60,
        )
        self.assertIn("flanks", filter_by_length.LIGAND_FLANK_NOTE.lower())


class TestMotifDriverWidthGuard(unittest.TestCase):
    def _run_motif(self, background: Path, signal: Path, seq_length: int = 40) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                "bash",
                str(ROOT / "motif" / "run_motif.sh"),
                "--background",
                str(background),
                "--signal",
                str(signal),
                "--seq-length",
                str(seq_length),
            ],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_refuses_mismatched_seq_width(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            background = tmp_path / "bg.seq"
            signal = tmp_path / "sig.seq"
            background.write_text(LIG147_SEQ[:40] + "\n", encoding="utf-8")
            # 21-base leftover as its own SEQ line (the old tile+pad tail).
            signal.write_text(LIG147_SEQ[:40] + "\n" + LIG147_SEQ[80:] + "\n", encoding="utf-8")
            proc = self._run_motif(background, signal, seq_length=40)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("length 21 != expected 40", proc.stderr)

    def test_refuses_empty_seq(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            background = tmp_path / "bg.seq"
            signal = tmp_path / "sig.seq"
            background.write_text(LIG147_SEQ[:40] + "\n", encoding="utf-8")
            signal.write_text("", encoding="utf-8")
            proc = self._run_motif(background, signal, seq_length=40)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("empty", proc.stderr)

    def test_matching_width_passes_guard(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            background = tmp_path / "bg.seq"
            signal = tmp_path / "sig.seq"
            background.write_text(LIG147_SEQ[:40] + "\n", encoding="utf-8")
            signal.write_text((LIG147_SEQ[10:50] + "\n") * 2, encoding="utf-8")
            proc = self._run_motif(background, signal, seq_length=40)
            # Width check passed; driver should then fail on missing spacek40,
            # not on SEQ length.
            self.assertNotIn("!= expected", proc.stderr)
            self.assertNotIn("SEQ is empty", proc.stderr)
            self.assertIn("spacek40", proc.stderr)


if __name__ == "__main__":
    raise SystemExit(unittest.main(verbosity=2))
