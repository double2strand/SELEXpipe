#!/usr/bin/env python3
"""Checks for the k-mer extract, the length filter, and the motif width guard.

trim_to_kmer defaults:
  * an exact-length read emits one SEQ line of the requested width
  * a short read is dropped
  * default output contains no random ATCG padding

filter_by_length.py (the shipped script):
  * drops reads outside the input-length window
  * pads shorter in-window reads and trims longer ones to the target
  * --seed makes that padding repeatable
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


trim_to_kmer = load_module("trim_to_kmer", "scripts/01_preprocess/trim_to_kmer.py")
filter_by_length = load_module("filter_by_length", "scripts/01_preprocess/filter_by_length.py")

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


def read_fasta_seqs(path: Path) -> list[str]:
    from Bio import SeqIO

    with path.open() as handle:
        return [str(rec.seq) for rec in SeqIO.parse(handle, "fasta")]


class TestFilterByLength(unittest.TestCase):
    def test_window_keeps_bounds_and_drops_outside(self) -> None:
        inside_short = "A" * 145
        exact = "C" * 154
        inside_long = "G" * 163
        too_short = "T" * 144
        too_long = "A" * 164
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            out = tmp_path / "out.fa"
            write_fasta(
                fasta,
                [
                    ("too_short", too_short),
                    ("inside_short", inside_short),
                    ("exact", exact),
                    ("inside_long", inside_long),
                    ("too_long", too_long),
                ],
            )
            stats = filter_by_length.normalize_fasta_to_fixed_length(
                str(fasta),
                str(out),
                target_length=154,
                min_input_len=145,
                max_input_len=163,
            )
            kept = read_fasta_seqs(out)
            self.assertEqual(stats["total"], 5)
            self.assertEqual(stats["discarded"], 2)
            self.assertEqual(stats["kept"], 3)
            self.assertEqual(stats["padded"], 1)
            self.assertEqual(stats["trimmed"], 1)
            self.assertEqual(len(kept), 3)
            self.assertTrue(all(len(seq) == 154 for seq in kept))
            self.assertTrue(kept[0].startswith(inside_short))
            self.assertEqual(kept[1], exact)
            self.assertEqual(kept[2], inside_long[:154])

    def test_padding_appends_random_bases(self) -> None:
        original = "A" * 150
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            out = tmp_path / "out.fa"
            write_fasta(fasta, [("short_in_window", original)])
            stats = filter_by_length.normalize_fasta_to_fixed_length(
                str(fasta), str(out), target_length=154, min_input_len=145, max_input_len=163
            )
            kept = read_fasta_seqs(out)
            self.assertEqual(stats["padded"], 1)
            self.assertEqual(stats["trimmed"], 0)
            self.assertEqual(len(kept), 1)
            self.assertEqual(len(kept[0]), 154)
            self.assertTrue(kept[0].startswith(original))
            self.assertTrue(set(kept[0][150:]).issubset(set("ATCG")))

    def test_trim_right_left_and_center(self) -> None:
        original = "ACGT" * 40  # 160 bp
        self.assertEqual(len(original), 160)
        modes = {
            "right": original[:154],
            "left": original[-154:],
            "center": original[3:157],  # extra = 6, left trim = 3
        }
        for mode, expected in modes.items():
            with self.subTest(mode=mode):
                with tempfile.TemporaryDirectory() as tmp:
                    tmp_path = Path(tmp)
                    fasta = tmp_path / "in.fa"
                    out = tmp_path / "out.fa"
                    write_fasta(fasta, [("long", original)])
                    stats = filter_by_length.normalize_fasta_to_fixed_length(
                        str(fasta),
                        str(out),
                        target_length=154,
                        min_input_len=145,
                        max_input_len=163,
                        trim_from=mode,
                    )
                    kept = read_fasta_seqs(out)
                    self.assertEqual(stats["trimmed"], 1)
                    self.assertEqual(stats["padded"], 0)
                    self.assertEqual(kept, [expected])
                    self.assertEqual(len(expected), 154)

    def test_seed_reproducible_via_cli(self) -> None:
        script = ROOT / "scripts" / "01_preprocess" / "filter_by_length.py"
        original = "A" * 150

        def run(seed: str, out_name: str, tmp_path: Path) -> str:
            fasta = tmp_path / "in.fa"
            out = tmp_path / out_name
            write_fasta(fasta, [("short_in_window", original)])
            proc = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    str(fasta),
                    str(out),
                    "--target_length",
                    "154",
                    "--min_input_len",
                    "145",
                    "--max_input_len",
                    "163",
                    "--seed",
                    seed,
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("[OK] Reads padded:", proc.stdout)
            self.assertIn("[OK] Reads dropped:", proc.stdout)
            seqs = read_fasta_seqs(out)
            self.assertEqual(len(seqs), 1)
            self.assertTrue(seqs[0].startswith(original))
            self.assertEqual(len(seqs[0]), 154)
            return seqs[0]

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            first = run("1", "a.fa", tmp_path)
            second = run("1", "b.fa", tmp_path)
            other = run("2", "c.fa", tmp_path)
            self.assertEqual(first, second)
            self.assertNotEqual(first, other)

    def test_filtered_154_read_yields_one_40mer(self) -> None:
        # 150 bp inside the default window. Padding is on the right, so the
        # default 40-mer is the original 5' end, then one line from trim_to_kmer.
        original = ("A" * 40) + ("C" * 110)
        self.assertEqual(len(original), 150)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "in.fa"
            filtered = tmp_path / "out.fa"
            seq_out = tmp_path / "out.seq"
            write_fasta(fasta, [("lig200ish", original)])
            stats = filter_by_length.normalize_fasta_to_fixed_length(
                str(fasta),
                str(filtered),
                target_length=154,
                min_input_len=145,
                max_input_len=163,
            )
            self.assertEqual(stats["kept"], 1)
            kept = read_fasta_seqs(filtered)
            self.assertEqual(len(kept[0]), 154)
            kmer_stats = trim_to_kmer.trim_to_kmers(filtered, seq_out, kmer_size=40)
            lines = read_seq_lines(seq_out)
            self.assertEqual(kmer_stats["written"], 1)
            self.assertEqual(kmer_stats["padded"], 0)
            self.assertEqual(kmer_stats["dropped"], 0)
            self.assertEqual(lines, [original[:40]])


class TestLigandWindowAssumption(unittest.TestCase):
    def test_driver_maps_ligand_lengths_with_nine_bp_window(self) -> None:
        text = (ROOT / "scripts" / "01_preprocess" / "run_preprocess.sh").read_text(encoding="utf-8")
        self.assertIn("lig147) echo 101", text)
        self.assertIn("lig200) echo 154", text)
        self.assertIn("ligN40) echo 187", text)
        self.assertIn("ligN70) echo 247", text)
        self.assertIn("TARGET_LEN - 9", text)
        self.assertIn("TARGET_LEN + 9", text)
        self.assertNotIn("--ligand", text.split("FILTER_ARGS=", 1)[-1])


class TestMotifDriverWidthGuard(unittest.TestCase):
    def _run_motif(self, background: Path, signal: Path, seq_length: int = 40) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                "bash",
                str(ROOT / "scripts" / "02_motif" / "run_motif.sh"),
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
