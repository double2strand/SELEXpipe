#!/usr/bin/env python3
"""Sanger-ID rename: dry-run, final paths, and a clash that changes nothing."""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "00_sample_renaming" / "selex_rename.py"
EXAMPLE = ROOT / "workflow" / "00_sample_renaming" / "examples" / "example_map.tsv"

# A few rows of examples/example_map.tsv, and where selex_rename.py puts them.
WANTED = {
    "8030STDY17011914": "HTS/SOX2R1/SOX2R1c1",
    "8030STDY17012491": "NCAP/GATA1_KLF1R1/GATA1_KLF1R1c1",
}


def example_rows() -> list[tuple[str, str]]:
    found: dict[str, str] = {}
    for line in EXAMPLE.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tag, name = line.split(None, 1)
        if tag in WANTED:
            found[tag] = name
    return [(tag, found[tag]) for tag in WANTED]


def write_fastqs(base: Path, tags: list[str]) -> None:
    for tag in tags:
        folder = base / tag
        folder.mkdir(parents=True)
        for lane, read in (("L1", "R1"), ("L2", "R2")):
            (folder / f"{tag}_{lane}_{read}.fastq.gz").write_bytes(b"fq")
        (folder / "notes.txt").write_bytes(b"note")


def files_under(base: Path) -> list[str]:
    return sorted(path.relative_to(base).as_posix() for path in base.rglob("*") if path.is_file())


def run_rename(base: Path, map_path: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["python3", str(SCRIPT), "--base-dir", str(base), "--map", str(map_path), *extra],
        check=False,
        capture_output=True,
        text=True,
    )


class TestSampleRenaming(unittest.TestCase):
    def test_example_map_has_the_rows_this_test_uses(self) -> None:
        rows = dict(example_rows())
        self.assertEqual(set(rows), set(WANTED))
        self.assertEqual(rows["8030STDY17011914"], "C1HTS200SOX2R1")
        self.assertEqual(rows["8030STDY17012491"], "C1NCAP200GATA1_KLF1R1")

    def test_dry_run_leaves_tree_and_real_run_builds_layout(self) -> None:
        rows = example_rows()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            base = tmp_path / "run"
            map_path = tmp_path / "map.tsv"
            map_path.write_text(
                "".join(f"{tag}\t{name}\n" for tag, name in rows),
                encoding="utf-8",
            )
            write_fastqs(base, [tag for tag, _ in rows])
            before = files_under(base)

            dry = run_rename(base, map_path, "--dry-run")
            self.assertEqual(dry.returncode, 0, dry.stderr)
            self.assertEqual(files_under(base), before)

            done = run_rename(base, map_path)
            self.assertEqual(done.returncode, 0, done.stderr)
            expected = []
            for rel in WANTED.values():
                stem = rel.rsplit("/", 1)[-1]
                expected.append(f"{rel}/{stem}_L1_R1.fastq.gz")
                expected.append(f"{rel}/{stem}_L2_R2.fastq.gz")
                expected.append(f"{rel}/notes.txt")
            self.assertEqual(files_under(base), sorted(expected))

    def test_flat_hts_ncap_pair_exits_2_and_changes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            base = tmp_path / "run"
            map_path = tmp_path / "map.tsv"
            map_path.write_text(
                "AAA\tC1HTS200FOOR1\nBBB\tC1NCAP200FOOR1\n",
                encoding="utf-8",
            )
            write_fastqs(base, ["AAA", "BBB"])
            before = files_under(base)

            proc = run_rename(base, map_path, "--no-assay-subdir")
            self.assertEqual(proc.returncode, 2, proc.stdout + proc.stderr)
            self.assertIn("nothing was changed", proc.stderr)
            self.assertEqual(files_under(base), before)


if __name__ == "__main__":
    raise SystemExit(unittest.main(verbosity=2))
