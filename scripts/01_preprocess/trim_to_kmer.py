#!/usr/bin/env python3
"""Extract one fixed-width window per FASTA record into a SEQ file.

Default contract (safe for motif input):
  * one k-mer per molecule (``--start`` or ``--center``)
  * short reads / windows that cannot fill ``--kmer-size`` are dropped
  * no random ATCG is written unless ``--allow-pad`` is set (requires ``--seed``)

Tiling a construct into consecutive k-mers with fake tails is not supported.
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path


def iter_fasta(path: Path):
    header: str | None = None
    seq_chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks).replace(" ", "").replace("\t", "")
                header = line[1:].strip() or line
                seq_chunks = []
            else:
                if header is None:
                    raise ValueError(f"Malformed FASTA in {path}: sequence before header")
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks).replace(" ", "").replace("\t", "")


def _pad_right(seq: str, length: int, rng: random.Random) -> str:
    missing = length - len(seq)
    if missing <= 0:
        return seq
    return seq + "".join(rng.choices("ATCG", k=missing))


def extract_window(
    sequence: str,
    kmer_size: int,
    start: int,
    center: bool,
    allow_pad: bool,
    rng: random.Random | None,
) -> str | None:
    """Return one k-mer or ``None`` if the read should be dropped."""
    seq = sequence.upper()
    if kmer_size <= 0:
        raise ValueError(f"k-mer size must be positive, got {kmer_size}")
    if start < 0:
        raise ValueError(f"start offset must be >= 0, got {start}")

    if center:
        if len(seq) >= kmer_size:
            window_start = (len(seq) - kmer_size) // 2
            return seq[window_start : window_start + kmer_size]
        if allow_pad:
            if rng is None:
                raise ValueError("rng is required when allow_pad is set")
            return _pad_right(seq, kmer_size, rng)
        return None

    window_end = start + kmer_size
    if start < len(seq) and window_end <= len(seq):
        return seq[start:window_end]
    if allow_pad and start < len(seq):
        if rng is None:
            raise ValueError("rng is required when allow_pad is set")
        return _pad_right(seq[start:], kmer_size, rng)
    return None


def trim_to_kmers(
    input_fasta: Path,
    output_path: Path,
    kmer_size: int,
    start: int = 0,
    center: bool = False,
    as_fasta: bool = False,
    allow_pad: bool = False,
    seed: int | None = None,
) -> dict[str, int]:
    if allow_pad and seed is None:
        raise ValueError("--allow-pad requires --seed so padding is reproducible")

    rng = random.Random(seed) if allow_pad else None

    total = 0
    written = 0
    dropped = 0
    padded = 0

    with output_path.open("w", encoding="utf-8") as out:
        for header, sequence in iter_fasta(input_fasta):
            total += 1
            raw_len = len(sequence.upper().replace(" ", "").replace("\t", ""))
            kmer = extract_window(sequence, kmer_size, start, center, allow_pad, rng)
            if kmer is None:
                dropped += 1
                continue
            if raw_len < kmer_size or (not center and start + kmer_size > raw_len):
                padded += 1
            if as_fasta:
                out.write(f">{header}\n")
                for i in range(0, len(kmer), 80):
                    out.write(kmer[i : i + 80] + "\n")
            else:
                out.write(kmer + "\n")
            written += 1

    return {
        "total": total,
        "written": written,
        "dropped": dropped,
        "padded": padded,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract one fixed-width k-mer per FASTA record. Short reads are "
            "dropped unless --allow-pad is set."
        )
    )
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("output_path", type=Path, help="Output .seq (or FASTA) path")
    parser.add_argument(
        "-k",
        "--kmer-size",
        "--length",
        dest="kmer_size",
        type=int,
        default=40,
        help="K-mer / SEQ line length to extract (default: 40)",
    )
    start_group = parser.add_mutually_exclusive_group()
    start_group.add_argument(
        "-s",
        "--start",
        type=int,
        default=0,
        help="0-based start offset for the k-mer window (default: 0)",
    )
    start_group.add_argument(
        "--center",
        action="store_true",
        help="Center the k-mer window within each sequence",
    )
    parser.add_argument(
        "--fasta",
        action="store_true",
        help="Write FASTA output instead of one sequence per line",
    )
    parser.add_argument(
        "--allow-pad",
        action="store_true",
        help="Opt-in: pad a short window on the right with random ATCG (requires --seed)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="RNG seed required when --allow-pad is set",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.input_fasta.is_file():
        print(f"error: input FASTA not found: {args.input_fasta}", file=sys.stderr)
        return 1

    try:
        stats = trim_to_kmers(
            args.input_fasta,
            args.output_path,
            args.kmer_size,
            args.start,
            args.center,
            args.fasta,
            args.allow_pad,
            args.seed,
        )
    except Exception as exc:  # noqa: BLE001 - surface parse/IO errors to CLI
        print(f"error: {exc}", file=sys.stderr)
        return 1

    mode = "centered" if args.center else f"start={args.start}"
    pad_note = "pad=on" if args.allow_pad else "pad=off"
    print(
        f"trim_to_kmer: wrote {stats['written']}/{stats['total']} "
        f"{args.kmer_size}-mers ({mode}, {pad_note}; "
        f"dropped {stats['dropped']}, padded {stats['padded']}) "
        f"-> {args.output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
