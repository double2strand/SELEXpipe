#!/usr/bin/env python3
"""Filter FASTA reads by merged-construct length.

Default contract (safe for motif input):
  * exact-length retention (min = max = target)
  * short reads are dropped
  * no random ATCG is written unless ``--allow-pad`` is set (requires ``--seed``)

Expected merged-read lengths for common ligand designs (do not "correct"):

| Ligand design | Expected length |
|---------------|-----------------|
| lig147        | 101 bp          |
| lig200        | 154 bp          |
| ligN40        | 187 bp          |
| ligN70        | 247 bp          |

ligN40 vs ligN70 constant flanks differ; the 60 bp (not 30 bp) length
delta is expected.

Use ``--ligand`` to select a design, or ``--target_length`` for an explicit value.
"""

from __future__ import annotations

import argparse
import random
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Expected merged-read lengths for supported ligand designs.
# ligN40 and ligN70 use different constant flanks; Δ60 vs a 30 bp N-region
# difference is expected and must not be "fixed".
LIGAND_LENGTHS: dict[str, int] = {
    "lig147": 101,
    "lig200": 154,
    "ligN40": 187,
    "ligN70": 247,
}

LIGAND_FLANK_NOTE = (
    "ligN40 (187 bp) and ligN70 (247 bp) use different constant flanks; "
    "the 60 bp merged-length delta is expected (not a 30 bp N-region typo)."
)


def normalize_fasta_to_fixed_length(
    input_fasta: str,
    output_fasta: str,
    target_length: int,
    min_input_len: int,
    max_input_len: int,
    trim_from: str = "right",
    allow_pad: bool = False,
    seed: int | None = None,
) -> dict[str, int]:
    """Keep reads in ``[min_input_len, max_input_len]``, then trim or (opt-in) pad.

    Short reads are dropped unless ``allow_pad`` is true. Padding always
    requires ``seed``.
    """
    if target_length <= 0:
        raise ValueError(f"target_length must be positive, got {target_length}")
    if min_input_len <= 0:
        raise ValueError(f"min_input_len must be positive, got {min_input_len}")
    if max_input_len < min_input_len:
        raise ValueError(
            f"max_input_len ({max_input_len}) must be >= min_input_len ({min_input_len})"
        )
    if not (min_input_len <= target_length <= max_input_len):
        raise ValueError(
            f"target_length ({target_length}) must lie within "
            f"[{min_input_len}, {max_input_len}]"
        )
    if trim_from not in {"right", "left", "center"}:
        raise ValueError("--trim_from must be one of: right, left, center")
    if allow_pad and seed is None:
        raise ValueError("--allow-pad requires --seed so padding is reproducible")

    rng = random.Random(seed) if allow_pad else None

    n_total = 0
    n_kept = 0
    n_padded = 0
    n_trimmed = 0
    n_discarded = 0

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            n_total += 1
            sequence = str(record.seq).upper()
            length = len(sequence)

            if length < min_input_len or length > max_input_len:
                n_discarded += 1
                continue

            if length < target_length:
                if not allow_pad:
                    n_discarded += 1
                    continue
                assert rng is not None
                padding_length = target_length - length
                sequence = sequence + "".join(rng.choices("ATCG", k=padding_length))
                n_padded += 1
            elif length > target_length:
                if trim_from == "right":
                    sequence = sequence[:target_length]
                elif trim_from == "left":
                    sequence = sequence[-target_length:]
                else:  # center
                    extra = length - target_length
                    left_trim = extra // 2
                    sequence = sequence[left_trim : left_trim + target_length]
                n_trimmed += 1

            assert len(sequence) == target_length

            adjusted_record = SeqRecord(
                Seq(sequence),
                id=record.id,
                description=record.description,
            )
            SeqIO.write(adjusted_record, outfile, "fasta")
            n_kept += 1

    return {
        "total": n_total,
        "kept": n_kept,
        "padded": n_padded,
        "trimmed": n_trimmed,
        "discarded": n_discarded,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ligand_help = ", ".join(
        f"{name}={length} bp" for name, length in LIGAND_LENGTHS.items()
    )
    parser = argparse.ArgumentParser(
        description=(
            "Filter FASTA reads by length. Default keeps exact target length "
            "only and never pads. Use --allow-pad (with --seed) to pad shorts "
            "when the input window is opened below the target."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Ligand design expected merged lengths:\n"
            "  lig147  101 bp\n"
            "  lig200  154 bp\n"
            "  ligN40  187 bp\n"
            "  ligN70  247 bp\n"
            f"\n{LIGAND_FLANK_NOTE}\n"
        ),
    )
    parser.add_argument("input_fasta", type=str, help="Input FASTA file")
    parser.add_argument("output_fasta", type=str, help="Output FASTA file")
    length_group = parser.add_mutually_exclusive_group(required=True)
    length_group.add_argument(
        "--target_length",
        type=int,
        default=None,
        help="Final required read length (no hard-coded default)",
    )
    length_group.add_argument(
        "--ligand",
        choices=sorted(LIGAND_LENGTHS),
        default=None,
        help=f"Ligand design key used to set target length ({ligand_help})",
    )
    parser.add_argument(
        "--min_input_len",
        type=int,
        default=None,
        help=(
            "Minimum input read length to retain "
            "(default: same as --target_length)"
        ),
    )
    parser.add_argument(
        "--max_input_len",
        type=int,
        default=None,
        help=(
            "Maximum input read length to retain "
            "(default: same as --target_length)"
        ),
    )
    parser.add_argument(
        "--trim_from",
        choices=["right", "left", "center"],
        default="right",
        help="How to trim reads longer than target_length (default: right)",
    )
    parser.add_argument(
        "--allow-pad",
        action="store_true",
        help="Opt-in: pad reads shorter than target_length with random ATCG (requires --seed)",
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

    if args.ligand is not None:
        target_length = LIGAND_LENGTHS[args.ligand]
        ligand_label = args.ligand
    else:
        target_length = args.target_length
        ligand_label = None

    min_input_len = (
        args.min_input_len if args.min_input_len is not None else target_length
    )
    max_input_len = (
        args.max_input_len if args.max_input_len is not None else target_length
    )

    try:
        stats = normalize_fasta_to_fixed_length(
            input_fasta=args.input_fasta,
            output_fasta=args.output_fasta,
            target_length=target_length,
            min_input_len=min_input_len,
            max_input_len=max_input_len,
            trim_from=args.trim_from,
            allow_pad=args.allow_pad,
            seed=args.seed,
        )
    except Exception as exc:  # noqa: BLE001 - surface errors to CLI
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if ligand_label is not None:
        print(f"[OK] Ligand design:  {ligand_label}")
    print(f"[OK] Target length:  {target_length}")
    print(f"[OK] Input window:   {min_input_len}-{max_input_len}")
    print(f"[OK] Padding:        {'on' if args.allow_pad else 'off (shorts dropped)'}")
    print(f"[OK] Total reads:    {stats['total']}")
    print(f"[OK] Reads kept:     {stats['kept']}")
    print(f"[OK] Reads padded:   {stats['padded']}")
    print(f"[OK] Reads trimmed:  {stats['trimmed']}")
    print(f"[OK] Reads dropped:  {stats['discarded']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
