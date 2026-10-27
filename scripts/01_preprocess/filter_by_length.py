#!/usr/bin/env python3

import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def normalize_fasta_to_fixed_length(
    input_fasta,
    output_fasta,
    target_length=154,
    min_input_len=145,
    max_input_len=163,
    trim_from="right",
):
    """
    Keep only reads within an allowed input-length window, then normalize all kept reads
    to exactly target_length by padding shorter reads or trimming longer reads.

    Args:
        input_fasta (str): Input FASTA file
        output_fasta (str): Output FASTA file
        target_length (int): Final required length for all retained reads
        min_input_len (int): Minimum input length to retain
        max_input_len (int): Maximum input length to retain
        trim_from (str): 'right', 'left', or 'center'
    """
    bases = ["A", "T", "C", "G"]

    n_total = 0
    n_kept = 0
    n_padded = 0
    n_trimmed = 0
    n_discarded = 0

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            n_total += 1
            sequence = str(record.seq).upper()
            L = len(sequence)

            # Discard reads outside allowed input range
            if L < min_input_len or L > max_input_len:
                n_discarded += 1
                continue

            # Pad shorter reads
            if L < target_length:
                padding_length = target_length - L
                sequence = sequence + "".join(random.choices(bases, k=padding_length))
                n_padded += 1

            # Trim longer reads
            elif L > target_length:
                if trim_from == "right":
                    sequence = sequence[:target_length]
                elif trim_from == "left":
                    sequence = sequence[-target_length:]
                elif trim_from == "center":
                    extra = L - target_length
                    left_trim = extra // 2
                    sequence = sequence[left_trim:left_trim + target_length]
                else:
                    raise ValueError("--trim_from must be one of: right, left, center")
                n_trimmed += 1

            # Keep exact-length reads unchanged
            assert len(sequence) == target_length

            adjusted_record = SeqRecord(
                Seq(sequence),
                id=record.id,
                description=record.description
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter FASTA reads by input-length window, then normalize all retained reads to exactly 154 bp."
    )
    parser.add_argument("input_fasta", type=str, help="Input FASTA file")
    parser.add_argument("output_fasta", type=str, help="Output FASTA file")
    parser.add_argument("--target_length", type=int, default=154, help="Final required read length (default: 154)")
    parser.add_argument("--min_input_len", type=int, default=145, help="Minimum input read length to retain (default: 145)")
    parser.add_argument("--max_input_len", type=int, default=163, help="Maximum input read length to retain (default: 163)")
    parser.add_argument("--trim_from", choices=["right", "left", "center"], default="right",
                        help="How to trim reads longer than target_length (default: right)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducible padding")

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    stats = normalize_fasta_to_fixed_length(
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta,
        target_length=args.target_length,
        min_input_len=args.min_input_len,
        max_input_len=args.max_input_len,
        trim_from=args.trim_from,
    )

    print(f"[OK] Total reads:     {stats['total']}")
    print(f"[OK] Reads kept:     {stats['kept']}")
    print(f"[OK] Reads padded:   {stats['padded']}")
    print(f"[OK] Reads trimmed:  {stats['trimmed']}")
    print(f"[OK] Reads dropped:  {stats['discarded']}")
