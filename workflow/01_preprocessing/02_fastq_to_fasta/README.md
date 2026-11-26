# FASTQ to FASTA

Second process in [preprocessing](../README.md). It runs after [lane merge](../01_merge_lanes/README.md) and before the [length filter](../03_length_filter/README.md).

This process has no script of its own. [`../run_preprocess.sh`](../run_preprocess.sh) runs both commands.

## Commands, in order

1. `gzip -d -f -c` writes an uncompressed copy and leaves `{SAMPLE}.fq.gz` in place.
2. `seqtk seq -a` converts that copy to FASTA.

## Inputs

| File | Role |
|------|------|
| `{SAMPLE}.fq.gz` | Concatenated lane merges from the previous process |

## Outputs

| File | Role |
|------|------|
| `{SAMPLE}.fq` | Uncompressed copy. `seqtk` reads this file |
| `{SAMPLE}.fa` | FASTA passed to `filter_by_length.py` |

```bash
COMBINED_GZ=processed/SAMPLE/SAMPLE.fq.gz
COMBINED_FQ=processed/SAMPLE/SAMPLE.fq
FASTA=processed/SAMPLE/SAMPLE.fa

gzip -d -f -c "$COMBINED_GZ" > "$COMBINED_FQ"
seqtk seq -a "$COMBINED_FQ" > "$FASTA"
```

`gzip` and `seqtk` need to be on `PATH`. Install notes: [envs/README.md](../../../envs/README.md).

Next process: [length filter](../03_length_filter/README.md).
