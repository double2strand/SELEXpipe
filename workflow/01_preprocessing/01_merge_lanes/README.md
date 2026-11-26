# Merge lanes

First process in [preprocessing](../README.md).

[fastp](https://github.com/opengene/fastp) quality-filters, adapter-trims, and overlap-merges each paired-end lane. The driver then concatenates the two lane files into one sample FASTQ.

## Commands, in order

[`../run_preprocess.sh`](../run_preprocess.sh) runs these unless `--skip-fastp` is set.

1. Lane 1. The function `run_fastp_lane` in the driver, or [`../run_fastp.sh`](../run_fastp.sh) on its own.
2. Lane 2. The same command.
3. `cat` of the two gzipped lane files into `{SAMPLE}.fq.gz`. The `cat` is only in the driver. `run_fastp.sh` stops after one lane.

`run_fastp.sh` takes `in1 in2 out_prefix [adapters] [threads]`.

## Inputs

Under `--input-dir`:

| File | Role |
|------|------|
| `{SAMPLE}_L1_R1.fastq.gz` | lane 1, read 1 |
| `{SAMPLE}_L1_R2.fastq.gz` | lane 1, read 2 |
| `{SAMPLE}_L2_R1.fastq.gz` | lane 2, read 1 |
| `{SAMPLE}_L2_R2.fastq.gz` | lane 2, read 2 |

Adapter FASTA: [`config/adapters.fa`](../../../config/adapters.fa), unless `--adapter-fasta` (driver) or the fourth argument (`run_fastp.sh`) is set. The sequences in that file are example adapters.

## Outputs

Under `--output-dir/{SAMPLE}/`:

| File | Step |
|------|------|
| `{SAMPLE}_L1.fq.gz`, `{SAMPLE}_L1.json`, `{SAMPLE}_L1.html` | fastp, lane 1 |
| `{SAMPLE}_L2.fq.gz`, `{SAMPLE}_L2.json`, `{SAMPLE}_L2.html` | fastp, lane 2 |
| `{SAMPLE}.fq.gz` | both lanes concatenated; the gzipped file is kept |

## fastp flags

`run_fastp.sh` uses the same flags as `run_fastp_lane`. The full table is in [preprocessing.md](../preprocessing.md).

| Flag | Value | Role |
|------|-------|------|
| `-q` | 10 | Qualified-base quality cutoff (Phred) |
| `--in1`, `--in2` | the lane pair | Paired reads |
| `-m`, `--merged_out` | | Merge overlapping pairs |
| `--adapter_fasta` | adapter FASTA | Adapter sequences |
| `--json`, `--html` | | QC reports |
| `--overlap_len_require` | 5 | Minimum overlap for merging |
| `-g`, `-x` | on | Trim polyG and polyX tails |
| `--length_required` | 15 | Drop merged reads shorter than 15 bp |
| `--n_base_limit` | 10 | Drop reads with more than 10 `N` bases |
| `-y`, `--complexity_threshold` | 5 | Low-complexity filter |
| `--thread` | `nproc`, or 4 | Worker threads |

From the repository root:

```bash
workflow/01_preprocessing/run_fastp.sh \
  raw/SAMPLE_L1_R1.fastq.gz \
  raw/SAMPLE_L1_R2.fastq.gz \
  processed/SAMPLE/SAMPLE_L1

workflow/01_preprocessing/run_fastp.sh \
  raw/SAMPLE_L2_R1.fastq.gz \
  raw/SAMPLE_L2_R2.fastq.gz \
  processed/SAMPLE/SAMPLE_L2
```

`--help` on `run_fastp.sh` still prints `scripts/01_preprocess/run_fastp.sh`. The commands above match this tree.

Next process: [FASTQ to FASTA](../02_fastq_to_fasta/README.md).
