# Preprocessing

Turns paired-end HT-SELEX and NCAP-SELEX FASTQs into one fixed-width SEQ line per molecule.

[`run_preprocess.sh`](run_preprocess.sh) runs the four processes below, in order. [`preprocessing.md`](preprocessing.md) has the flag tables and those commands written out on their own.

## Processes

| Order | Process | What runs |
|-------|---------|-----------|
| 1 | [Merge lanes](01_merge_lanes/README.md) | `run_fastp.sh` once per lane, then `cat` of the two `.fq.gz` files |
| 2 | [FASTQ to FASTA](02_fastq_to_fasta/README.md) | `gzip -d -c`, then `seqtk seq -a`, inside the driver |
| 3 | [Length filter](03_length_filter/README.md) | `filter_by_length.py` |
| 4 | [SEQ extract](04_seq_extract/README.md) | `trim_to_kmer.py`, default one 40 bp window per read |

## Scripts in this directory

| File | Role |
|------|------|
| `run_preprocess.sh` | Driver for all four processes |
| `run_fastp.sh` | One lane of the fastp merge in process 1 |
| `filter_by_length.py` | Process 3 |
| `trim_to_kmer.py` | Process 4 |

`run_preprocess.sh` calls `filter_by_length.py` and `trim_to_kmer.py` from its own directory. It also sets the repository root to two levels above the script and opens `config/adapters.fa` from there. The numbered folders describe each process. The scripts stay here so those paths still resolve.

## Inputs

FASTQs come from [sample renaming](../00_sample_renaming/README.md). A Sanger folder such as `8030STDY17011914/` is now `HTS/SOX2R1/SOX2R1c1/` (assay, TF plus replicate, then cycle). Pass that cycle directory as `--input-dir` and `SOX2R1c1` as `--sample`. The read suffixes are `_R1` and `_R2`. The `R1` inside `SOX2R1` is the replicate from the sample name, not the read.

Four gzipped FASTQs per sample, under `--input-dir`:

| File | Role |
|------|------|
| `{SAMPLE}_L1_R1.fastq.gz` | lane 1, read 1 |
| `{SAMPLE}_L1_R2.fastq.gz` | lane 1, read 2 |
| `{SAMPLE}_L2_R1.fastq.gz` | lane 2, read 1 |
| `{SAMPLE}_L2_R2.fastq.gz` | lane 2, read 2 |

Pass exactly one of `--ligand` (`lig147`, `lig200`, `ligN40`, `ligN70`) or `--target-length`. The driver has no numeric default. The adapter file defaults to [`config/adapters.fa`](../../config/adapters.fa).

## Outputs

Written under `--output-dir/{SAMPLE}/`. The names below use lig200, target 154:

| File | Process |
|------|---------|
| `{SAMPLE}_L1.fq.gz`, `{SAMPLE}_L1.json`, `{SAMPLE}_L1.html` | lane 1 merge |
| `{SAMPLE}_L2.fq.gz`, `{SAMPLE}_L2.json`, `{SAMPLE}_L2.html` | lane 2 merge |
| `{SAMPLE}.fq.gz`, `{SAMPLE}.fq` | concatenated lanes, and the uncompressed copy |
| `{SAMPLE}.fa` | FASTA |
| `{SAMPLE}_154.fa` | length-filtered FASTA |
| `{SAMPLE}_154.seq` | one 40 bp line per kept read |

## Run

From the repository root:

```bash
workflow/01_preprocessing/run_preprocess.sh \
  --sample SOX2R1c1 \
  --ligand lig200 \
  --input-dir ./HTS/SOX2R1/SOX2R1c1 \
  --output-dir ./processed
```

`--skip-fastp` resumes from the lane-merged gzip files. `--kmer-start` or `--center-kmer` moves the extracted window. `--allow-pad --seed N` pads only a short k-mer window. The length filter pads shorter in-window reads on its own.

Tools: [fastp](https://github.com/opengene/fastp), seqtk, Python ≥ 3.9, Biopython, bash, gzip, and cat. Install notes: [envs/README.md](../../envs/README.md).
