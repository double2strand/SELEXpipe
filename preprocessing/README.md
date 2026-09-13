# Preprocessing

This module converts raw paired-end SELEX FASTQ reads into analysis-ready
fixed-width SEQ files for motif discovery.

**Default contract (safe):** one variable-region extract per molecule, exact-length
merged-read filter, short reads dropped, **no random ATCG** in motif input.
Padding exists only behind `--allow-pad` (requires `--seed`).

## What it does

For each sample (two sequencing lanes, paired-end):

1. **Adapter trim + read merge (`fastp`)** — quality filter, TruSeq (or other) adapter removal, and overlap-based merging of R1/R2 into a single read per fragment, independently for lane 1 and lane 2.
2. **Lane concatenation** — concatenate the two lane-merged gzipped FASTQ files into one sample-level FASTQ.gz.
3. **Decompress** — write an uncompressed FASTQ **copy** with `gzip -d -f -c` (the `.fq.gz` is kept). Do not use in-place `gzip -d` here.
4. **FASTQ → FASTA (`seqtk`)** — convert the combined FASTQ to FASTA.
5. **Length filter (`filter_by_length.py`)** — keep reads whose length equals `--target_length` (or the length implied by `--ligand`). Shorts are dropped. Opening `--min-input-len` / `--max-input-len` still does **not** pad unless `--allow-pad` is set.
6. **One-window extract (`trim_to_kmer.py`)** — write **one** SEQ line per retained molecule (`--start` or `--center`, default **40 bp**). Short windows are dropped. This is not a tiled 40-mer library.

## Ligand design expected lengths

These are **merged construct** lengths. They are correct; do not "fix" them.

| Ligand design | Expected merged length |
|---------------|------------------------|
| lig147        | 101 bp                 |
| lig200        | 154 bp                 |
| ligN40        | 187 bp                 |
| ligN70        | 247 bp                 |

**Flanks differ:** ligN40 vs ligN70 constant flanks are not the same. The
merged-length delta is **60 bp, not 30 bp**. That is expected (N-region 40 vs 70
plus different flanks). Do not change 187 / 247 to force a 30 bp gap.

Use `--ligand lig147` (etc.) as a shorthand for the corresponding `--target-length`, or pass an explicit length when needed.

## Generic naming convention

Replace experiment-specific IDs with a sample prefix `SAMPLE`. Expected inputs under `--input-dir`:

| Role | Filename |
|------|----------|
| Lane 1, read 1 | `{SAMPLE}_L1_R1.fastq.gz` |
| Lane 1, read 2 | `{SAMPLE}_L1_R2.fastq.gz` |
| Lane 2, read 1 | `{SAMPLE}_L2_R1.fastq.gz` |
| Lane 2, read 2 | `{SAMPLE}_L2_R2.fastq.gz` |

Outputs written under `--output-dir/{SAMPLE}/` (example with lig147 / 101 bp):

| Stage | Filename |
|-------|----------|
| Lane 1 merged | `{SAMPLE}_L1.fq.gz` |
| Lane 2 merged | `{SAMPLE}_L2.fq.gz` |
| Combined lanes | `{SAMPLE}.fq.gz`, `{SAMPLE}.fq` |
| FASTA | `{SAMPLE}.fa` |
| Length-filtered FASTA | `{SAMPLE}_{N}.fa` |
| One k-mer per read (plain text) | `{SAMPLE}_{N}.seq` |

Adapter FASTA is supplied separately (e.g. `Adapters_TruSeq3PE.fa`).

## Dependencies

- [`fastp`](https://github.com/OpenGene/fastp)
- [`seqtk`](https://github.com/lh3/seqtk)
- `gzip`, `bash`, `python3` (≥ 3.9 recommended)
- [Biopython](https://biopython.org/) (`pip install -r requirements.txt`)

## Usage

```bash
./preprocessing/run_preprocess.sh \
  --sample SAMPLE \
  --adapter-fasta ./Adapters_TruSeq3PE.fa \
  --ligand lig147 \
  --input-dir ./raw \
  --output-dir ./processed
```

Provide exactly one of `--ligand` or `--target-length` (no hard-coded default such as 154 or 101). The default retention window is exact (`min = max = target`). To extract a window that is not the leftmost k bases, pass `--kmer-start` (skip a 5′ flank) or `--center-kmer`.

Opt-in padding (not the default; requires a seed):

```bash
./preprocessing/run_preprocess.sh \
  --sample SAMPLE \
  --adapter-fasta ./Adapters_TruSeq3PE.fa \
  --ligand lig200 \
  --min-input-len 145 \
  --max-input-len 163 \
  --allow-pad \
  --seed 1 \
  --input-dir ./raw \
  --output-dir ./processed
```

Equivalent standalone steps (mirroring the driver; generic names):

```bash
fastp -q 10 \
  --in1 SAMPLE_L1_R1.fastq.gz --in2 SAMPLE_L1_R2.fastq.gz \
  -m --merged_out SAMPLE_L1.fq.gz \
  --adapter_fasta=./Adapters_TruSeq3PE.fa \
  --json=./SAMPLE_L1.json --html=./SAMPLE_L1.html \
  --overlap_len_require=5 -gx --length_required=15 \
  --n_base_limit=10 -y --complexity_threshold=5

fastp -q 10 \
  --in1 SAMPLE_L2_R1.fastq.gz --in2 SAMPLE_L2_R2.fastq.gz \
  -m --merged_out SAMPLE_L2.fq.gz \
  --adapter_fasta=./Adapters_TruSeq3PE.fa \
  --json=./SAMPLE_L2.json --html=./SAMPLE_L2.html \
  --overlap_len_require=5 -gx --length_required=15 \
  --n_base_limit=10 -y --complexity_threshold=5

cat SAMPLE_L1.fq.gz SAMPLE_L2.fq.gz > SAMPLE.fq.gz
gzip -d -f -c SAMPLE.fq.gz > SAMPLE.fq
seqtk seq -a SAMPLE.fq > SAMPLE.fa
python3 preprocessing/filter_by_length.py SAMPLE.fa SAMPLE_101.fa --ligand lig147
python3 preprocessing/trim_to_kmer.py SAMPLE_101.fa SAMPLE_101.seq --kmer-size 40
```

### Helper scripts

`filter_by_length.py` keeps FASTA records at a resolved target length. Provide either `--target_length` or `--ligand` (mutually exclusive; no hard-coded default). If `--min_input_len` / `--max_input_len` are omitted, both default to the target (exact-length retention). Shorts are **dropped**. Random ATCG padding happens only with `--allow-pad --seed`. Longer reads in an opened window are trimmed (`--trim_from right|left|center`).

`trim_to_kmer.py` writes **one** SEQ line per record. `--length` is accepted as an alias for `--kmer-size`. There is no chunk/tile mode.

```bash
# Exact-length filter via ligand design key
python3 preprocessing/filter_by_length.py input.fa output_101.fa --ligand lig147

# Or set an explicit length
python3 preprocessing/filter_by_length.py input.fa output_154.fa --target_length 154

# Windowed retention; shorts still dropped unless --allow-pad --seed
python3 preprocessing/filter_by_length.py input.fa output_154.fa \
  --ligand lig200 \
  --min_input_len 145 \
  --max_input_len 163 \
  --trim_from right \
  --allow-pad \
  --seed 1

# One 40-base window per read (default: start=0). Shorts dropped, no padding.
python3 preprocessing/trim_to_kmer.py input_101.fa output.seq --kmer-size 40

# Centered window
python3 preprocessing/trim_to_kmer.py input_101.fa output.seq --kmer-size 40 --center

# Opt-in short-window padding (seed required)
python3 preprocessing/trim_to_kmer.py input.fa output.seq --kmer-size 40 --allow-pad --seed 1
```

## `fastp` parameters (current defaults)

| Flag | Value | Role |
|------|-------|------|
| `-q` | 10 | Base quality threshold |
| `-m` / `--merged_out` | — | Merge overlapping PE reads |
| `--adapter_fasta` | user-supplied | Adapter sequences |
| `--overlap_len_require` | 5 | Minimum overlap for merging |
| `-g` / `-x` | on | PolyG trimming / length filtering helpers as configured by fastp |
| `--length_required` | 15 | Discard very short merges |
| `--n_base_limit` | 10 | Max ambiguous bases |
| `-y` / `--complexity_threshold` | 5 | Low-complexity filter |

These defaults match the established SELEX preprocessing recipe used for this project and can be parameterised further as the pipeline evolves.
