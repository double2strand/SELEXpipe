# Preprocessing

Paths below are relative to the repository root. Scripts live in
`scripts/01_preprocess/`.

This module converts raw paired-end SELEX FASTQ reads into analysis-ready
fixed-width SEQ files for motif discovery.

**Default contract:** one variable-region extract per molecule. The length filter
keeps an input-length window and normalizes every kept read to the target by
padding (random A/T/C/G on the right) or trimming. The k-mer step still drops a
short window unless `--allow-pad` and `--seed` are set.

## What it does

For each sample (two sequencing lanes, paired-end):

1. **Adapter trim + read merge (`fastp`)** — quality filter, TruSeq (or other) adapter removal, and overlap-based merging of R1/R2 into a single read per fragment, independently for lane 1 and lane 2.
2. **Lane concatenation** — concatenate the two lane-merged gzipped FASTQ files into one sample-level FASTQ.gz.
3. **Decompress** — write an uncompressed FASTQ **copy** with `gzip -d -f -c` (the `.fq.gz` is kept). Do not use in-place `gzip -d` here.
4. **FASTQ → FASTA (`seqtk`)** — convert the combined FASTQ to FASTA.
5. **Length filter (`filter_by_length.py`)** — keep reads inside `[--min_input_len, --max_input_len]`, then pad or trim them to `--target_length`. See the options table below. The script has no `--ligand` flag.
6. **One-window extract (`trim_to_kmer.py`)** — write **one** SEQ line per retained molecule (`--start` or `--center`, default **40 bp**). Short windows are dropped. This is not a tiled 40-mer library. A 154 bp filtered read is long enough for that 40-mer, so this step does not pad.

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

`run_preprocess.sh` accepts `--ligand` or `--target-length` and maps a ligand design onto `--target_length` for `filter_by_length.py`. If you omit the input window, the driver uses target±9 bp. That matches the filter script's own defaults for lig200 (target 154, window 145–163). The same ±9 bp tolerance is applied to the other designs; those windows were not measured separately. Set `--min-input-len` or `--max-input-len` to override.

| Ligand design | Target | Window used by the driver |
|---------------|--------|---------------------------|
| lig147        | 101    | 92–110                    |
| lig200        | 154    | 145–163                   |
| ligN40        | 187    | 178–196                   |
| ligN70        | 247    | 238–256                   |

## Generic naming convention

Replace experiment-specific IDs with a sample prefix `SAMPLE`. Expected inputs under `--input-dir`:

| Role | Filename |
|------|----------|
| Lane 1, read 1 | `{SAMPLE}_L1_R1.fastq.gz` |
| Lane 1, read 2 | `{SAMPLE}_L1_R2.fastq.gz` |
| Lane 2, read 1 | `{SAMPLE}_L2_R1.fastq.gz` |
| Lane 2, read 2 | `{SAMPLE}_L2_R2.fastq.gz` |

Outputs written under `--output-dir/{SAMPLE}/` (example with lig200 / 154 bp):

| Stage | Filename |
|-------|----------|
| Lane 1 merged | `{SAMPLE}_L1.fq.gz` |
| Lane 2 merged | `{SAMPLE}_L2.fq.gz` |
| Combined lanes | `{SAMPLE}.fq.gz`, `{SAMPLE}.fq` |
| FASTA | `{SAMPLE}.fa` |
| Length-filtered FASTA | `{SAMPLE}_{N}.fa` |
| One k-mer per read (plain text) | `{SAMPLE}_{N}.seq` |

`config/adapters.fa` is the default adapter FASTA. These are example adapters. If you use custom adapters, edit config/adapters.fa or pass your own file.

`data/examples/lig147_synthetic.fa` is a three-record synthetic FASTA for the length filter and k-mer extract only. It is not a sequencing run.

## Dependencies

- [`fastp`](https://github.com/opengene/fastp)
- [`seqtk`](https://github.com/lh3/seqtk)
- `gzip`, `bash`, `python3` (≥ 3.9 recommended)
- [Biopython](https://biopython.org/) (`pip install -r envs/requirements.txt`)

## Usage

```bash
scripts/01_preprocess/run_preprocess.sh \
  --sample SAMPLE \
  --ligand lig200 \
  --seed 1 \
  --input-dir ./raw \
  --output-dir ./processed
```

The command above uses `config/adapters.fa`. These are example adapters. If you use custom adapters, edit config/adapters.fa or pass your own file.

Provide exactly one of `--ligand` or `--target-length`. The driver has no numeric default. With the window omitted, lig200 is target 154 and window 145–163, which is what `filter_by_length.py` does on its own. `--seed` is forwarded to the length filter so its padding repeats. `--allow-pad` is only for a short k-mer window and still requires `--seed`. To extract a k-mer that is not the leftmost k bases, pass `--kmer-start` or `--center-kmer`. Pass `--adapter-fasta` to override the default adapters.

Equivalent standalone steps (mirroring the driver; generic names):

```bash
fastp -q 10 \
  --in1 SAMPLE_L1_R1.fastq.gz --in2 SAMPLE_L1_R2.fastq.gz \
  -m --merged_out SAMPLE_L1.fq.gz \
  --adapter_fasta=config/adapters.fa \
  --json=./SAMPLE_L1.json --html=./SAMPLE_L1.html \
  --overlap_len_require=5 -gx --length_required=15 \
  --n_base_limit=10 -y --complexity_threshold=5

fastp -q 10 \
  --in1 SAMPLE_L2_R1.fastq.gz --in2 SAMPLE_L2_R2.fastq.gz \
  -m --merged_out SAMPLE_L2.fq.gz \
  --adapter_fasta=config/adapters.fa \
  --json=./SAMPLE_L2.json --html=./SAMPLE_L2.html \
  --overlap_len_require=5 -gx --length_required=15 \
  --n_base_limit=10 -y --complexity_threshold=5

cat SAMPLE_L1.fq.gz SAMPLE_L2.fq.gz > SAMPLE.fq.gz
gzip -d -f -c SAMPLE.fq.gz > SAMPLE.fq
seqtk seq -a SAMPLE.fq > SAMPLE.fa
python3 scripts/01_preprocess/filter_by_length.py SAMPLE.fa SAMPLE_154.fa \
  --target_length 154 \
  --min_input_len 145 \
  --max_input_len 163 \
  --trim_from right \
  --seed 1
python3 scripts/01_preprocess/trim_to_kmer.py SAMPLE_154.fa SAMPLE_154.seq --kmer-size 40
```

The 154 bp FASTA is a valid input to `trim_to_kmer.py`. The default 40-mer is the first 40 bases of each read. Filter padding sits on the right, so that 40-mer is original sequence when the retained read was at least 40 bp.

### `filter_by_length.py` options

| Argument | Default | Role |
|----------|---------|------|
| `input_fasta` | required positional | Input FASTA |
| `output_fasta` | required positional | Output FASTA |
| `--target_length` | 154 | Length of every kept read |
| `--min_input_len` | 145 | Shortest input read to keep |
| `--max_input_len` | 163 | Longest input read to keep |
| `--trim_from` | `right` | `right`, `left`, or `center` |
| `--seed` | none | Seeds padding. Omit it for an unseeded draw |

Reads shorter than the minimum or longer than the maximum are dropped. A kept read shorter than the target is padded on the right with random A/T/C/G. A kept read longer than the target is trimmed from the side named by `--trim_from`. `center` drops half of the extra bases from the left (`extra // 2`) and keeps the next `target_length` bases.

```bash
python3 scripts/01_preprocess/filter_by_length.py \
  input.fa output_154.fa \
  --target_length 154 \
  --min_input_len 145 \
  --max_input_len 163 \
  --trim_from right \
  --seed 1
```

Stats printed by the script:

```text
[OK] Total reads:     N
[OK] Reads kept:     N
[OK] Reads padded:   N
[OK] Reads trimmed:  N
[OK] Reads dropped:  N
```

### Helper scripts

`trim_to_kmer.py` writes **one** SEQ line per record. `--length` is accepted as an alias for `--kmer-size`. There is no chunk/tile mode.

```bash
# One 40-base window per read (default: start=0). Shorts dropped, no padding.
python3 scripts/01_preprocess/trim_to_kmer.py input_154.fa output.seq --kmer-size 40

# Centered window
python3 scripts/01_preprocess/trim_to_kmer.py input_101.fa output.seq --kmer-size 40 --center

# Opt-in short-window padding (seed required)
python3 scripts/01_preprocess/trim_to_kmer.py input.fa output.seq --kmer-size 40 --allow-pad --seed 1
```

## `fastp` parameters (current defaults)

| Flag | Value | Role |
|------|-------|------|
| `-q` | 10 | Base quality threshold |
| `-m` / `--merged_out` | — | Merge overlapping PE reads |
| `--adapter_fasta` | `config/adapters.fa` | Example adapter sequences (override with your own file) |
| `--overlap_len_require` | 5 | Minimum overlap for merging |
| `-g` / `-x` | on | PolyG trimming / length filtering helpers as configured by fastp |
| `--length_required` | 15 | Discard very short merges |
| `--n_base_limit` | 10 | Max ambiguous bases |
| `-y` / `--complexity_threshold` | 5 | Low-complexity filter |

These defaults match the established SELEX preprocessing recipe used for this project and can be parameterised further as the pipeline evolves.
