# Motif processing

This module runs motif discovery on SEQ files produced by [`preprocessing/trim_to_kmer.py`](../preprocessing/trim_to_kmer.py) (one fixed-width line per molecule), using **`totalautomodelIC`** (v0.22) with **`spacek40`**.

`run_motif.sh` **refuses to run** if the background or signal SEQ is empty, or if any line length ≠ `--seq-length` (default **40**). That guard exists so tiled/padded/variable-width libraries cannot silently feed `totalautomodelIC`.

## What it does

Given a background SEQ library and a signal (enriched) SEQ library:

1. Find local-maximum k-mer seeds (`spacek40 -local`).
2. Build hit logos / PFMs for each seed (`spacek40 --f`).
3. Filter by complexity and pairwise logo similarity.
4. Refine accepted seeds (information content between 10 and 20 bits).
5. Write accepted / unrefined / rejected models into an SVG logo file.

The shipped `totalautomodelIC` script is the established analysis recipe and is **not modified** in this repository. It must be run with `spacek40` available as `./spacek40` in the working directory (the driver handles this).

## Generic naming

Original experiment-specific commands:

```bash
./totalautomodelIC -40N HTS147xlSOX3DBDc1_101.seq HTS147xlSOX3DBDc4_101.seq \
  1 6 6 0.25 - 20 50 0.1 HTS147xlSOX3DBDc4v1

./totalautomodelIC -40N HTS147xlSOX3DBDc1_101.seq HTS147xlSOX3DBDc3_101.seq \
  1 6 6 0.25 - 20 50 0.1 HTS147xlSOX3DBDc3v1
```

Generic equivalents:

| Role | Original | Generic |
|------|----------|---------|
| Background SEQ | `HTS147xlSOX3DBDc1_101.seq` | `BACKGROUND_101.seq` |
| Signal SEQ (round 4) | `HTS147xlSOX3DBDc4_101.seq` | `SIGNAL_r4_101.seq` |
| Signal SEQ (round 3) | `HTS147xlSOX3DBDc3_101.seq` | `SIGNAL_r3_101.seq` |
| Output prefix (r4) | `HTS147xlSOX3DBDc4v1` | `SIGNAL_r4v1` |
| Output prefix (r3) | `HTS147xlSOX3DBDc3v1` | `SIGNAL_r3v1` |

```bash
./totalautomodelIC -40N BACKGROUND_101.seq SIGNAL_r4_101.seq \
  1 6 6 0.25 - 20 50 0.1 SIGNAL_r4v1

./totalautomodelIC -40N BACKGROUND_101.seq SIGNAL_r3_101.seq \
  1 6 6 0.25 - 20 50 0.1 SIGNAL_r3v1
```

SEQ inputs are one-sequence-per-line k-mers from preprocessing (`{SAMPLE}_{N}.seq`), one extract per molecule. The `-40N` flag must match that width; the driver checks every line before invoking `totalautomodelIC`.

## `totalautomodelIC` argument order

```text
./totalautomodelIC -[seq length] background_seq signal_seq multinomial \
  min_seed_length max_seed_length length_cutoff iupac_cutoff \
  local_max_cutoff_count logo_number_cutoff logo_similarity_cutoff \
  output_prefix
```

Defaults used by this pipeline (matching the original recipe):

| Argument | Default |
|----------|---------|
| seq length flag | `-40N` |
| multinomial | `1` |
| min / max seed length | `6` / `6` |
| length cutoff | `0.25` |
| IUPAC cutoff | `-` (disabled) |
| local-max cutoff count | `20` |
| logo number cutoff | `50` |
| logo similarity cutoff | `0.1` |

## Dependencies

- `totalautomodelIC` — included in this directory (bash; do not edit)
- `spacek40` — external binary invoked as `./spacek40` by `totalautomodelIC`; place it in `motif/` (or pass `--bin-dir`)
- `bash`, `bc`, `awk`, `sed`, `grep`, `coreutils`

## Usage

Place `spacek40` next to `totalautomodelIC`, then:

```bash
./motif/run_motif.sh \
  --background ./processed/BACKGROUND/BACKGROUND_101.seq \
  --signal ./processed/SIGNAL_r4/SIGNAL_r4_101.seq \
  --output-prefix SIGNAL_r4v1 \
  --work-dir ./motif_out/SIGNAL_r4
```

```bash
./motif/run_motif.sh \
  --background ./processed/BACKGROUND/BACKGROUND_101.seq \
  --signal ./processed/SIGNAL_r3/SIGNAL_r3_101.seq \
  --output-prefix SIGNAL_r3v1 \
  --work-dir ./motif_out/SIGNAL_r3
```

Primary output: `{output_prefix}_logos.svg` in `--work-dir`, plus refined PFMs / seed lists written by `totalautomodelIC`.

## Files

| File | Role |
|------|------|
| `totalautomodelIC` | Upstream motif automation script (unchanged) |
| `run_motif.sh` | Generic-named driver around `totalautomodelIC`; refuses SEQ lines whose length ≠ `--seq-length` |
| `spacek40` | Required external binary (not shipped; supply locally) |
