# Motif processing

Paths below are relative to the repository root. Scripts live in
`scripts/02_motif/`.

This module runs motif discovery on SEQ files produced by [`scripts/01_preprocess/trim_to_kmer.py`](../scripts/01_preprocess/trim_to_kmer.py) (one fixed-width line per molecule), using **`totalautomodelIC`** (v0.22) with **`spacek40`**.

`run_motif.sh` **refuses to run** if the background or signal SEQ is empty, or if any line length ≠ `--seq-length` (default **40**). That guard exists so tiled/padded/variable-width libraries cannot silently feed `totalautomodelIC`.

## What it does

Given a background SEQ library and a signal (enriched) SEQ library:

1. Find local-maximum k-mer seeds (`spacek40 -local`).
2. Build hit logos / PFMs for each seed (`spacek40 --f`).
3. Filter by complexity and pairwise logo similarity.
4. Refine accepted seeds (information content between 10 and 20 bits).
5. Write accepted / unrefined / rejected models into an SVG logo file.

The shipped `totalautomodelIC` script (`scripts/02_motif/totalautomodelIC`) is the established analysis recipe and is **not modified** in this repository. It calls `./spacek40` in the working directory. `run_motif.sh` symlinks the shipped binary `scripts/02_motif/spacek40` there. `--bin-dir` defaults to that directory.

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

- `totalautomodelIC` — included at `scripts/02_motif/totalautomodelIC` (bash; do not edit)
- `spacek40` — included at `scripts/02_motif/spacek40`. `totalautomodelIC` invokes it as `./spacek40`. The driver symlinks this copy into `--work-dir`. Pass `--bin-dir` only to use a different copy.
- `bash`, `bc`, `awk`, `sed`, `grep`, `coreutils`

## Usage

`spacek40` is already next to `totalautomodelIC`. From the repository root:

```bash
scripts/02_motif/run_motif.sh \
  --background ./processed/BACKGROUND/BACKGROUND_101.seq \
  --signal ./processed/SIGNAL_r4/SIGNAL_r4_101.seq \
  --output-prefix SIGNAL_r4v1 \
  --work-dir ./motif_out/SIGNAL_r4
```

```bash
scripts/02_motif/run_motif.sh \
  --background ./processed/BACKGROUND/BACKGROUND_101.seq \
  --signal ./processed/SIGNAL_r3/SIGNAL_r3_101.seq \
  --output-prefix SIGNAL_r3v1 \
  --work-dir ./motif_out/SIGNAL_r3
```

Primary output: `{output_prefix}_logos.svg` in `--work-dir`, plus refined PFMs / seed lists written by `totalautomodelIC`.

## Running steps standalone

These are the two steps inside `run_motif.sh`, in order. `spacek40` is `scripts/02_motif/spacek40`. The values below are that script's defaults.

The driver refuses a SEQ file before it looks for the motif binaries. This is `validate_seq_width`, called first on the background file and then on the signal file:

```bash
BACKGROUND=processed/BACKGROUND/BACKGROUND_101.seq
SIGNAL=processed/SIGNAL_r4/SIGNAL_r4_101.seq
SEQ_LENGTH=40

validate_seq_width() {
  local file="$1"
  local expected="$2"
  local label="$3"
  python3 - "$file" "$expected" "$label" <<'PY'
import sys

path, expected_s, label = sys.argv[1], sys.argv[2], sys.argv[3]
try:
    expected = int(expected_s)
except ValueError:
    sys.exit(f"[motif] error: --seq-length must be an integer, got {expected_s!r}")
if expected <= 0:
    sys.exit(f"[motif] error: --seq-length must be positive, got {expected}")

n = 0
with open(path, encoding="utf-8") as handle:
    for lineno, raw in enumerate(handle, 1):
        seq = raw.rstrip("\n\r")
        if not seq:
            sys.exit(f"[motif] error: {label} SEQ line {lineno} is empty ({path})")
        if len(seq) != expected:
            sys.exit(
                f"[motif] error: {label} SEQ line {lineno} length {len(seq)} "
                f"!= expected {expected} ({path})"
            )
        n += 1
if n == 0:
    sys.exit(f"[motif] error: {label} SEQ is empty ({path})")
PY
}

validate_seq_width "$BACKGROUND" "$SEQ_LENGTH" "background"
validate_seq_width "$SIGNAL" "$SEQ_LENGTH" "signal"
```

A line counts after the newline is stripped. The check exits if `--seq-length` is not a positive integer, if the file has no sequence lines, or if any line is empty or the wrong length.

`totalautomodelIC` calls `./spacek40` and writes relative to the current directory. The driver turns the SEQ paths and `BIN_DIR` into absolute paths, links `totalautomodelIC` and `spacek40` into the work directory, then runs this command. `LENGTH_FLAG` is one shell argument. `--all-gaps` changes that argument from `-${SEQ_LENGTH}N` to the single string `-${SEQ_LENGTH}N -allgaps`. If `OUTPUT_PREFIX` is empty, the script uses the signal filename without `.seq`, plus `v1`.

```bash
BIN_DIR=scripts/02_motif
WORK_DIR=motif_out/SIGNAL_r4
OUTPUT_PREFIX=SIGNAL_r4v1
MULTINOMIAL=1
MIN_SEED=6
MAX_SEED=6
LENGTH_CUTOFF=0.25
IUPAC_CUTOFF="-"
LOCAL_MAX_CUTOFF=20
LOGO_NUMBER=50
LOGO_SIMILARITY=0.1

BACKGROUND="$(cd "$(dirname "$BACKGROUND")" && pwd)/$(basename "$BACKGROUND")"
SIGNAL="$(cd "$(dirname "$SIGNAL")" && pwd)/$(basename "$SIGNAL")"
BIN_DIR="$(cd "$BIN_DIR" && pwd)"
mkdir -p "$WORK_DIR"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"

TOTALAUTOMODELIC="${BIN_DIR}/totalautomodelIC"
SPACEK40="${BIN_DIR}/spacek40"
chmod +x "$TOTALAUTOMODELIC" "$SPACEK40" 2>/dev/null || true
LENGTH_FLAG="-${SEQ_LENGTH}N"

cd "$WORK_DIR"
ln -sfn "$TOTALAUTOMODELIC" ./totalautomodelIC
ln -sfn "$SPACEK40" ./spacek40

./totalautomodelIC \
  "$LENGTH_FLAG" \
  "$BACKGROUND" \
  "$SIGNAL" \
  "$MULTINOMIAL" \
  "$MIN_SEED" \
  "$MAX_SEED" \
  "$LENGTH_CUTOFF" \
  "$IUPAC_CUTOFF" \
  "$LOCAL_MAX_CUTOFF" \
  "$LOGO_NUMBER" \
  "$LOGO_SIMILARITY" \
  "$OUTPUT_PREFIX"
```

| Argument | Default | Role |
|----------|---------|------|
| `$LENGTH_FLAG` | `-40N` | Sequence length passed as one argument (`-${SEQ_LENGTH}N`). `--all-gaps` makes it `-${SEQ_LENGTH}N -allgaps`, still one argument |
| `$BACKGROUND` | — | Background SEQ path |
| `$SIGNAL` | — | Enriched SEQ path |
| `$MULTINOMIAL` | `1` | Multinomial setting |
| `$MIN_SEED` | `6` | Minimum seed length |
| `$MAX_SEED` | `6` | Maximum seed length |
| `$LENGTH_CUTOFF` | `0.25` | Length cutoff |
| `$IUPAC_CUTOFF` | `-` | IUPAC cutoff; `-` disables it |
| `$LOCAL_MAX_CUTOFF` | `20` | Minimum count for a local-maximum k-mer |
| `$LOGO_NUMBER` | `50` | Maximum number of logos |
| `$LOGO_SIMILARITY` | `0.1` | Logo similarity cutoff |
| `$OUTPUT_PREFIX` | `<signal_basename>v1` | Prefix for `${OUTPUT_PREFIX}_logos.svg` and the other motif files |

## Files

| File | Role |
|------|------|
| `scripts/02_motif/totalautomodelIC` | Upstream motif automation script (unchanged) |
| `scripts/02_motif/run_motif.sh` | Driver around `totalautomodelIC`; refuses SEQ lines whose length ≠ `--seq-length`. `--bin-dir` defaults to this directory |
| `scripts/02_motif/spacek40` | Background binary used by `totalautomodelIC` |
