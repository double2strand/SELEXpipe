#!/usr/bin/env bash
# SELEXpipe preprocessing driver
#
# For each sample:
#   1. fastp: quality/adapter trim + merge paired-end reads (per lane)
#   2. concatenate lane-merged FASTQ.gz files
#   3. decompress (keep .fq.gz) and convert FASTQ -> FASTA (seqtk)
#   4. filter FASTA to exact merged-read length (drop shorts; no padding)
#   5. extract one fixed-width window per molecule (default 40-mer)
#
# Default path never writes random ATCG. Padding is opt-in (--allow-pad + --seed).
#
# Expected per-sample lane inputs (under --input-dir):
#   ${SAMPLE}_L1_R1.fastq.gz  ${SAMPLE}_L1_R2.fastq.gz
#   ${SAMPLE}_L2_R1.fastq.gz  ${SAMPLE}_L2_R2.fastq.gz
#
# Outputs (under --output-dir / ${SAMPLE}/):
#   ${SAMPLE}_L1.fq.gz   ${SAMPLE}_L2.fq.gz
#   ${SAMPLE}.fq.gz      ${SAMPLE}.fq
#   ${SAMPLE}.fa
#   ${SAMPLE}_${TARGET_LEN}.fa
#   ${SAMPLE}_${TARGET_LEN}.seq
#
# Ligand design expected merged lengths (do not "correct"):
#   lig147  101 bp
#   lig200  154 bp
#   ligN40  187 bp
#   ligN70  247 bp
# ligN40 vs ligN70 constant flanks differ; the 60 bp delta is expected.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SAMPLE=""
INPUT_DIR="."
OUTPUT_DIR="."
ADAPTER_FASTA=""
TARGET_LEN=""
LIGAND=""
MIN_INPUT_LEN=""
MAX_INPUT_LEN=""
TRIM_FROM="right"
SEED=""
ALLOW_PAD=0
KMER_SIZE=40
KMER_START=0
CENTER_KMER=0
THREADS="$(nproc 2>/dev/null || echo 4)"
SKIP_FASTP=0

# Expected merged-read lengths for supported ligand designs.
ligand_length() {
  case "$1" in
    lig147) echo 101 ;;
    lig200) echo 154 ;;
    ligN40) echo 187 ;;
    ligN70) echo 247 ;;
    *) return 1 ;;
  esac
}

usage() {
  cat <<'EOF'
Usage: run_preprocess.sh --sample SAMPLE --adapter-fasta ADAPTERS.fa \
                         (--target-length N | --ligand DESIGN) [options]

Required:
  --sample NAME              Sample prefix (generic; replaces experiment-specific IDs)
  --adapter-fasta PATH       Adapter FASTA for fastp (e.g. Adapters_TruSeq3PE.fa)
  --target-length N          Exact merged-read length for filter_by_length.py
  --ligand DESIGN            Ligand design key (sets target length; see table below)

  Provide exactly one of --target-length or --ligand.

Ligand design expected merged lengths (do not "correct"):
  lig147  101 bp
  lig200  154 bp
  ligN40  187 bp
  ligN70  247 bp

  ligN40 vs ligN70 constant flanks differ; the 60 bp (not 30 bp) delta is expected.

Optional:
  --input-dir DIR            Directory containing lane FASTQ.gz files (default: .)
  --output-dir DIR           Directory for outputs (default: .); writes into DIR/SAMPLE/
  --min-input-len N          Min input length to retain (default: same as target length)
  --max-input-len N          Max input length to retain (default: same as target length)
  --trim-from MODE           right|left|center trim for long reads (default: right)
  --kmer-size K              One SEQ line / k-mer per molecule (default: 40)
  --kmer-start S             0-based start of the extract window (default: 0)
  --center-kmer              Center the extract window instead of --kmer-start
  --allow-pad                Opt-in: pad short reads/windows with random ATCG
  --seed N                   RNG seed (required if --allow-pad is set)
  --threads N                Threads hint for fastp (default: nproc)
  --skip-fastp               Resume from existing lane-merged FASTQ.gz outputs
  -h, --help                 Show this help

Default path: exact-length filter, one window per read, drop shorts, no padding.

Example:
  ./run_preprocess.sh \
      --sample SAMPLE \
      --adapter-fasta ./Adapters_TruSeq3PE.fa \
      --ligand lig147 \
      --input-dir ./raw \
      --output-dir ./processed
EOF
}

log() {
  printf '[preprocess] %s\n' "$*"
}

die() {
  printf '[preprocess] error: %s\n' "$*" >&2
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="${2:-}"; shift 2 ;;
    --adapter-fasta) ADAPTER_FASTA="${2:-}"; shift 2 ;;
    --input-dir) INPUT_DIR="${2:-}"; shift 2 ;;
    --output-dir) OUTPUT_DIR="${2:-}"; shift 2 ;;
    --target-length|--expected-length) TARGET_LEN="${2:-}"; shift 2 ;;
    --ligand) LIGAND="${2:-}"; shift 2 ;;
    --min-input-len) MIN_INPUT_LEN="${2:-}"; shift 2 ;;
    --max-input-len) MAX_INPUT_LEN="${2:-}"; shift 2 ;;
    --trim-from) TRIM_FROM="${2:-}"; shift 2 ;;
    --seed) SEED="${2:-}"; shift 2 ;;
    --allow-pad) ALLOW_PAD=1; shift ;;
    --kmer-size) KMER_SIZE="${2:-}"; shift 2 ;;
    --kmer-start) KMER_START="${2:-}"; shift 2 ;;
    --center-kmer) CENTER_KMER=1; shift ;;
    --threads) THREADS="${2:-}"; shift 2 ;;
    --skip-fastp) SKIP_FASTP=1; shift ;;
    --kmer-seed)
      die "--kmer-seed is removed. Default extract does not pad; use --allow-pad --seed"
      ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ -n "$SAMPLE" ]] || { usage; die "--sample is required"; }
[[ -n "$ADAPTER_FASTA" ]] || { usage; die "--adapter-fasta is required"; }
[[ -f "$ADAPTER_FASTA" ]] || die "adapter FASTA not found: $ADAPTER_FASTA"
[[ -d "$INPUT_DIR" ]] || die "input directory not found: $INPUT_DIR"

if [[ -n "$LIGAND" && -n "$TARGET_LEN" ]]; then
  die "provide only one of --ligand or --target-length"
fi
if [[ -n "$LIGAND" ]]; then
  TARGET_LEN="$(ligand_length "$LIGAND")" || die "unknown ligand design: $LIGAND (expected lig147, lig200, ligN40, ligN70)"
elif [[ -z "$TARGET_LEN" ]]; then
  usage
  die "provide --target-length N or --ligand DESIGN"
fi

case "$TRIM_FROM" in
  right|left|center) ;;
  *) die "--trim-from must be one of: right, left, center" ;;
esac

if [[ "$ALLOW_PAD" -eq 1 && -z "$SEED" ]]; then
  die "--allow-pad requires --seed so padding is reproducible"
fi

require_cmd fastp
require_cmd gzip
require_cmd seqtk
require_cmd python3
require_cmd cat

OUT="${OUTPUT_DIR%/}/${SAMPLE}"
mkdir -p "$OUT"

L1_R1="${INPUT_DIR%/}/${SAMPLE}_L1_R1.fastq.gz"
L1_R2="${INPUT_DIR%/}/${SAMPLE}_L1_R2.fastq.gz"
L2_R1="${INPUT_DIR%/}/${SAMPLE}_L2_R1.fastq.gz"
L2_R2="${INPUT_DIR%/}/${SAMPLE}_L2_R2.fastq.gz"

for f in "$L1_R1" "$L1_R2" "$L2_R1" "$L2_R2"; do
  [[ -f "$f" ]] || die "missing input FASTQ: $f"
done

L1_MERGED="${OUT}/${SAMPLE}_L1.fq.gz"
L2_MERGED="${OUT}/${SAMPLE}_L2.fq.gz"
COMBINED_GZ="${OUT}/${SAMPLE}.fq.gz"
COMBINED_FQ="${OUT}/${SAMPLE}.fq"
FASTA="${OUT}/${SAMPLE}.fa"
FASTA_FILTERED="${OUT}/${SAMPLE}_${TARGET_LEN}.fa"
KMER_OUT="${OUT}/${SAMPLE}_${TARGET_LEN}.seq"

run_fastp_lane() {
  local in1="$1" in2="$2" merged="$3" lane_tag="$4"
  local json="${OUT}/${SAMPLE}_${lane_tag}.json"
  local html="${OUT}/${SAMPLE}_${lane_tag}.html"

  log "fastp merge ${lane_tag}: $(basename "$in1") + $(basename "$in2")"
  fastp \
    -q 10 \
    --in1 "$in1" \
    --in2 "$in2" \
    -m \
    --merged_out "$merged" \
    --adapter_fasta="$ADAPTER_FASTA" \
    --json="$json" \
    --html="$html" \
    --overlap_len_require=5 \
    -gx \
    --length_required=15 \
    --n_base_limit=10 \
    -y \
    --complexity_threshold=5 \
    --thread "$THREADS"
}

if [[ "$SKIP_FASTP" -eq 0 ]]; then
  run_fastp_lane "$L1_R1" "$L1_R2" "$L1_MERGED" "L1"
  run_fastp_lane "$L2_R1" "$L2_R2" "$L2_MERGED" "L2"
else
  [[ -f "$L1_MERGED" && -f "$L2_MERGED" ]] || die "--skip-fastp requires existing lane merges"
  log "skipping fastp; using existing lane merges"
fi

log "concatenating lane-merged FASTQ.gz -> $(basename "$COMBINED_GZ")"
cat "$L1_MERGED" "$L2_MERGED" > "$COMBINED_GZ"

# Keep the gzipped combined file; write an uncompressed copy for seqtk.
log "decompressing (keep .fq.gz) -> $(basename "$COMBINED_FQ")"
gzip -d -f -c "$COMBINED_GZ" > "$COMBINED_FQ"

log "FASTQ -> FASTA (seqtk) -> $(basename "$FASTA")"
seqtk seq -a "$COMBINED_FQ" > "$FASTA"

if [[ -n "$LIGAND" ]]; then
  log "filtering FASTA for ligand ${LIGAND} -> exact length ${TARGET_LEN} (filter_by_length.py)"
else
  log "filtering FASTA to exact length ${TARGET_LEN} (filter_by_length.py)"
fi
FILTER_ARGS=(
  "$FASTA"
  "$FASTA_FILTERED"
  --trim_from "$TRIM_FROM"
)
if [[ -n "$LIGAND" ]]; then
  FILTER_ARGS+=(--ligand "$LIGAND")
else
  FILTER_ARGS+=(--target_length "$TARGET_LEN")
fi
[[ -n "$MIN_INPUT_LEN" ]] && FILTER_ARGS+=(--min_input_len "$MIN_INPUT_LEN")
[[ -n "$MAX_INPUT_LEN" ]] && FILTER_ARGS+=(--max_input_len "$MAX_INPUT_LEN")
if [[ "$ALLOW_PAD" -eq 1 ]]; then
  FILTER_ARGS+=(--allow-pad --seed "$SEED")
fi
python3 "${SCRIPT_DIR}/filter_by_length.py" "${FILTER_ARGS[@]}"

if [[ "$CENTER_KMER" -eq 1 ]]; then
  log "extracting one centered ${KMER_SIZE}-mer per read -> $(basename "$KMER_OUT")"
else
  log "extracting one ${KMER_SIZE}-mer per read (start=${KMER_START}) -> $(basename "$KMER_OUT")"
fi
TRIM_ARGS=("$FASTA_FILTERED" "$KMER_OUT" --kmer-size "$KMER_SIZE")
if [[ "$CENTER_KMER" -eq 1 ]]; then
  TRIM_ARGS+=(--center)
else
  TRIM_ARGS+=(--start "$KMER_START")
fi
if [[ "$ALLOW_PAD" -eq 1 ]]; then
  TRIM_ARGS+=(--allow-pad --seed "$SEED")
fi
python3 "${SCRIPT_DIR}/trim_to_kmer.py" "${TRIM_ARGS[@]}"

log "done"
log "  lane merges : $L1_MERGED , $L2_MERGED"
log "  combined    : $COMBINED_GZ / $COMBINED_FQ"
log "  fasta       : $FASTA"
log "  length filt : $FASTA_FILTERED"
log "  k-mers      : $KMER_OUT"
