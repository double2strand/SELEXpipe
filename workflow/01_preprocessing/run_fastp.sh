#!/usr/bin/env bash
# One lane of the fastp merge used by run_preprocess.sh.
# Keep the fastp invocation identical to run_fastp_lane in that driver.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DEFAULT_ADAPTER_FASTA="${REPO_ROOT}/config/adapters.fa"

usage() {
  cat <<'EOF'
Usage: run_fastp.sh in1 in2 out_prefix [adapters] [threads]

Merge one paired-end lane with the same fastp flags as run_preprocess.sh.

  in1         Read 1 FASTQ (typically .fastq.gz)
  in2         Read 2 FASTQ
  out_prefix  Output path prefix:
                ${out_prefix}.fq.gz   merged reads
                ${out_prefix}.json    fastp JSON report
                ${out_prefix}.html    fastp HTML report
  adapters    Adapter FASTA (default: config/adapters.fa)
  threads     Value for fastp --thread (default: nproc, or 4)

These are example adapters. If you use custom adapters, edit config/adapters.fa or pass your own file.

Run once per lane. Example for lane 1, from the repository root:

  scripts/01_preprocess/run_fastp.sh \
      raw/SAMPLE_L1_R1.fastq.gz \
      raw/SAMPLE_L1_R2.fastq.gz \
      processed/SAMPLE/SAMPLE_L1

  scripts/01_preprocess/run_fastp.sh \
      raw/SAMPLE_L2_R1.fastq.gz \
      raw/SAMPLE_L2_R2.fastq.gz \
      processed/SAMPLE/SAMPLE_L2
EOF
}

die() {
  printf '[fastp] error: %s\n' "$*" >&2
  exit 1
}

if [[ $# -eq 0 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 3 || $# -gt 5 ]]; then
  usage >&2
  die "expected in1 in2 out_prefix [adapters] [threads]"
fi

in1="$1"
in2="$2"
out_prefix="$3"
ADAPTER_FASTA="${4:-$DEFAULT_ADAPTER_FASTA}"
if [[ $# -ge 5 ]]; then
  THREADS="$5"
else
  THREADS="$(nproc 2>/dev/null || echo 4)"
fi

merged="${out_prefix}.fq.gz"
json="${out_prefix}.json"
html="${out_prefix}.html"

[[ -f "$in1" ]] || die "missing in1: $in1"
[[ -f "$in2" ]] || die "missing in2: $in2"
[[ -f "$ADAPTER_FASTA" ]] || die "adapter FASTA not found: $ADAPTER_FASTA"
[[ "$THREADS" =~ ^[0-9]+$ ]] || die "threads must be a positive integer, got: $THREADS"
command -v fastp >/dev/null 2>&1 || die "required command not found: fastp"

mkdir -p "$(dirname "$out_prefix")"

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
