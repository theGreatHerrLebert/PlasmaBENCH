#!/usr/bin/env bash
# Isolate the 0-ID culprit AND produce the real A+B result.
#
# Confirmed so far:
#   v2  (0 IDs):      prebuilt --lib + --reanalyse + --mass-acc 0/0  [A,B]
#   diag(26441 IDs):  prebuilt --lib + (no reanalyse) + --mass-acc 15/15  [A only]
# This run reuses the prebuilt lib but keeps --reanalyse (MBR, needed for the
# A/B quant matrix) WITH the fixed 15 ppm mass accuracy. If it returns IDs,
# the culprit was auto mass-acc (--mass-acc 0); if 0, --reanalyse is implicated.
#
# Usage: bash scripts/diag_diann_AB_mbr.sh
set -euo pipefail
REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO_ROOT"

BIN="tools/diann/2.5.0/diann-2.5.0/diann-linux"
LIB="results/diann-2.5-AB/report-lib.predicted.speclib"
DA="simulations/stage1/dia/sampleA_imposed/PLB-S1-DIA-sampleA-imposed/PLB-S1-DIA-sampleA-imposed.d"
DB="simulations/stage1/dia/sampleB_imposed/PLB-S1-DIA-sampleB-imposed/PLB-S1-DIA-sampleB-imposed.d"
OUT_DIR="results/diag-diann-2.5-AB-mbr-15ppm"
THREADS="${THREADS:-14}"

mkdir -p "$OUT_DIR"
echo "==> A+B, MBR on, mass-acc fixed 15/15 (reusing prebuilt lib)"
"$BIN" \
  --lib "$LIB" \
  --fasta "data/raw/fasta/plasma_yeast_ecoli.fasta" \
  --f "$DA" --f "$DB" \
  --out "$OUT_DIR/report.tsv" \
  --temp "$OUT_DIR" \
  --threads "$THREADS" \
  --qvalue 0.01 \
  --matrix-qvalue 0.01 \
  --reanalyse \
  --mass-acc 15 \
  --mass-acc-ms1 15

echo "==> done: $OUT_DIR/report.tsv"
