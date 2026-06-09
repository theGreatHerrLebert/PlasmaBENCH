#!/usr/bin/env bash
# Diagnostic: DIA-NN 2.5 the way 2.x actually wants it — TWO STEPS.
#
# The one-shot run (configs/dia-nn/run.sh) triggered DIA-NN 2.5's warning
# "the in silico-predicted library must be generated in a separate pipeline
# step and then used to process the raw data" and returned 0 IDs. In DIA-NN
# 2.x, library-free must be split: (1) predict the library from FASTA with no
# raw files, (2) search the raw with --lib and NO --fasta-search.
#
# Step 1 was already done by the previous run, so we reuse its predicted
# library and only run step 2 here (the cheap part). We also FIX the mass
# accuracies to 15/15 ppm (DIA-NN's own timsTOF recommendation) instead of
# auto, and drop --reanalyse (single-file diagnostic — MBR needs >=2 and isn't
# what we're testing).
#
# Usage: bash scripts/diag_diann_twostep.sh
set -euo pipefail
REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO_ROOT"

BIN="tools/diann/2.5.0/diann-2.5.0/diann-linux"
LIB="results/diann-2.5-AB/report-lib.predicted.speclib"   # reuse step-1 output
DCONV="simulations/stage1/dia/sampleA_imposed/PLB-S1-DIA-sampleA-imposed/PLB-S1-DIA-sampleA-imposed.d"
OUT_DIR="results/diag-diann-2.5-twostep-A"
THREADS="${THREADS:-14}"

mkdir -p "$OUT_DIR"
echo "==> Step 2 only (reusing predicted lib): search sampleA against --lib, fixed 15 ppm"
"$BIN" \
  --lib "$LIB" \
  --fasta "data/raw/fasta/plasma_yeast_ecoli.fasta" \
  --f "$DCONV" \
  --out "$OUT_DIR/report.tsv" \
  --temp "$OUT_DIR" \
  --threads "$THREADS" \
  --qvalue 0.01 \
  --matrix-qvalue 0.01 \
  --mass-acc 15 \
  --mass-acc-ms1 15 \
  --min-pep-len 7 --max-pep-len 30 \
  --min-pr-charge 1 --max-pr-charge 4 \
  --min-pr-mz 300 --max-pr-mz 1800 \
  --min-fr-mz 200 --max-fr-mz 1800 \
  --unimod4 --var-mods 0

echo "==> done: $OUT_DIR/report.tsv"
