#!/usr/bin/env bash
# Run a DIA-NN version on PlasmaBENCH SIM .d files with the UNIFIED PYE1 search
# settings (configs/dia-nn/diann-pye1.cfg). Both versions get the identical
# settings file — the only intended differences are the binary and its built-in
# defaults (QuantUMS, DL models): that IS the version comparison.
#
# Pass Sample A and Sample B .d together so DIA-NN's cross-run / MBR step
# (--reanalyse, in the cfg) can build the A-vs-B quant matrix.
#
# Usage:
#   bash configs/dia-nn/run.sh <1.8|2.5> <fasta> <out_dir> <runA.d> <runB.d> [...]
# Env:
#   THREADS   (default: nproc)
#
# Example (once the SIM .d's and FASTA exist):
#   bash configs/dia-nn/run.sh 1.8 data/raw/fasta/plasma_yeast_ecoli.fasta \
#     results/diann-1.8 \
#     simulations/stage1/dia/sampleA_imposed/*.d \
#     simulations/stage1/dia/sampleB_imposed/*.d
set -euo pipefail

REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)

[[ $# -ge 4 ]] || { echo "usage: run.sh <1.8|2.5> <fasta> <out_dir> <d1> [d2 ...]" >&2; exit 1; }
VER="$1"; FASTA="$2"; OUT_DIR="$3"; shift 3

case "$VER" in
  1.8|1.8.1) BIN="$REPO_ROOT/tools/diann/1.8.1/diann-1.8.1" ;;
  2.5|2.5.0) BIN="$REPO_ROOT/tools/diann/2.5.0/diann-2.5.0/diann-linux" ;;
  2.6|2.6.0) BIN="$REPO_ROOT/tools/diann/2.6.0/diann-2.6.0/diann-linux" ;;
  *) echo "unknown DIA-NN version '$VER' (use 1.8, 2.5 or 2.6)" >&2; exit 1 ;;
esac
[[ -x "$BIN" ]] || { echo "missing/!exec $BIN — run: bash scripts/get_diann.sh" >&2; exit 2; }
[[ -e "$FASTA" ]] || { echo "FASTA not found: $FASTA (TODO: plasma+yeast+ecoli reference)" >&2; exit 2; }

CFG="$REPO_ROOT/configs/dia-nn/diann-pye1.cfg"
[[ -f "$CFG" ]] || { echo "missing settings file: $CFG" >&2; exit 2; }
# Tokenize the cfg into an array. `set -f` disables pathname expansion so values
# like `--cut K*,R*` are NOT glob-expanded against the cwd (an unquoted $SETTINGS
# would corrupt the command if any file matched K*/R*).
set -f
SETTINGS=( $(grep -vE '^\s*(#|$)' "$CFG") )
set +f

# NO_MBR=1 strips --reanalyse (DIA-NN's match-between-runs / double-pass) — control
# to test whether low-abundance A/B ratio behaviour is MBR transfer vs intrinsic.
if [[ "${NO_MBR:-0}" == "1" ]]; then
  FILTERED=(); for s in "${SETTINGS[@]}"; do [[ "$s" == "--reanalyse" ]] || FILTERED+=("$s"); done
  SETTINGS=( "${FILTERED[@]}" )
  echo "==> NO_MBR set: --reanalyse stripped"
fi

THREADS="${THREADS:-$(nproc)}"
mkdir -p "$OUT_DIR"

# Validate every .d input up front; MBR (--reanalyse) needs >=2.
F_ARGS=()
for d in "$@"; do
  [[ -e "$d" ]] || { echo "input .d not found: $d" >&2; exit 2; }
  F_ARGS+=(--f "$d")
done
[[ $# -ge 2 ]] || echo "WARNING: only $# run(s) given; MBR (--reanalyse) needs >=2 (pass A and B together)." >&2

echo "==> DIA-NN $VER  ($BIN)"
echo "    fasta=$FASTA  out=$OUT_DIR  threads=$THREADS  runs=$#"
set -x
"$BIN" "${SETTINGS[@]}" \
  --fasta "$FASTA" \
  "${F_ARGS[@]}" \
  --out "$OUT_DIR/report.tsv" \
  --temp "$OUT_DIR" \
  --threads "$THREADS"
set +x

# DIA-NN 1.8 honours `--out report.tsv`; 2.x ALWAYS writes the main report as
# report.parquet regardless of the requested extension. Accept whichever exists.
REPORT=""
for cand in "$OUT_DIR/report.parquet" "$OUT_DIR/report.tsv"; do
  if [[ -s "$cand" ]]; then REPORT="$cand"; break; fi
done
if [[ -z "$REPORT" ]]; then
  echo "ERROR: DIA-NN produced no (or empty) report (report.parquet/.tsv) in $OUT_DIR" >&2
  exit 4
fi
echo "==> done: $REPORT"
