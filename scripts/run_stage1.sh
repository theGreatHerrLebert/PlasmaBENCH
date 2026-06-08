#!/usr/bin/env bash
# Stage 1 orchestration: build seeds (if reports given) → render configs →
# run timsim for PYE1 Sample A and Sample B.
#
# Run from a from-source rustims venv until the published Docker image contains
# from_findings (see docs/venv-setup.md), or inside `make shell`.
#
# Pre-reqs:
#   - simulations/stage1/seeds/seed_sample{A,B}.csv exist, OR pass
#     REPORT_A / REPORT_B (+ ENGINE) to build them via build_seed_from_report.py.
#   - rustims submodule on feature/simulate-from-findings OR `make venv-source` done.
#
# Usage:
#   bash scripts/run_stage1.sh                                  # seeds already built
#   ENGINE=diann REPORT_A=.../A/report.tsv REPORT_B=.../B/report.tsv \
#     bash scripts/run_stage1.sh                                # build seeds first
#   DRY_RUN=1 bash scripts/run_stage1.sh                        # plan only
set -euo pipefail

REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO_ROOT"

DRY_RUN="${DRY_RUN:-0}"
PYTHON="${PYTHON:-python3}"
TIMSIM="${TIMSIM:-timsim}"
ENGINE="${ENGINE:-diann}"
SAMPLES=(A B)

SEEDS_DIR="simulations/stage1/seeds"
OUT_DIR="simulations/stage1/dia"
BASE_CFG="simulations/stage1/configs/base-dia.toml"

# Optionally build seeds from reports.
declare -A REPORTS=( [A]="${REPORT_A:-}" [B]="${REPORT_B:-}" )
for s in "${SAMPLES[@]}"; do
  if [[ -n "${REPORTS[$s]}" ]]; then
    echo "==> build seed for sample $s from ${REPORTS[$s]} (engine=$ENGINE)"
    $PYTHON scripts/build_seed_from_report.py \
      --engine "$ENGINE" --report "${REPORTS[$s]}" --sample "$s" \
      --out "$SEEDS_DIR/seed_sample$s.csv"
  fi
done

for s in "${SAMPLES[@]}"; do
  if [[ ! -f "$SEEDS_DIR/seed_sample$s.csv" ]]; then
    echo "Missing $SEEDS_DIR/seed_sample$s.csv — build it (build_seed_from_report.py) "\
         "or pass REPORT_$s=..." >&2
    exit 1
  fi
done

# Capability check: stage 1 needs the from_findings TimSim mode.
if [[ "$DRY_RUN" != "1" && "${SKIP_TIMSIM_CHECK:-0}" != "1" ]]; then
  if ! command -v "$TIMSIM" >/dev/null 2>&1; then
    echo "ERROR: '$TIMSIM' not on PATH. Build it: make venv-source && source .venv/bin/activate" >&2
    exit 2
  fi
  if ! "$TIMSIM" --help 2>&1 | grep -q -- "--findings-path"; then
    echo "ERROR: '$TIMSIM' lacks --findings-path (older build without from_findings)." >&2
    echo "       Build the pinned submodule: make venv-source" >&2
    exit 2
  fi
fi

echo "==> render configs for samples: ${SAMPLES[*]}"
PATH_MODE="${PATH_MODE:-host}"
$PYTHON scripts/render_configs.py \
  --base "$BASE_CFG" --seeds-dir "$SEEDS_DIR" --out-dir "$OUT_DIR" \
  --samples "${SAMPLES[@]}" --path-mode "$PATH_MODE"

# Provenance for evidence.
PROV="$OUT_DIR/RUN_PROVENANCE.txt"
{
  echo "Run started:    $(date -Is)"
  echo "Repo HEAD:      $(git rev-parse HEAD 2>/dev/null || echo n/a)"
  echo "rustims HEAD:   $(git -C rustims rev-parse HEAD 2>/dev/null || echo n/a)"
  echo "Engine (seed):  $ENGINE"
  echo "Samples:        ${SAMPLES[*]}"
} > "$PROV"
echo "==> wrote $PROV"

for s in "${SAMPLES[@]}"; do
  cfg="$OUT_DIR/sample$s/config.toml"
  echo
  echo "==> sample $s — $TIMSIM $cfg"
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "    (dry-run; skipping)"
  else
    $TIMSIM "$cfg"
  fi
done

echo
echo "Stage 1 done. Outputs under $OUT_DIR/sample{A,B}/"
