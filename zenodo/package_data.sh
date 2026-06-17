#!/usr/bin/env bash
# Package the PlasmaBENCH data archive for Zenodo.
#
# Produces three .tar.zst tarballs (+ SHA256SUMS) under $OUT_DIR that unpack
# INTO the repo root, recreating data/raw, simulations/, results/.
#
# Archive scope (decided with the maintainer):
#   - data/raw/        : FASTA, blank .d, real G PYE1 .d's, FragPipe DDA inputs
#   - simulations/     : simulated .d's + synthetic_data.db blueprints + seeds + configs
#   - results/         : DIA-NN reports  MINUS  *.predicted.speclib  (the ~20 GB
#                        predicted spectral libraries are EXCLUDED by request)
#
# Engine binaries are never included (users bring their own DIA-NN / FragPipe).
#
# Usage:
#   bash zenodo/package_data.sh                 # writes to ./_zenodo_pkg
#   OUT_DIR=/big/disk bash zenodo/package_data.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

OUT_DIR="${OUT_DIR:-$REPO_ROOT/_zenodo_pkg}"
ZSTD_LEVEL="${ZSTD_LEVEL:-19}"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 4)}"
mkdir -p "$OUT_DIR"

# zstd with multithreading; --long improves ratio on the big .d trees.
ZSTD=(zstd "-${ZSTD_LEVEL}" "-T${THREADS}" --long=27 -q)

say() { printf '\033[1;36m==>\033[0m %s\n' "$*"; }

# --- 1. data/raw -------------------------------------------------------------
say "Packing data/raw  ($(du -sh data/raw 2>/dev/null | cut -f1))"
tar -cf - data/raw | "${ZSTD[@]}" -o "$OUT_DIR/plasmabench-data-raw.tar.zst" -f

# --- 2. simulations ----------------------------------------------------------
say "Packing simulations  ($(du -sh simulations 2>/dev/null | cut -f1))"
tar -cf - simulations | "${ZSTD[@]}" -o "$OUT_DIR/plasmabench-simulations.tar.zst" -f

# --- 3. results (excluding predicted spectral libraries) ---------------------
# *.predicted.speclib is the ~20 GB DIA-NN predicted library — regenerable, excluded.
say "Packing results MINUS *.predicted.speclib"
tar -cf - --exclude='*.predicted.speclib' results \
  | "${ZSTD[@]}" -o "$OUT_DIR/plasmabench-results.tar.zst" -f

# --- checksums + manifest ----------------------------------------------------
say "Writing checksums"
( cd "$OUT_DIR" && sha256sum plasmabench-*.tar.zst > SHA256SUMS )

say "Done. Archive contents:"
( cd "$OUT_DIR" && du -h plasmabench-*.tar.zst && echo && cat SHA256SUMS )

cat <<EOF

Next:
  export ZENODO_TOKEN=...        # personal access token (deposit:write scope)
  bash zenodo/upload_to_zenodo.sh
EOF
