# PlasmaBENCH analysis-replay image (CPU, slim).
#
# The EVIDENT replay surface for the analysis claims: the CODE + environment are
# baked in (pinned), and the RAW DATA (gitignored DIA-NN reports + SIM blueprints)
# is provided EXTERNALLY at run time via -v mounts onto results/ simulations/ data/.
# This is NOT the timsim runtime (no GPU, no rustims) — it only scores/plots reports.
#
# Build:  make analysis-image
# Run  :  make stage1-eval-docker   (mounts the external data, runs `make stage1-eval`)
FROM python:3.12-slim

# Analysis deps only (no jupyter/sagepy — the DIA-NN scorer + plots don't need them).
RUN pip install --no-cache-dir \
    "pandas>=2.0" "numpy>=1.24" "pyarrow>=14.0" "matplotlib>=3.7" "PyYAML>=6.0" \
 && apt-get update && apt-get install -y --no-install-recommends make \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /work
# Pinned code + manifest (data dirs are mounted, never copied).
COPY plasmabench/ ./plasmabench/
COPY scripts/ ./scripts/
COPY configs/ ./configs/
COPY Makefile evident.yaml ./

# `make stage1-eval` runs system python here (no .venv); raw data arrives via mounts.
CMD ["make", "stage1-eval"]
