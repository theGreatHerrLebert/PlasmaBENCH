# PlasmaBENCH owns this image. It is built from the pinned rustims submodule
# (feature/simulate-from-findings) by .github/workflows/docker-publish.yml and
# pushed to ghcr.io. We point at it instead of upstream rustims:latest so
# stage 1 is not blocked by main lagging behind from_findings.
IMAGE := ghcr.io/thegreatherrlebert/plasmabench-timsim:from-findings
# Pinned digest for benchmark-grade runs. After the first successful run of
# .github/workflows/docker-publish.yml, replace the value below with the digest
# the workflow logs (`Pushed image digest: ...`). Until then, `make pull-pinned`
# fails loudly — use `make pull` to grab the convenience tag.
IMAGE_DIGEST := ghcr.io/thegreatherrlebert/plasmabench-timsim@sha256:527ee78f48cc0431c2abf7adf9807c4efe9c39c9b843fc351977c6a746a7be21

.PHONY: help pull pull-pinned shell timsim-help compose-up compose-down venv venv-source submodules stage1-eval

help:
	@echo "PlasmaBENCH targets:"
	@echo "  make pull          – docker pull $(IMAGE)"
	@echo "  make pull-pinned   – docker pull $(IMAGE_DIGEST) (release-grade)"
	@echo "  make shell         – interactive shell in the rustims container (GPU)"
	@echo "  make timsim-help   – run 'timsim --help' inside the container"
	@echo "  make compose-up    – docker compose up (named services in docker/)"
	@echo "  make compose-down  – tear down compose stack"
	@echo "  make submodules    – init/update evident and rustims submodules"
	@echo "  make venv          – create local Python 3.12 venv (lightweight, no timsim)"
	@echo "  make venv-source   – build timsim from the pinned rustims submodule into .venv"
	@echo "  make stage1-eval   – (stub; exits 64) entrypoint declared in evident.yaml"

pull:
	docker pull $(IMAGE)

pull-pinned:
	@if echo "$(IMAGE_DIGEST)" | grep -q UNPINNED; then \
		echo "ERROR: IMAGE_DIGEST is still UNPINNED in Makefile." >&2; \
		echo "Run the docker-publish workflow once, then paste the" >&2; \
		echo "'Pushed image digest' value into Makefile's IMAGE_DIGEST." >&2; \
		exit 64; \
	fi
	docker pull $(IMAGE_DIGEST)

# PLB_DATA / PLB_OUT mirror docker/docker-compose.yml so `make shell` and
# `docker compose run` see the same /data and /sim contents. Override on the
# command line, e.g. `make shell PLB_DATA=/mnt/big/raw`.
PLB_DATA ?= $(CURDIR)/data
PLB_OUT  ?= $(CURDIR)/simulations

shell:
	docker run --rm -it --gpus all \
		-v $(CURDIR):/work \
		-v $(PLB_DATA):/data \
		-v $(PLB_OUT):/sim \
		-w /work \
		$(IMAGE) /bin/bash

timsim-help:
	docker run --rm --gpus all $(IMAGE) timsim --help

compose-up:
	cd docker && docker compose up -d

compose-down:
	cd docker && docker compose down

venv:
	python3.12 -m venv .venv
	. .venv/bin/activate && pip install --upgrade pip && pip install -r requirements-dev.txt
	@echo "Venv created at .venv (dev tools only). See docs/venv-setup.md to build timsim from source."

# --- EVIDENT replay commands -------------------------------------------------------
# Each *-eval target is the turnkey replay for one evident.yaml claim: it runs against
# concrete (gitignored, externally-provided) inputs and writes the cited artifact.
# Pure-python (pandas + plasmabench/plots), so they also run in the analysis image.
# RUNPY picks the host .venv if present, else the system python (the image case).
# Prefer the host .venv; else python3 (universally present + the slim image); else python.
RUNPY = PY=$$([ -x .venv/bin/python ] && echo .venv/bin/python || command -v python3 || command -v python) && "$$PY"

# plasmabench-stage1-simulated-truth — ratio recovery + FDR/recall vs blueprint.
stage1-eval:
	@$(RUNPY) scripts/stage1_eval.py
# plasmabench-stage2-superimpose-truth — VERIFIES blueprint B/A within 5% on both
# backgrounds (sim-QC gate) + writes metrics.json and the ratio panel.
stage2-eval:
	@$(RUNPY) scripts/stage2_eval.py
# plasmabench-diann-crossrun-normalization — the normalization sweep figure + the
# shared-precursor / common-column / residual-vs-truth controls (diag_quant_method).
quant-eval:
	@$(RUNPY) scripts/plot_quant_18v25.py
	@$(RUNPY) scripts/diag_quant_method.py --report results/diann-2.5-s2fix-r080/report.parquet \
	  --sim-dir simulations/stage2/dia-fix --experiment stage2 --human-scope background_unknown
# plasmabench-real-vs-sim-validation — the real-vs-SIM compression capstone figure.
real-vs-sim-eval:
	@$(RUNPY) scripts/plot_real_vs_sim_compression.py
# plasmabench-lowinput-overseparation — the per-member under-quantification figure.
lowinput-eval:
	@$(RUNPY) scripts/plot_lowabund_underquant.py
# plasmabench-ideoverlap-real — cross-engine real-data ID overlap (1.8 vs 2.5 AND 1.8 vs 2.6).
idoverlap-eval:
	@$(RUNPY) scripts/plot_id_overlap.py --report-a results/diann-1.8-realG/report.tsv \
	  --label-a "DIA-NN 1.8" --report-b results/diann-2.5-realG/report.parquet \
	  --label-b "DIA-NN 2.5" --title "Real PYE1 ID overlap — 1.8 vs 2.5" \
	  --out results/figures/real_overlap_18_vs_25.png
	@$(RUNPY) scripts/plot_id_overlap.py --report-a results/diann-1.8-realG/report.tsv \
	  --label-a "DIA-NN 1.8" --report-b results/diann-2.6-realG/report.parquet \
	  --label-b "DIA-NN 2.6" --title "Real PYE1 ID overlap — 1.8 vs 2.6" \
	  --out results/figures/real_overlap_18_vs_26.png

# --- Analysis-replay image (EVIDENT replay surface; code pinned, raw data external) ---
# Default = the CI-published image (.github/workflows/analysis-image.yml). Build locally
# instead with `make analysis-image ANALYSIS_IMAGE=plasmabench-analysis:local`
# (note: a snap-confined docker cannot read build contexts/mounts outside $HOME).
ANALYSIS_IMAGE ?= ghcr.io/thegreatherrlebert/plasmabench-analysis:latest
analysis-image:
	docker build -f docker/analysis.Dockerfile -t $(ANALYSIS_IMAGE) .

# Replay ANY *-eval target INSIDE the container with raw data mounted externally
# (results/ simulations/ data/ — gitignored, provided by you). Override roots with
# PLB_OUT (simulations) / PLB_DATA (data).  Usage: make eval-docker EVAL=lowinput-eval
EVAL ?= stage1-eval
EVAL_TARGETS = stage1-eval stage2-eval quant-eval real-vs-sim-eval lowinput-eval idoverlap-eval
eval-docker:
	@case " $(EVAL_TARGETS) " in *" $(EVAL) "*) ;; *) \
	  echo "eval-docker: unknown EVAL=$(EVAL); choose one of: $(EVAL_TARGETS)" >&2; exit 2 ;; esac
	docker run --rm \
	  -v "$(CURDIR)/results:/work/results" \
	  -v "$${PLB_OUT:-$(CURDIR)/simulations}:/work/simulations" \
	  -v "$${PLB_DATA:-$(CURDIR)/data}:/work/data" \
	  $(ANALYSIS_IMAGE) make $(EVAL)
# Back-compat alias.
stage1-eval-docker:
	@$(MAKE) eval-docker EVAL=stage1-eval

submodules:
	git submodule update --init --recursive

venv-source: submodules
	test -d .venv || python3.12 -m venv .venv
	. .venv/bin/activate && pip install --upgrade pip "maturin>=1.2,<2.0"
	. .venv/bin/activate && cd rustims && maturin build --release -m imspy_connector/Cargo.toml --out dist/
	. .venv/bin/activate && pip install rustims/dist/imspy_connector-*.whl
	. .venv/bin/activate && pip install ./rustims/packages/imspy-core ./rustims/packages/imspy-predictors ./rustims/packages/imspy-search ./rustims/packages/imspy-simulation
	. .venv/bin/activate && pip install -r requirements-dev.txt
	. .venv/bin/activate && timsim --help | head -5
