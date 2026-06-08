# PlasmaBENCH owns this image. It is built from the pinned rustims submodule
# (feature/simulate-from-findings) by .github/workflows/docker-publish.yml and
# pushed to ghcr.io. We point at it instead of upstream rustims:latest so
# stage 1 is not blocked by main lagging behind from_findings.
IMAGE := ghcr.io/thegreatherrlebert/plasmabench-timsim:from-findings
# Pinned digest for benchmark-grade runs. After the first successful run of
# .github/workflows/docker-publish.yml, replace the value below with the digest
# the workflow logs (`Pushed image digest: ...`). Until then, `make pull-pinned`
# fails loudly — use `make pull` to grab the convenience tag.
IMAGE_DIGEST := ghcr.io/thegreatherrlebert/plasmabench-timsim@sha256:UNPINNED

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

# Evaluation entrypoint declared in evident.yaml. The actual eval pipeline
# (species-resolved ratio recovery, FDR, recall against the TimSim blueprint)
# is not implemented yet — this stub fails loudly so the manifest's
# reproducible command does not pretend to work.
stage1-eval:
	@echo "stage1-eval is not implemented yet. Tracking issue: TODO.md (Stage 1)." >&2
	@exit 64

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
