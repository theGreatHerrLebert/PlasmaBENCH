# Local venv setup (no Docker)

Use this when iterating on the simulator source, or when Docker isn't
available. `imspy-simulation` is not on PyPI — it builds from the rustims
monorepo (vendored as the `rustims/` submodule) and depends on a Rust+PyO3
wheel (`imspy_connector`).

## Prerequisites
- Python 3.12 (3.11–3.13 supported)
- Rust toolchain (`rustup`, stable)
- `maturin` (installed into the venv)

## Steps

```bash
cd /path/to/PlasmaBENCH

# Check out the rustims submodule at the pinned (from_findings) commit.
git submodule update --init --recursive rustims

python3.12 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip "maturin>=1.2,<2.0"

# Build the PyO3 connector wheel from the pinned submodule.
cd rustims
maturin build --release -m imspy_connector/Cargo.toml --out dist/

# Install connector + simulation packages.
pip install dist/imspy_connector-*.whl
pip install ./packages/imspy-core ./packages/imspy-predictors \
            ./packages/imspy-search ./packages/imspy-simulation

cd ..
pip install -r requirements-dev.txt   # pandas, pyarrow, sagepy (seed builder)
timsim --help | head -5               # should list --findings-path
```

`make venv-source` runs all of the above.

## Bumping the rustims pin

```bash
cd rustims && git fetch origin && git checkout <new-tag-or-sha> && cd ..
git add rustims && git commit -m "Bump rustims submodule to <new-tag-or-sha>"
```

Bump the Docker image's `IMAGE_DIGEST` in the Makefile at the same time so the
from-source and Docker paths stay in lockstep.
