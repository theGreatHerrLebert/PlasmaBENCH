# Zenodo upload kit

Prepares and uploads the PlasmaBENCH data archive to Zenodo. The data is too
large for git (~44 GB) and is hosted as a Zenodo record instead.

**Target record:** [`10.5281/zenodo.20733913`](https://doi.org/10.5281/zenodo.20733913) (deposition id `20733913`).

## Files

| File | What it does |
|---|---|
| `package_data.sh` | Tars `data/raw/`, `simulations/`, and `results/` (minus `*.predicted.speclib`) into three `.tar.zst` + `SHA256SUMS` under `_zenodo_pkg/`. |
| `upload_to_zenodo.sh` | Uploads those files to the existing Zenodo draft `20733913` via the bucket API; optionally sets metadata and publishes. |
| `zenodo-metadata.json` | Deposition metadata (title, description, creators, license, related identifiers). **Has `FILL_ME` fields** — complete them before publishing. |
| `ARCHIVE_README.md` | The README to include *inside* the Zenodo record (data layout + how it maps to the repo). |

## Archive scope

Decided with the maintainer: **replay inputs + our DIA-NN reports, minus
predicted spectral libraries.** Engine binaries are never included — users bring
their own DIA-NN / FragPipe (see `REPRODUCE.md`).

| Tarball | ~size | Unpacks to |
|---|---|---|
| `plasmabench-data-raw.tar.zst` | ~17 GB | `data/raw/` |
| `plasmabench-simulations.tar.zst` | ~24 GB | `simulations/` |
| `plasmabench-results.tar.zst` | ~3 GB | `results/` (no `*.predicted.speclib`) |

## How to run it

```bash
# 1. complete the FILL_ME fields in zenodo-metadata.json (creators/affiliations)

# 2. package (writes ./_zenodo_pkg). Put OUT_DIR on a disk with ~45 GB free.
OUT_DIR=/big/disk/_zenodo_pkg bash zenodo/package_data.sh

# 3. get a Zenodo token: zenodo.org → Applications → Personal access tokens
#    (scope: deposit:write). Then:
export ZENODO_TOKEN=xxxxxxxx
PKG_DIR=/big/disk/_zenodo_pkg SET_METADATA=1 bash zenodo/upload_to_zenodo.sh

# 4. review at https://zenodo.org/deposit/20733913 and Publish in the UI
#    (or re-run with ZENODO_PUBLISH=1). Publishing is irreversible.
```

### Dry-run on the sandbox first (recommended)

```bash
# create a throwaway deposition on sandbox.zenodo.org, note its id, then:
ZENODO_SANDBOX=1 DEPOSITION_ID=<sandbox-id> ZENODO_TOKEN=<sandbox-token> \
  SET_METADATA=1 bash zenodo/upload_to_zenodo.sh
```

## Notes

- The upload uses Zenodo's **bucket API** (`--upload-file`), which streams and
  handles multi-GB files; no per-file 100 MB limit.
- The script **never publishes** unless you pass `ZENODO_PUBLISH=1` — by default
  it leaves a reviewable draft.
- Re-uploading a file with the same name overwrites it in the draft. To add a new
  *version* to an already-published record, create a new-version draft in the UI
  first, then point `DEPOSITION_ID` at it.
