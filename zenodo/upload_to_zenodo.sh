#!/usr/bin/env bash
# Upload the packaged PlasmaBENCH data to an EXISTING Zenodo deposition and
# (optionally) set its metadata. Leaves the deposition as a DRAFT — publishing
# (which mints the DOI permanently) is a separate, explicit step you do in the
# Zenodo web UI after reviewing, OR by passing ZENODO_PUBLISH=1.
#
# The deposition already exists: DOI 10.5281/zenodo.20733913 (id 20733913).
#
# Prereqs:
#   - bash zenodo/package_data.sh   (creates ./_zenodo_pkg/*.tar.zst + SHA256SUMS)
#   - a Zenodo personal access token with scope deposit:write
#
# Usage:
#   export ZENODO_TOKEN=xxxxxxxx
#   bash zenodo/upload_to_zenodo.sh                 # upload files to draft 20733913
#   SET_METADATA=1 bash zenodo/upload_to_zenodo.sh  # also push zenodo-metadata.json
#   ZENODO_SANDBOX=1 ... bash zenodo/upload_to_zenodo.sh   # target sandbox.zenodo.org
#   ZENODO_PUBLISH=1 ...                            # publish after upload (IRREVERSIBLE)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-$HERE/../_zenodo_pkg}"
DEPOSITION_ID="${DEPOSITION_ID:-20733913}"

: "${ZENODO_TOKEN:?Set ZENODO_TOKEN (personal access token, scope deposit:write)}"

if [[ "${ZENODO_SANDBOX:-0}" == "1" ]]; then
  BASE="https://sandbox.zenodo.org"
else
  BASE="https://zenodo.org"
fi
API="$BASE/api/deposit/depositions/$DEPOSITION_ID"
AUTH=(-H "Authorization: Bearer $ZENODO_TOKEN")

say()  { printf '\033[1;36m==>\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null || die "missing dependency: $1"; }
need curl

# --- 1. fetch the deposition, read its upload bucket -------------------------
say "Fetching deposition $DEPOSITION_ID from $BASE"
DEP_JSON="$(curl -sf "${AUTH[@]}" "$API")" \
  || die "could not fetch deposition (token wrong / not your deposition / wrong instance?)"

# Prefer python for robust JSON parsing; fall back to grep.
read_field() {
  local key="$1"
  if command -v python3 >/dev/null; then
    printf '%s' "$DEP_JSON" | python3 -c "import sys,json;d=json.load(sys.stdin);print(d$key)"
  else
    printf '%s' "$DEP_JSON" | grep -oE "\"bucket\"[^,]*" | head -1 | sed -E 's/.*"(https[^"]+)".*/\1/'
  fi
}
BUCKET="$(read_field "['links']['bucket']")"
[[ "$BUCKET" == https* ]] || die "no bucket link on deposition (is it published/locked? create a new version draft first)"
say "Upload bucket: $BUCKET"

# --- 2. optional metadata ----------------------------------------------------
if [[ "${SET_METADATA:-0}" == "1" ]]; then
  META="$HERE/zenodo-metadata.json"
  [[ -f "$META" ]] || die "metadata file not found: $META"
  say "Setting metadata from $META"
  curl -sf "${AUTH[@]}" -H "Content-Type: application/json" \
       -X PUT "$API" -d @"$META" >/dev/null \
    && say "metadata updated" || die "metadata PUT failed (check zenodo-metadata.json)"
fi

# --- 3. upload files via the bucket API (handles multi-GB) -------------------
shopt -s nullglob
FILES=("$PKG_DIR"/plasmabench-*.tar.zst "$PKG_DIR"/SHA256SUMS)
[[ ${#FILES[@]} -gt 0 ]] || die "no files in $PKG_DIR — run zenodo/package_data.sh first"

for f in "${FILES[@]}"; do
  name="$(basename "$f")"
  size="$(du -h "$f" | cut -f1)"
  say "Uploading $name ($size) ..."
  curl --fail-with-body -# "${AUTH[@]}" \
       --upload-file "$f" "$BUCKET/$name" -o /dev/null \
    || die "upload of $name failed"
  say "  ✓ $name"
done

# --- 4. optional publish (irreversible) -------------------------------------
if [[ "${ZENODO_PUBLISH:-0}" == "1" ]]; then
  say "Publishing deposition (this mints the DOI permanently) ..."
  curl -sf "${AUTH[@]}" -X POST "$API/actions/publish" >/dev/null \
    && say "PUBLISHED → $BASE/record/$DEPOSITION_ID" \
    || die "publish failed"
else
  cat <<EOF

Files uploaded to the DRAFT deposition. Nothing is public yet.
Review and publish here:  $BASE/deposit/$DEPOSITION_ID
(or re-run with ZENODO_PUBLISH=1 to publish from the CLI.)
EOF
fi
