#!/usr/bin/env bash
# Fetch the DIA-NN Linux builds used as software-under-test (versions 1.8 and 2.5).
#
# Binaries are NOT committed (tools/ is gitignored); this script makes the
# download reproducible. Academic builds from the official GitHub releases
# (https://github.com/vdemichev/DiaNN). Review the bundled LICENSE.txt before use.
#
# Usage:  bash scripts/get_diann.sh
# Result: tools/diann/1.8.1/  and  tools/diann/2.5.0/  with the diann binary inside.
set -euo pipefail

REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
DEST="$REPO_ROOT/tools/diann"
mkdir -p "$DEST"

# version label -> release asset URL
declare -A URLS=(
  ["1.8.1"]="https://github.com/vdemichev/DiaNN/releases/download/1.8.1/diann_1.8.1.tar.gz"
  ["2.5.0"]="https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.5.0-Academia-Linux.zip"
)

# Expected CLI binary per version (validated before marking the install OK).
declare -A BINS=(
  ["1.8.1"]="1.8.1/diann-1.8.1"
  ["2.5.0"]="2.5.0/diann-2.5.0/diann-linux"
)

for ver in "${!URLS[@]}"; do
  url="${URLS[$ver]}"
  out="$DEST/$ver"
  if [[ -e "$out/.ok" ]]; then echo "==> $ver already present, skipping"; continue; fi
  mkdir -p "$out"
  archive="$DEST/$(basename "$url")"
  echo "==> downloading DIA-NN $ver"
  curl -L --fail --retry 3 -o "$archive" "$url"
  echo "==> extracting into $out"
  case "$archive" in
    *.tar.gz) tar -xzf "$archive" -C "$out" ;;
    *.zip)    unzip -q -o "$archive" -d "$out" ;;
  esac
  rm -f "$archive"

  # Validate: the expected binary must exist and actually print its version
  # banner. Some archives ship the binary without the exec bit (e.g. 2.5.0's
  # diann-linux) — fix it, then probe. Only then mark the install OK so a
  # broken/changed layout fails loudly instead of caching as "successful".
  bin="$DEST/${BINS[$ver]}"
  [[ -f "$bin" ]] || { echo "ERROR: expected binary not found: $bin" >&2; exit 3; }
  chmod +x "$bin"
  if ! timeout 60 "$bin" 2>&1 | grep -q "DIA-NN $ver"; then
    echo "ERROR: $bin did not report 'DIA-NN $ver' on launch — install rejected" >&2
    exit 3
  fi
  echo "    verified: $("$bin" 2>&1 | grep -m1 'DIA-NN')"
  touch "$out/.ok"
done

echo
echo "Installed & verified:"
for ver in "${!BINS[@]}"; do echo "  $ver -> $DEST/${BINS[$ver]}"; done
