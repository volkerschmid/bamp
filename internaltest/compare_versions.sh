#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"

VERSIONS=("bamp-1" "bamp")
OUT_DIR="$SCRIPT_DIR/compare_out"
mkdir -p "$OUT_DIR"

for v in "${VERSIONS[@]}"; do
  pkg_path="$ROOT/$v"
  out_file="$OUT_DIR/${v}.txt"

  if [[ ! -d "$pkg_path" ]]; then
    echo "ERROR: $pkg_path not found" >&2
    exit 1
  fi

  echo "=== Running $v ==="
  start=$(date +%s%3N)
  "$RSCRIPT" --vanilla "$SCRIPT_DIR/run_one_version.R" "$pkg_path" \
    > "$out_file" 2>&1
  end=$(date +%s%3N)
  wall=$(( end - start ))

  echo "  Wall time: ${wall} ms"
  grep "^TIME" "$out_file" || true
  echo "  Output saved to: $out_file"
  echo
done

echo "=== Diff (bamp vs bamp-1) ==="
diff --unified \
  <(grep -v "^TIME" "$OUT_DIR/bamp.txt") \
  <(grep -v "^TIME" "$OUT_DIR/bamp-1.txt") \
  || true

echo
echo "=== Timing summary ==="
for v in "${VERSIONS[@]}"; do
  echo -n "  $v: "
  grep "^TIME" "$OUT_DIR/${v}.txt" | tr '\n' '  ' || true
  echo
done
