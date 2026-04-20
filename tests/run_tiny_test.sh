#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

out_dir="tests/tmp"
out_file="$out_dir/tiny_top_hits.tsv"
PYTHON="${PYTHON:-python3}"

mkdir -p "$out_dir"

"$PYTHON" -B scripts/summarize_gtcheck_top_hits.py \
    -i tests/fixtures/tiny.gtcheck.tsv \
    -n 2 \
    --min-sites 10 \
    --sample-summary none \
    -o "$out_file" >/dev/null

diff -u tests/expected/tiny_top_hits.tsv "$out_file"

echo "Tiny gtcheck parser test passed."
