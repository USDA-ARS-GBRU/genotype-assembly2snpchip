#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

out_dir="tests/tmp/grin"
out_tsv="$out_dir/tiny_grin_enriched.tsv"
out_xlsx="$out_dir/tiny_grin_enriched.xlsx"
out_cache="$out_dir/tiny_grin_cache.json"
PYTHON="${PYTHON:-python3}"

mkdir -p "$out_dir"

cp tests/fixtures/tiny_grin_cache.json "$out_cache"

"$PYTHON" -B scripts/enrich_gtcheck_top_hits_with_grin.py \
    --input tests/fixtures/tiny_grin_input.tsv \
    --cache-json "$out_cache" \
    --output-tsv "$out_tsv" \
    --output-xlsx "$out_xlsx" \
    --offline >/dev/null

diff -u tests/expected/tiny_grin_enriched.tsv "$out_tsv"
test -s "$out_xlsx"

echo "Tiny GRIN enrichment test passed."
