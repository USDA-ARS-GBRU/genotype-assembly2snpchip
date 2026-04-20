#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

out_dir="tests/tmp"
top_hits="$out_dir/tiny_top_hits.tsv"
summary="$out_dir/tiny_top_hits.sample_summary.tsv"
fig_dir="$out_dir/figures"
PYTHON="${PYTHON:-python3}"

mkdir -p "$out_dir" "$fig_dir"
export MPLCONFIGDIR="$out_dir/matplotlib"
mkdir -p "$MPLCONFIGDIR"

"$PYTHON" -B scripts/summarize_gtcheck_top_hits.py \
    -i tests/fixtures/tiny.gtcheck.tsv \
    -n 2 \
    --min-sites 10 \
    -o "$top_hits" >/dev/null

"$PYTHON" -B scripts/plot_gtcheck_summary.py \
    --top-hits "$top_hits" \
    --sample-summary "$summary" \
    --out-dir "$fig_dir" \
    --prefix tiny \
    --formats svg >/dev/null

for fig in \
    tiny_top_hits_lollipop.svg \
    tiny_rank1_rank2_gap.svg \
    tiny_match_fraction_vs_sites.svg \
    tiny_top_hit_heatmap.svg
do
    test -s "$fig_dir/$fig"
    grep -q '<svg' "$fig_dir/$fig"
done

echo "Tiny gtcheck plot test passed."
