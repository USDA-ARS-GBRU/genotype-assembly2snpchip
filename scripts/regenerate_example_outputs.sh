#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ -n "${PYTHON:-}" ]]; then
    python_bin="$PYTHON"
elif [[ -x ".venv/bin/python" ]]; then
    python_bin=".venv/bin/python"
else
    python_bin="python3"
fi

results_dir="examples/results"
figures_dir="examples/figures"
prefix="soy50k_example"
top_hits="$results_dir/example_gtcheck_top3.tsv"
sample_summary="$results_dir/example_gtcheck_top3.sample_summary.tsv"
enriched_tsv="$results_dir/example_gtcheck_top3.grin_enriched.tsv"
enriched_xlsx="$results_dir/example_gtcheck_top3.grin_enriched.xlsx"
cache_json="$results_dir/example_gtcheck_top3.grin_cache.json"

mkdir -p "$results_dir" "$figures_dir" tests/tmp/matplotlib
export MPLCONFIGDIR="tests/tmp/matplotlib"

"$python_bin" -B scripts/summarize_gtcheck_top_hits.py \
    -i 'examples/results/*.soy50k.fixploidy.gtcheck.tsv' \
    -n 3 \
    --min-sites 1000 \
    -o "$top_hits"

"$python_bin" -B scripts/enrich_gtcheck_top_hits_with_grin.py \
    --input "$top_hits" \
    --cache-json "$cache_json" \
    --output-tsv "$enriched_tsv" \
    --output-xlsx "$enriched_xlsx" \
    --offline

"$python_bin" -B scripts/plot_gtcheck_summary.py \
    --top-hits "$top_hits" \
    --sample-summary "$sample_summary" \
    --out-dir "$figures_dir" \
    --prefix "$prefix"

echo "Regenerated example outputs:"
echo "  $top_hits"
echo "  $sample_summary"
echo "  $enriched_tsv"
echo "  $enriched_xlsx"
echo "  $cache_json"
echo "  $figures_dir/${prefix}_top_hits_lollipop.png"

