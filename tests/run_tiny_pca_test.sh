#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

out_dir="tests/tmp/pca_figures"
PYTHON="${PYTHON:-python3}"

mkdir -p "$out_dir"
export MPLCONFIGDIR="tests/tmp/matplotlib"
mkdir -p "$MPLCONFIGDIR"

"$PYTHON" -B scripts/plot_panel_pca_mds.py \
    --panel-vcf tests/fixtures/tiny_panel.vcf \
    --query-vcfs tests/fixtures/tiny_query_A.vcf tests/fixtures/tiny_query_C.vcf \
    --metadata tests/fixtures/tiny_pca_metadata.tsv \
    --out-dir "$out_dir" \
    --prefix tiny \
    --method both \
    --min-site-call-rate 0.5 \
    --min-maf 0.01 \
    --formats svg >/dev/null

for file in \
    tiny_pca_coordinates.tsv \
    tiny_pca_pc1_pc2.svg \
    tiny_mds_coordinates.tsv \
    tiny_mds1_mds2.svg \
    tiny_matrix_summary.tsv
do
    test -s "$out_dir/$file"
done

grep -q 'query_A' "$out_dir/tiny_pca_coordinates.tsv"
grep -q '<svg' "$out_dir/tiny_pca_pc1_pc2.svg"
grep -q '<svg' "$out_dir/tiny_mds1_mds2.svg"

echo "Tiny PCA/MDS plot test passed."
