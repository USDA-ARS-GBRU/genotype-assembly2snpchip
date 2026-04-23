# Visualization

## Plot the `gtcheck` Summary

After generating the ranked hit table and sample summary, create summary figures:

```bash
python scripts/plot_gtcheck_summary.py \
  --top-hits results/gtcheck_top10.tsv \
  --sample-summary results/gtcheck_top10.sample_summary.tsv \
  --out-dir figures \
  --prefix gtcheck
```

Typical outputs:

```text
figures/gtcheck_top_hits_lollipop.png
figures/gtcheck_rank1_rank2_gap.png
figures/gtcheck_match_fraction_vs_sites.png
figures/gtcheck_top_hit_heatmap.png
```

These plots help answer:

- is there a clean rank-1 hit?
- how large is the rank-1 vs rank-2 gap?
- how many sites support the match?
- do multiple queries share the same top-hit neighborhood?

## PCA and MDS Context

You can also place assembly-derived query genotypes into the context of the full SNP-chip panel:

```bash
python scripts/plot_panel_pca_mds.py \
  --panel-vcf work/panel/panel.biallelic.snps.vcf.gz \
  --query-vcfs results/*.panel.filtered.diploid.vcf.gz \
  --out-dir figures \
  --prefix panel_context \
  --method pca
```

For smaller focused subsets, you can request both PCA and MDS:

```bash
python scripts/plot_panel_pca_mds.py \
  --panel-vcf work/panel/panel.biallelic.snps.vcf.gz \
  --query-vcfs results/*.panel.filtered.diploid.vcf.gz \
  --out-dir figures \
  --prefix panel_context \
  --method both \
  --max-mds-samples 500
```

## PCA or MDS?

- Use PCA as the default for most projects
- Use MDS for smaller, focused subsets when a distance-style view is helpful
- Treat both as supportive context, not the primary identity call

## Common Plot Issues

| Symptom | Likely cause | What to check |
| --- | --- | --- |
| Query is a strong outlier | weak marker overlap, wrong reference, contamination, mislabeled sample | compare `sites_compared`, inspect call rate, confirm the panel reference genome |
| Query does not appear | no overlapping sites, filtered markers, unexpected VCF sample name | inspect the coordinates TSV and VCF header |
| PCA is dominated by one sample | low-quality or very divergent sample | inspect PC2/PC3, use a smaller subset |
| MDS is slow | pairwise scaling problem on large panels | use PCA or subset the panel |
| Colors or labels are missing | metadata sample names do not match | compare metadata names to the coordinates TSV |
