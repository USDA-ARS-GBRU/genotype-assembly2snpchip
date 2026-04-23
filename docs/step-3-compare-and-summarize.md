# Step 3: Compare to the Panel and Summarize Hits

## `bcftools gtcheck`

After panel-site genotyping, the workflow compares each query VCF to the panel:

```bash
bcftools gtcheck \
  -u GT,GT \
  -E 0 \
  --keep-refs \
  -O t \
  -o "$gtcheck_tsv" \
  -g "$panel_snps" \
  "$diploid_vcf"
```

Key points:

- `-u GT,GT` forces genotype-vs-genotype comparison
- `-E 0` keeps discordance interpretation simple
- `--keep-refs` preserves reference-only sites instead of limiting output to distinctive sites

## Summarize the Best Matches

```bash
python scripts/summarize_gtcheck_top_hits.py \
  -i 'results/*.gtcheck.tsv' \
  -n 10 \
  --min-sites 1000 \
  -o results/gtcheck_top10.tsv
```

This script ranks the best panel matches for each query sample.

Useful columns include:

- `query_sample`
- `panel_sample`
- `rank`
- `discordance`
- `sites_compared`
- `matching_genotypes`
- `match_fraction`
- `confidence_score`

The companion sample summary file is written by default to:

```text
results/gtcheck_top10.sample_summary.tsv
```

## Optional GRIN Metadata Enrichment

If your panel sample names include PI accessions or named cultivars, you can enrich the ranked table with GRIN metadata:

```bash
python scripts/enrich_gtcheck_top_hits_with_grin.py \
  --input results/gtcheck_top10.tsv \
  --crop soybean
```

This script can normalize accession formatting such as:

- `PI647962 -> PI 647962`
- `PI424405B -> PI 424405 B`

It writes:

```text
results/gtcheck_top10.grin_enriched.tsv
results/gtcheck_top10.grin_enriched.xlsx
results/gtcheck_top10.grin_cache.json
```

The enrichment output adds:

- `genotyped_sample`
- `PLANT NAME`
- `TAXONOMY`
- `ORIGIN`
- `GRIN ID`
- lookup-status columns
