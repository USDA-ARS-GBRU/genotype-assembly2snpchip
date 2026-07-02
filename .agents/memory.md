# Repo Memory

## Durable facts

- Repo purpose: validate genome assembly identity against SNP-chip or marker-panel VCF genotypes.
- Core pipeline:
  1. map assemblies to the panel reference
  2. call genotypes only at panel sites
  3. fix ploidy
  4. run `bcftools gtcheck`
  5. summarize, enrich, and visualize

## Hard constraints

- `bcftools >= 1.23`
- `bcftools gtcheck --keep-refs`
- `bcftools call -i` with no integer
- `bcftools +fixploidy` before `gtcheck`
- reference genome used for mapping must match the SNP-chip panel reference space

## Canonical outputs

- `results/gtcheck_top*.tsv`
- `results/gtcheck_top*.sample_summary.tsv`
- `results/gtcheck_top*.grin_enriched.tsv`
- `results/gtcheck_top*.grin_enriched.xlsx`
- `figures/gtcheck_*.png`
- `figures/panel_context_*.png`

## Repo maintenance helpers

- `scripts/check_workflow_invariants.sh`
- `scripts/regenerate_example_outputs.sh`
- `.agents/pr-review-agent.md`
- `.agents/release-checklist.md`

## Canonical scripts

- `scripts/summarize_gtcheck_top_hits.py`
- `scripts/enrich_gtcheck_top_hits_with_grin.py`
- `scripts/plot_gtcheck_summary.py`
- `scripts/plot_panel_pca_mds.py`

## Change-coupling reminders

- Workflow flag changes usually require updates in:
  - `sbatch/`
  - `docs/`
  - `examples/`
- Plot behavior changes usually require:
  - regenerated example figures
  - updated captions in `README.md` and `docs/visualization.md`

## Common pitfalls

- old cluster modules, especially old `bcftools`
- soybean-specific language leaking into generic docs
- changing only one sbatch variant
- forgetting to update example docs after workflow changes
- interpreting strong match fraction without checking `sites_compared`
