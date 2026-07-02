# Repo Map

## Core workflow

- `sbatch/map_assemblies_to_reference.sbatch`
- `sbatch/map_assemblies_to_reference_array.sbatch`
- `sbatch/call_panel_variants_and_gtcheck.sbatch`
- `sbatch/call_panel_variants_and_gtcheck_array.sbatch`

## Python utilities

- `scripts/summarize_gtcheck_top_hits.py`
  - Ranks gtcheck hits and writes a per-query sample summary table.
- `scripts/enrich_gtcheck_top_hits_with_grin.py`
  - Normalizes accession names and enriches tables with GRIN metadata.
- `scripts/plot_gtcheck_summary.py`
  - Builds lollipop, rank-gap, match-vs-sites, and heatmap summary figures.
- `scripts/plot_panel_pca_mds.py`
  - Builds PCA and optional MDS context figures directly from VCF input.
- `scripts/check_workflow_invariants.sh`
  - Verifies repo-level workflow invariants across sbatch scripts and docs.
- `scripts/regenerate_example_outputs.sh`
  - Rebuilds the bundled soybean example summaries, enrichment outputs, and figures.

## Documentation

- `README.md`
  - Brief front door, install, quick start, and example outputs.
- `docs/`
  - Long-form workflow, setup, interpretation, and visualization docs.

## Soybean case-study material

- `examples/`
  - Historical soybean-specific notes and example outputs. Keep this secondary to the generic workflow.

## Test fixtures

- `tests/fixtures/`
  - Tiny VCF, gtcheck, and GRIN input fixtures.
- `tests/run_tiny_*.sh`
  - Small smoke tests for parser, plotting, PCA/MDS, and GRIN enrichment.

## Agent-facing priorities

When a task crosses multiple surfaces, prefer this update order:

1. fix the canonical script or sbatch file
2. update the generic docs
3. update example or soybean-specific notes if they mention the same behavior
