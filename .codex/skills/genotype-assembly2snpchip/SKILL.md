---
name: genotype-assembly2snpchip
description: Use when working in the genotype-assembly2snpchip repository to edit the SNP-chip identity workflow, update sbatch scripts or Python utilities, maintain docs/examples, or interpret repo-specific workflow constraints for bcftools/gtcheck, GRIN enrichment, and PCA/MDS plotting.
---

# Genotype Assembly to SNP-Chip Panel

Use this skill when the task is specific to this repository's workflow, docs, examples, or HPC execution templates.

## Start Here

Read only what you need:

- `README.md` for the short repo contract and quick start
- `docs/` for the long-form workflow
- `sbatch/` for the supported SLURM templates
- `scripts/` for the canonical species-agnostic utilities
- `tests/` for tiny regression fixtures and smoke tests

Prefer `rg` for locating commands, options, and previously documented behavior.

## Repo Invariants

These are the important workflow constraints. Preserve them unless the user explicitly wants a workflow change.

- `bcftools >= 1.23` is required because `bcftools gtcheck` uses `--keep-refs`.
- The `bcftools call` stage uses bare `-i` (`--insert-missed`), not `-i 1`.
- Run `bcftools +fixploidy` before `bcftools gtcheck` so non-diploid GT rows are not skipped.
- The reference genome used to map assemblies must match, or be coordinate-compatible with, the SNP-chip panel reference.
- The Python scripts are intended to stay species-agnostic and filename-agnostic.
- The top-level repo docs are generic. Soybean-specific history belongs in `examples/`.

## Files That Matter

- `sbatch/map_assemblies_to_reference*.sbatch`
- `sbatch/call_panel_variants_and_gtcheck*.sbatch`
- `scripts/summarize_gtcheck_top_hits.py`
- `scripts/enrich_gtcheck_top_hits_with_grin.py`
- `scripts/plot_gtcheck_summary.py`
- `scripts/plot_panel_pca_mds.py`
- `docs/setup.md`
- `docs/step-2-prepare-panel-and-call-sites.md`
- `docs/step-3-compare-and-summarize.md`
- `docs/visualization.md`
- `examples/README*.md`
- `examples/soy50k_gtcheck_from_asm20.sbatch`

## Editing Rules for This Repo

- When changing core workflow behavior, update both loop and array sbatch scripts.
- If you change a command-line flag or workflow explanation, update the matching docs and soybean example notes in the same pass.
- If you change plotting logic, regenerate the affected example figures when practical.
- Do not reintroduce soybean-only assumptions into the generic scripts or docs.
- Keep the quick-start README terse and the deeper explanation in `docs/`.

## Validation Shortcuts

Use the smallest useful check:

- Shell syntax:
  - `bash -n sbatch/call_panel_variants_and_gtcheck.sbatch`
  - `bash -n sbatch/call_panel_variants_and_gtcheck_array.sbatch`
- Tiny tests:
  - `bash tests/run_tiny_test.sh`
  - `bash tests/run_tiny_plot_test.sh`
  - `bash tests/run_tiny_pca_test.sh`
  - `bash tests/run_tiny_grin_test.sh`
- Pixi tasks if preferred:
  - `pixi run test-parser`
  - `pixi run test-plots`
  - `pixi run test-pca`
  - `pixi run test-grin`

If a change only affects docs or agent scaffolding, lightweight review is enough.

## When You Need More Detail

Read these references only if the task needs them:

- `references/repo-map.md` for a quick file map and task routing
- `references/common-workflows.md` for the usual agent tasks in this repo

