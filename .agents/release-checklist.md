# Release Checklist

Use this before tagging a stable version or announcing a substantial workflow update.

## Workflow integrity

- Run `bash scripts/check_workflow_invariants.sh`
- Run:
  - `bash tests/run_tiny_test.sh`
  - `bash tests/run_tiny_plot_test.sh`
  - `bash tests/run_tiny_pca_test.sh`
  - `bash tests/run_tiny_grin_test.sh`
- If workflow behavior changed, verify both:
  - `sbatch/call_panel_variants_and_gtcheck.sbatch`
  - `sbatch/call_panel_variants_and_gtcheck_array.sbatch`

## Documentation integrity

- Confirm `README.md` still matches the actual commands and expected outputs
- Confirm `docs/setup.md`, `docs/step-2-prepare-panel-and-call-sites.md`, and `docs/step-3-compare-and-summarize.md` match the current scripts
- Confirm `docs/visualization.md` matches current plotting behavior and filenames
- Confirm any environment or version requirement changes are reflected in:
  - `environment.yml`
  - `pixi.toml`
  - `requirements.txt` if relevant
  - `docs/hpc_notes.md`

## Example integrity

- Run `bash scripts/regenerate_example_outputs.sh` if changes affect:
  - summary outputs
  - GRIN enrichment
  - plotting
- Review whether any regenerated tracked example files should be committed deliberately
- Confirm soybean-specific notes in `examples/` still agree with the generic workflow where they overlap

## Release hygiene

- Check `git status`
- Review `git diff --stat`
- Make sure commits are intentional and scoped
- If publishing a release, prepare short notes summarizing:
  - workflow changes
  - doc changes
  - breaking assumptions or environment changes
