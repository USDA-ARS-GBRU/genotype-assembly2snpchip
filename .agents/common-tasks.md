# Common Agent Tasks

## Sanity check toolchain

```bash
bcftools --version
samtools --version
minimap2 --version
python -c "import pandas, matplotlib, seaborn, sklearn, openpyxl; print('Python packages OK')"
```

## Smoke tests

```bash
bash scripts/check_workflow_invariants.sh
bash tests/run_tiny_test.sh
bash tests/run_tiny_plot_test.sh
bash tests/run_tiny_pca_test.sh
bash tests/run_tiny_grin_test.sh
```

## Rebuild bundled example outputs

```bash
bash scripts/regenerate_example_outputs.sh
```

## Common edit surfaces

- Workflow flags: `sbatch/`, `docs/step-2...`, `docs/step-3...`, `examples/`
- Summary outputs: `scripts/summarize_gtcheck_top_hits.py`, `README.md`, `docs/`
- Plotting: `scripts/plot_gtcheck_summary.py`, `docs/visualization.md`, `examples/figures/`
- PCA/MDS: `scripts/plot_panel_pca_mds.py`, `docs/visualization.md`, `tests/run_tiny_pca_test.sh`

## PR review checklist

- Does the change preserve species-agnostic behavior?
- Are both loop and array sbatch scripts still aligned?
- Are docs consistent with actual script behavior?
- Does the change affect example figures or bundled outputs?
- Is there a lightweight validation step to run?

For a review-specific prompt, use `.agents/pr-review-agent.md`.
