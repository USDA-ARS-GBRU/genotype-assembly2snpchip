# Contributing

Thanks for contributing to `genotype-assembly2snpchip`.

This repository mixes workflow scripts, teaching documentation, bundled examples, and agent-maintainer tooling. The most helpful contributions keep those surfaces synchronized.

## Before opening a pull request

Run the smallest relevant checks for your change.

At minimum:

```bash
bash scripts/check_workflow_invariants.sh
```

Common local checks:

```bash
bash tests/run_tiny_test.sh
bash tests/run_tiny_plot_test.sh
bash tests/run_tiny_pca_test.sh
bash tests/run_tiny_grin_test.sh
```

If you changed bundled soybean examples, you may also want:

```bash
bash scripts/regenerate_example_outputs.sh
```

## Contribution guidelines

- Keep the generic workflow species-agnostic.
- Keep soybean-specific history or teaching notes in `examples/`.
- If you change a core workflow flag or command, update:
  - the relevant `sbatch/` script
  - the matching docs in `docs/`
  - the matching soybean example notes in `examples/` when applicable
- If you change panel genotyping or `gtcheck`, update both:
  - `sbatch/call_panel_variants_and_gtcheck.sbatch`
  - `sbatch/call_panel_variants_and_gtcheck_array.sbatch`
- If you change plotting behavior, check whether bundled example figures or captions should be updated.

## Important workflow assumptions

Please do not change these casually:

- `bcftools >= 1.23`
- `bcftools call -i` with no integer
- `bcftools gtcheck --keep-refs`
- `bcftools +fixploidy` before `gtcheck`
- the mapping reference must match the SNP-chip panel reference space

## Pull requests

Use the PR template. It asks for:

- what changed
- why it changed
- what checks you ran
- whether docs/examples/sbatch variants were updated

## Issues

Please use the GitHub issue templates for:

- bug reports
- workflow questions
- documentation suggestions

## Maintainers and agents

There is repo-local maintainer scaffolding in:

- `.agents/`
- `.codex/skills/`

Helpful entry points:

- `.agents/common-tasks.md`
- `.agents/pr-review-agent.md`
- `.agents/release-checklist.md`

