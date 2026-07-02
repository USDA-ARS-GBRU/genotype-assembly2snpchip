# Maintainer Notes

This page is for repository maintainers and agents, not first-time workflow users.

## Purpose

The goal is to keep workflow logic, docs, examples, and maintainer tooling aligned.

This page also serves as a lightweight Pages rebuild surface when maintainer-facing docs need to change without touching the user workflow pages.

## Most useful files

- `scripts/check_workflow_invariants.sh`
- `scripts/regenerate_example_outputs.sh`
- `.agents/common-tasks.md`
- `.agents/pr-review-agent.md`
- `.agents/release-checklist.md`
- `.agents/memory.md`
- `.codex/skills/genotype-assembly2snpchip/`

## Quick maintenance commands

### Invariant check

```bash
bash scripts/check_workflow_invariants.sh
```

### Tiny smoke tests

```bash
bash tests/run_tiny_test.sh
bash tests/run_tiny_plot_test.sh
bash tests/run_tiny_pca_test.sh
bash tests/run_tiny_grin_test.sh
```

### Rebuild bundled example outputs

```bash
bash scripts/regenerate_example_outputs.sh
```

## Contribution and review reminders

- Keep the generic workflow species-agnostic.
- Update both loop and array sbatch scripts when changing the core panel workflow.
- When changing workflow commands or flags, update the matching docs and soybean example notes.
- When changing plotting behavior, decide deliberately whether bundled example figures should be regenerated and committed.

## Important repo invariants

- `bcftools >= 1.23`
- `bcftools call -i` with no integer
- `bcftools gtcheck --keep-refs`
- `bcftools +fixploidy` before `gtcheck`
- the mapping reference must match the SNP-chip panel reference space

## CI and GitHub scaffolding

The repository now includes:

- GitHub Actions checks in `.github/workflows/repo-checks.yml`
- issue templates in `.github/ISSUE_TEMPLATE/`
- a pull request template in `.github/pull_request_template.md`

## Recommended release prep

Use:

- `.agents/release-checklist.md`

That checklist is the shortest path to a clean release or tagged workflow update.
