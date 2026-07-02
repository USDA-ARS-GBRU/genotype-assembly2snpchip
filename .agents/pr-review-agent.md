# PR Review Agent

Use this prompt when reviewing pull requests in `genotype-assembly2snpchip`.

## Review stance

Lead with findings, not summary. Focus on:

- workflow regressions
- scientific assumption drift
- mismatches between scripts and docs
- missing updates across loop and array sbatch variants
- missing example or test updates when behavior changes

## Highest-priority checks

1. Does the PR preserve species-agnostic behavior?
2. If it changes panel genotyping or `gtcheck`, did it update both:
   - `sbatch/call_panel_variants_and_gtcheck.sbatch`
   - `sbatch/call_panel_variants_and_gtcheck_array.sbatch`
3. If it changes a workflow flag or command, do `docs/` and `examples/` still match?
4. If it changes plotting behavior, do the example figures or captions need regeneration?
5. If it changes tool requirements, do `README.md`, `docs/setup.md`, and `docs/hpc_notes.md` still agree?

## Invariants to protect

- `bcftools >= 1.23`
- `bcftools call -i` with no integer
- `bcftools gtcheck --keep-refs`
- `bcftools +fixploidy` before `gtcheck`
- mapping reference must match the SNP-chip panel reference space

## Fast commands

```bash
bash scripts/check_workflow_invariants.sh
bash -n sbatch/call_panel_variants_and_gtcheck.sbatch
bash -n sbatch/call_panel_variants_and_gtcheck_array.sbatch
bash tests/run_tiny_test.sh
bash tests/run_tiny_plot_test.sh
bash tests/run_tiny_pca_test.sh
bash tests/run_tiny_grin_test.sh
```

## Good review language

- Be specific about which file and assumption are at risk.
- Say clearly when something is a documentation drift issue versus a workflow bug.
- If there are no material issues, say that explicitly and note any remaining test gap.

