# Genotype Assembly2SNPchip Agent

You are operating a bioinformatics workflow repository for validating assembly identity against SNP-chip or marker-panel genotypes.

## Mission

Help maintain and extend the generic workflow while protecting its key scientific and operational assumptions.

## Default behavior

- Prefer the generic workflow over soybean-only conventions.
- Keep documentation synchronized with workflow behavior.
- Treat `sbatch/` as the canonical HPC execution surface.
- Treat `scripts/` as the canonical post-processing surface.
- Keep fixes conservative unless the user asks for a workflow redesign.

## Scientific guardrails

- Reference compatibility is pivotal: the mapping reference must match, or be coordinate-compatible with, the SNP-chip panel reference.
- `gtcheck` should be interpreted as an identity-comparison workflow, not whole-genome discovery.
- Missingness and low site overlap matter; do not oversell weak comparisons.
- Preserve support for any species and any SNP-chip-style panel VCF.

## Operational guardrails

- `bcftools >= 1.23` is required.
- `bcftools call` uses bare `-i`.
- `bcftools gtcheck` uses `--keep-refs`.
- `bcftools +fixploidy` runs before `gtcheck`.
- Update both loop and array sbatch scripts when changing the core panel genotyping workflow.

## Useful first reads

- `README.md`
- `docs/setup.md`
- `docs/step-2-prepare-panel-and-call-sites.md`
- `docs/step-3-compare-and-summarize.md`
- `scripts/`
- `.agents/memory.md`

