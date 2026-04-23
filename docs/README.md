# Genotype Assembly to SNP-Chip Panel Docs

Welcome to the documentation site for `genotype-assembly2snpchip`.

This guide is organized for two kinds of readers:

- a new graduate student running the workflow for the first time
- an experienced bioinformatician who wants the exact inputs, steps, and outputs quickly

The workflow answers a practical identity question:

> Does this genome assembly match the accession, line, cultivar, or individual we think it is?

It does that by mapping one or more assemblies to a reference genome, genotyping them at SNP-chip marker positions, comparing those genotypes to a multi-sample panel VCF with `bcftools gtcheck`, and summarizing the closest matches.

## Start Here

- [Overview](overview.md)
- [Setup and Environment](setup.md)
- [Inputs](inputs.md)

## Workflow

- [Step 1: Map Assemblies to the Reference](step-1-map-assemblies.md)
- [Step 2: Prepare the Panel and Call Panel Sites](step-2-prepare-panel-and-call-sites.md)
- [Step 3: Compare to the Panel and Summarize Hits](step-3-compare-and-summarize.md)
- [Visualization](visualization.md)

## Results and Support Pages

- [Interpreting Results](interpreting-results.md)
- [Testing](testing.md)
- [HPC Notes](hpc_notes.md)

## Suggested Reading Path

For a first pass through the workflow:

1. [Overview](overview.md)
2. [Setup and Environment](setup.md)
3. [Inputs](inputs.md)
4. [Step 1: Map Assemblies to the Reference](step-1-map-assemblies.md)
5. [Step 2: Prepare the Panel and Call Panel Sites](step-2-prepare-panel-and-call-sites.md)
6. [Step 3: Compare to the Panel and Summarize Hits](step-3-compare-and-summarize.md)
7. [Visualization](visualization.md)
8. [Interpreting Results](interpreting-results.md)

## Repository Front Door

The repository README stays intentionally short and practical:

- [Repository README](https://github.com/USDA-ARS-GBRU/genotype-assembly2snpchip#readme)

That page is the quick intro, install guide, and quick start.
