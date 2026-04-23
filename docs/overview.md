# Overview

## The Question

The biological question is simple:

> Does this assembly look like the accession, line, cultivar, or individual that we think it is?

This comes up often in plant genomics when assemblies are renamed, transferred between projects, compared against public SNP-chip databases, or derived from seed lots with uncertain labels.

This repository is an identity-checking workflow, not a whole-genome variant-discovery workflow.

## Why Not Use a Sparse Variant VCF?

A standard variant-calling VCF usually records sites where the query differs from the reference. That is often the wrong data shape for SNP-chip identity checking.

If the query assembly is reference-like, many informative SNP-chip sites may be omitted from a variant-only VCF. For identity matching, omitted is not the same thing as confidently reference.

The key design choice in this workflow is:

> Genotype the assembly at the exact SNP-chip marker positions.

That preserves reference-state calls when they are supported and makes the assembly directly comparable to the panel.

## Workflow Overview

```text
query assembly FASTA
  -> minimap2 assembly-to-reference alignment
  -> sorted and indexed BAM
  -> bcftools mpileup/call restricted to SNP-chip marker sites
  -> filtered single-sample query VCF
  -> bcftools gtcheck against the panel VCF
  -> Python summary, metadata enrichment, and plots
```

## Critical Requirement

The reference genome used to map the whole-genome assembly must match the reference genome coordinate system used by the SNP-chip panel VCF.

At minimum, the mapping reference and panel VCF must agree on:

- chromosome or contig names
- marker positions
- REF alleles
- reference assembly version / coordinate system

Best practice is to map assemblies to the exact FASTA used to create or coordinate the panel VCF. If the panel was built on another reference assembly, lift over or rebuild the panel first.

## Related Pages

- [Setup and Environment](setup.md)
- [Inputs](inputs.md)
- [Step 1: Map Assemblies to the Reference](step-1-map-assemblies.md)
