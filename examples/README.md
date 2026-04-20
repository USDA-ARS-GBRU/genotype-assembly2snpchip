Here is a clean, GitHub-ready README.md that pulls everything together into a coherent, instructive workflow.

You can drop this directly into your repo.

⸻

SoySNP50K Assembly Identity Validation Pipeline

Overview

This repository contains a reproducible workflow for validating the identity of whole-genome soybean assemblies using the SoySNP50K genotype panel (~20,000 accessions).

The core question this pipeline answers:

Does this assembled genome match the accession label we think it is?

⸻

Conceptual Foundation

The Problem

Whole-genome assemblies aligned to a reference and converted to VCF:
    •   only report variant sites
    •   differ across samples
    •   perform poorly for near-reference genomes

This leads to a major issue:

Identity comparison tools (like bcftools gtcheck) often have very few overlapping sites to compare — sometimes zero.

⸻

The Solution

Instead of comparing arbitrary variant sites:

Genotype every assembly at the exact SoySNP50K marker positions

This creates a standardized genotype profile across all samples.

⸻

Workflow Summary

Assembly FASTA
    ↓
minimap2 (assembly → reference)
    ↓
BAM
    ↓
bcftools mpileup + call (targeted to SoySNP50K sites)
    ↓
SoySNP50K-only VCF (query)
    ↓
bcftools gtcheck
    ↓
Identity comparison vs 20k accessions


⸻

Directory Structure

00_log/                 # SLURM logs
01_sbatch/              # SLURM scripts
02_scripts/             # helper scripts (optional)
03_genomes/             # reference genome (Wm82a6)
04b_bams-asm20/         # assembly → genome BAMs
05_vcfs_soy50k/         # output VCFs + gtcheck results
50/                     # SoySNP50K VCF


⸻

Step 1: Align Assemblies to Reference

Use minimap2 to align assemblies to Wm82a6.

Script: minimap-A2G.sbatch

#!/bin/bash
#SBATCH --job-name=minimap_A2G
#SBATCH --output=00_log/minimap_A2G_%j.out
#SBATCH --error=00_log/minimap_A2G_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

module load minimap2
module load samtools

ref="03_genomes/Wm82a6.fa"
indir="03_genomes/assemblies"
outdir="04b_bams-asm20"
mkdir -p ${outdir}

for asm in ${indir}/*.fa; do
    sample=$(basename ${asm} .fa)

    minimap2 -ax asm20 -t 8 ${ref} ${asm} \
    | samtools sort -@ 4 -o ${outdir}/${sample}-asm20.sorted.bam

    samtools index ${outdir}/${sample}-asm20.sorted.bam
done


⸻

Step 2: Prepare SoySNP50K Panel

bcftools view \
    -m2 -M2 \
    -v snps \
    -Oz \
    -o 50k_targeted/soy50k.biallelic.snps.vcf.gz \
    50/50k.renamed.vcf.gz

bcftools index -t 50k_targeted/soy50k.biallelic.snps.vcf.gz


⸻

Step 3: Genotype Assemblies at SoySNP50K Sites

Script: soy50k_gtcheck_from_asm20.sbatch

Key operations:
    •   Restrict genotyping to SoySNP50K sites
    •   Call genotypes using known alleles
    •   Retain reference matches (critical)
    •   Set weak calls to missing (not removed)

⸻

Important Parameter Choices

bcftools mpileup

-q 5     # filter low-confidence mappings
-Q 0     # keep all bases (assembly BAMs don't use BQ meaningfully)
--no-BAQ # avoid BAQ adjustment

Interpretation:
    •   -q: removes ambiguous contig placements
    •   -Q: irrelevant for assemblies → keep permissive

⸻

bcftools call

-C alleles -T panel.vcf.gz

→ genotype only known SoySNP50K alleles

-i 1

→ include all target sites (even weak ones)

⸻

bcftools gtcheck

-u GT,GT -E 0 --keep-refs

→ ensures:
    •   genotype vs genotype comparison
    •   discordance = actual mismatch count
    •   reference-only / monoallelic sites are retained in the comparison

Before gtcheck, run bcftools +fixploidy with -- -f 2 on the filtered query VCF so haploid GT fields are converted to diploid GT fields. This prevents gtcheck from skipping sites with the alert "only diploid FORMAT/GT fields supported."

⸻

Output Files

Per sample:

sample.soy50k.raw.vcf.gz
sample.soy50k.named.vcf.gz
sample.soy50k.filtered.vcf.gz
sample.soy50k.filtered.diploid.vcf.gz
sample.soy50k.filtered.gtcheck.tsv
sample.soy50k.qc.txt


⸻

Understanding gtcheck Output

Each row compares:

query_sample  panel_sample  discordance  HWE_score  sites_compared  matching_genotypes

Key columns

sites_compared
Number of SoySNP50K markers used in comparison
    •   high → strong evidence
    •   low → weak / unreliable

matching_genotypes
Number of sites where genotypes match

match_fraction
(custom metric)

matching_genotypes / sites_compared

Most intuitive identity measure.

⸻

Interpreting Results

Strong Identity Match
    •   match_fraction ≥ 0.98
    •   sites_compared ≥ several thousand
    •   clear separation from next best hit

Example:

Rank 1: 0.995
Rank 2: 0.960


⸻

Ambiguous / Related Accessions

Rank 1: 0.992
Rank 2: 0.989
Rank 3: 0.988

→ likely closely related lines or shared ancestry

⸻

Weak / Problematic Sample
    •   low match_fraction
    •   few sites_compared
    •   no clear best hit

Possible causes:
    •   poor mapping
    •   low coverage
    •   contamination
    •   mislabeled sample

⸻

Special Cases

Wm82a6 mapped to itself
    •   often yields 0 comparable sites in variant-only workflows
    •   NOT informative
    •   fixed by targeted genotyping

⸻

Wm82v4 vs Wm82a6
    •   very few differences
    •   variant-only workflows underrepresent overlap
    •   targeted workflow resolves this

⸻

Quality Control Checks

Call rate

called vs missing genotypes

Heterozygosity

Soybean is highly inbred:
    •   high heterozygosity → warning sign

⸻

Detecting Sample Mix-ups

Look for:
    •   top hit ≠ expected accession
    •   strong match to a different accession
    •   consistent mismatch pattern

⸻

Best Practices
    •   Always use targeted marker genotyping
    •   Avoid variant-only query VCFs
    •   Interpret results using:
    •   match_fraction
    •   sites_compared
    •   rank separation
    •   Validate pipeline using known controls (Wm82)

⸻

Key Insight

Identity validation depends more on consistency across many known markers than on discovering new variants.

⸻

Future Extensions
    •   automated mismatch flagging
    •   PCA validation plots
    •   clustering of top hits
    •   report generation (PDF summaries)

⸻

Bottom Line

This workflow transforms genome assemblies into panel-compatible genotype profiles, enabling robust and interpretable identity validation against a large soybean reference population.
