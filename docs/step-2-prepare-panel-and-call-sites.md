# Step 2: Prepare the Panel and Call Panel Sites

## Submit the Example Job

```bash
sbatch sbatch/call_panel_variants_and_gtcheck.sbatch
```

## Variables to Review

The most important workflow variables in this script are:

```text
REFERENCE_FASTA=reference/genome.fa
PANEL_VCF=data/snp_chip_panel.vcf.gz
BAM_DIR=bams
RESULTS_DIR=results
PANEL_DIR=work/panel
MIN_MQ=5
MIN_GQ=20
THREADS=16
```

Also update the `#SBATCH` lines and the cluster-specific module or environment block.

## What This Step Produces

Reusable panel files:

```text
work/panel/panel.biallelic.snps.vcf.gz
work/panel/panel.positions.tsv
work/panel/panel.alleles.tsv.gz
```

Per-sample query VCFs:

```text
results/<sample>.panel.raw.vcf.gz
results/<sample>.panel.filtered.vcf.gz
results/<sample>.panel.filtered.diploid.vcf.gz
```

Per-sample comparison output:

```text
results/<sample>*.gtcheck.tsv
```

## What the Script Is Doing

### 1. Prepare a reusable SNP-only panel

The workflow filters the panel to biallelic SNPs and extracts known alleles for `bcftools call -C alleles`.

### 2. Genotype only the panel sites

This is the central logic:

```bash
bcftools mpileup -f "$REFERENCE_FASTA" -T "$panel_snps" "$bam" \
| bcftools call -m -C alleles -T "$panel_alleles" -i 1 -V indels -a GQ
```

The important idea is that the workflow is not doing generic variant discovery. It is asking for genotypes at the predefined SNP-chip markers.

### 3. Set weak calls to missing

The script keeps marker records but converts weak genotypes to missing:

```bash
bcftools filter -S . -e "FMT/DP<1 || FMT/GQ<${MIN_GQ}"
```

That keeps the panel structure intact and makes missingness easier to interpret.

### 4. Force diploid GT values before `gtcheck`

The workflow runs:

```bash
bcftools +fixploidy -- -f 2
```

That prevents `gtcheck` from skipping useful rows because only diploid `FORMAT/GT` fields are supported.

## Array Version

For larger projects, use:

```text
sbatch/call_panel_variants_and_gtcheck_array.sbatch
```

Typical pattern:

```bash
find bams -maxdepth 1 -name '*.bam' | sort > bams.fofn
N=$(wc -l < bams.fofn)
sbatch --array=1-"$N" sbatch/call_panel_variants_and_gtcheck_array.sbatch
```
