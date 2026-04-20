#!/bin/bash
#SBATCH --job-name=panel_gtcheck
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Portable module loading. Edit these lines for your HPC, or comment them out
# if bcftools, bgzip, and tabix are already available in your environment.
if command -v module >/dev/null 2>&1; then
    module purge || true
    module load bcftools || module load BCFtools/1.21-GCC-13.3.0 || true
    module load htslib || module load HTSlib || true
fi

REFERENCE_FASTA="${REFERENCE_FASTA:-reference/genome.fa}"
PANEL_VCF="${PANEL_VCF:-data/snp_chip_panel.vcf.gz}"
BAM_DIR="${BAM_DIR:-bams}"
RESULTS_DIR="${RESULTS_DIR:-results}"
PANEL_DIR="${PANEL_DIR:-work/panel}"
MIN_MQ="${MIN_MQ:-5}"
MIN_GQ="${MIN_GQ:-20}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
VCF_SUFFIX="${VCF_SUFFIX:-panel}"

mkdir -p logs "$RESULTS_DIR" "$PANEL_DIR"

echo "Reference:     $REFERENCE_FASTA"
echo "Panel VCF:     $PANEL_VCF"
echo "BAM dir:       $BAM_DIR"
echo "Results dir:   $RESULTS_DIR"
echo "Panel workdir: $PANEL_DIR"
echo "MIN_MQ:        $MIN_MQ"
echo "MIN_GQ:        $MIN_GQ"

[[ -s "$REFERENCE_FASTA" ]] || { echo "ERROR: missing reference FASTA: $REFERENCE_FASTA"; exit 1; }
[[ -s "$PANEL_VCF" ]] || { echo "ERROR: missing panel VCF: $PANEL_VCF"; exit 1; }
[[ -s "${PANEL_VCF}.tbi" || -s "${PANEL_VCF}.csi" ]] || { echo "ERROR: missing panel VCF index"; exit 1; }
[[ -d "$BAM_DIR" ]] || { echo "ERROR: missing BAM directory: $BAM_DIR"; exit 1; }

panel_snps="$PANEL_DIR/panel.biallelic.snps.vcf.gz"
panel_alleles="$PANEL_DIR/panel.alleles.tsv.gz"
panel_positions="$PANEL_DIR/panel.positions.tsv"

if [[ ! -s "$panel_snps" || ! -s "${panel_snps}.tbi" ]]; then
    echo "[panel] Keeping only biallelic SNP markers"
    bcftools view -m2 -M2 -v snps -Oz -o "$panel_snps" "$PANEL_VCF"
    bcftools index -f -t "$panel_snps"
fi

if [[ ! -s "$panel_positions" ]]; then
    bcftools query -f '%CHROM\t%POS\n' "$panel_snps" > "$panel_positions"
fi

if [[ ! -s "$panel_alleles" || ! -s "${panel_alleles}.tbi" ]]; then
    echo "[panel] Writing bcftools call allele target file"
    bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' "$panel_snps" | bgzip -c > "$panel_alleles"
    tabix -f -s1 -b2 -e2 "$panel_alleles"
fi

shopt -s nullglob
bams=("$BAM_DIR"/*.bam)
[[ ${#bams[@]} -gt 0 ]] || { echo "ERROR: no BAM files found in $BAM_DIR"; exit 1; }

for bam in "${bams[@]}"; do
    [[ "$bam" == *.sorted.bam || "$bam" == *.bam ]] || continue
    sample=$(basename "$bam")
    sample=${sample%.sorted.bam}
    sample=${sample%.bam}
    sample=${sample%.asm20}
    sample=${sample%.asm10}
    sample=${sample%.asm5}

    if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
        echo "ERROR: missing BAM index for $bam"
        exit 1
    fi

    raw_vcf="$RESULTS_DIR/${sample}.${VCF_SUFFIX}.raw.vcf.gz"
    named_vcf="$RESULTS_DIR/${sample}.${VCF_SUFFIX}.named.vcf.gz"
    filtered_vcf="$RESULTS_DIR/${sample}.${VCF_SUFFIX}.filtered.vcf.gz"
    diploid_vcf="$RESULTS_DIR/${sample}.${VCF_SUFFIX}.filtered.diploid.vcf.gz"
    gtcheck_tsv="$RESULTS_DIR/${sample}.gtcheck.tsv"
    sample_name_file="$RESULTS_DIR/${sample}.sample_name.txt"
    qc_tsv="$RESULTS_DIR/${sample}.${VCF_SUFFIX}.qc.tsv"

    echo "============================================================"
    echo "Sample: $sample"
    echo "BAM:    $bam"
    echo "============================================================"

    if [[ ! -s "$raw_vcf" || ! -s "${raw_vcf}.tbi" ]]; then
        bcftools mpileup \
            -Ou \
            --threads "$THREADS" \
            -f "$REFERENCE_FASTA" \
            -T "$panel_snps" \
            --no-BAQ \
            -Q 0 \
            -q "$MIN_MQ" \
            --max-depth 100 \
            -a FORMAT/AD,FORMAT/DP \
            "$bam" \
        | bcftools call \
            -Ou \
            -m \
            -C alleles \
            -T "$panel_alleles" \
            -i 1 \
            -V indels \
            -a GQ \
        | bcftools norm \
            -Ou \
            -f "$REFERENCE_FASTA" \
        | bcftools annotate \
            -x INFO,^FORMAT/GT,FORMAT/DP,FORMAT/AD,FORMAT/GQ \
            -Oz \
            -o "$raw_vcf"
        bcftools index -f -t "$raw_vcf"
    fi

    if [[ ! -s "$named_vcf" || ! -s "${named_vcf}.tbi" ]]; then
        printf '%s\n' "$sample" > "$sample_name_file"
        bcftools reheader -s "$sample_name_file" -o "$named_vcf" "$raw_vcf"
        bcftools index -f -t "$named_vcf"
    fi

    if [[ ! -s "$filtered_vcf" || ! -s "${filtered_vcf}.tbi" ]]; then
        bcftools filter \
            -S . \
            -e "FMT/DP<1 || FMT/GQ<${MIN_GQ}" \
            -Oz \
            -o "$filtered_vcf" \
            "$named_vcf"
        bcftools index -f -t "$filtered_vcf"
    fi

    if [[ ! -s "$diploid_vcf" || ! -s "${diploid_vcf}.tbi" ]]; then
        bcftools +fixploidy \
            "$filtered_vcf" \
            -Oz \
            -o "$diploid_vcf" \
            -- \
            -f 2
        bcftools index -f -t "$diploid_vcf"
    fi

    if [[ ! -s "$gtcheck_tsv" ]]; then
        bcftools gtcheck \
            -u GT,GT \
            -E 0 \
            --keep-refs \
            -O t \
            -o "$gtcheck_tsv" \
            -g "$panel_snps" \
            "$diploid_vcf"
    fi

    {
        echo -e "metric\tvalue"
        echo -ne "records\t"
        bcftools view -H "$diploid_vcf" | wc -l
        echo -ne "called\t"
        bcftools query -f '[%GT\n]' "$diploid_vcf" | awk '$1!="./." && $1!="."{n++} END{print n+0}'
        echo -ne "missing\t"
        bcftools query -f '[%GT\n]' "$diploid_vcf" | awk '$1=="./." || $1=="."{n++} END{print n+0}'
        echo -ne "heterozygous\t"
        bcftools query -f '[%GT\n]' "$diploid_vcf" | awk '$1=="0/1" || $1=="1/0" || $1=="0|1" || $1=="1|0"{n++} END{print n+0}'
    } > "$qc_tsv"

    echo "[$sample] Wrote $gtcheck_tsv"
done

echo "Finished targeted panel genotyping and gtcheck for ${#bams[@]} BAM files."
