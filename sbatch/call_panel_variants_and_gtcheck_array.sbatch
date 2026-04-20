#!/bin/bash
#SBATCH --job-name=panel_gtcheck_array
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Submit after creating a BAM manifest, for example:
#   find bams -maxdepth 1 -name '*.bam' | sort > bams.fofn
#   N=$(wc -l < bams.fofn)
#   sbatch --array=1-"$N" sbatch/call_panel_variants_and_gtcheck_array.sbatch

if command -v module >/dev/null 2>&1; then
    module purge || true
    module load bcftools || module load BCFtools/1.21-GCC-13.3.0 || true
    module load htslib || module load HTSlib || true
fi

REFERENCE_FASTA="${REFERENCE_FASTA:-reference/genome.fa}"
PANEL_VCF="${PANEL_VCF:-data/snp_chip_panel.vcf.gz}"
BAM_MANIFEST="${BAM_MANIFEST:-bams.fofn}"
RESULTS_DIR="${RESULTS_DIR:-results}"
PANEL_DIR="${PANEL_DIR:-work/panel}"
MIN_MQ="${MIN_MQ:-5}"
MIN_GQ="${MIN_GQ:-20}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
VCF_SUFFIX="${VCF_SUFFIX:-panel}"

mkdir -p logs "$RESULTS_DIR" "$PANEL_DIR"

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "ERROR: submit this script as a SLURM array job"; exit 1; }
[[ -s "$REFERENCE_FASTA" ]] || { echo "ERROR: missing reference FASTA: $REFERENCE_FASTA"; exit 1; }
[[ -s "$PANEL_VCF" ]] || { echo "ERROR: missing panel VCF: $PANEL_VCF"; exit 1; }
[[ -s "${PANEL_VCF}.tbi" || -s "${PANEL_VCF}.csi" ]] || { echo "ERROR: missing panel VCF index"; exit 1; }
[[ -s "$BAM_MANIFEST" ]] || { echo "ERROR: missing BAM manifest: $BAM_MANIFEST"; exit 1; }

panel_snps="$PANEL_DIR/panel.biallelic.snps.vcf.gz"
panel_alleles="$PANEL_DIR/panel.alleles.tsv.gz"
panel_positions="$PANEL_DIR/panel.positions.tsv"
panel_lock="$PANEL_DIR/.panel_build.lock"

prepare_panel() {
    if [[ -s "$panel_snps" && -s "${panel_snps}.tbi" && -s "$panel_alleles" && -s "${panel_alleles}.tbi" && -s "$panel_positions" ]]; then
        return
    fi

    if mkdir "$panel_lock" 2>/dev/null; then
        trap 'rmdir "$panel_lock" 2>/dev/null || true' EXIT
        echo "[panel] Preparing reusable biallelic SNP panel files"

        if [[ ! -s "$panel_snps" || ! -s "${panel_snps}.tbi" ]]; then
            bcftools view -m2 -M2 -v snps -Oz -o "$panel_snps" "$PANEL_VCF"
            bcftools index -f -t "$panel_snps"
        fi

        if [[ ! -s "$panel_positions" ]]; then
            bcftools query -f '%CHROM\t%POS\n' "$panel_snps" > "$panel_positions"
        fi

        if [[ ! -s "$panel_alleles" || ! -s "${panel_alleles}.tbi" ]]; then
            bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' "$panel_snps" | bgzip -c > "$panel_alleles"
            tabix -f -s1 -b2 -e2 "$panel_alleles"
        fi

        rmdir "$panel_lock"
        trap - EXIT
    else
        echo "[panel] Another array task is preparing panel files; waiting"
        while [[ -d "$panel_lock" ]]; do
            sleep 10
        done
    fi

    [[ -s "$panel_snps" && -s "${panel_snps}.tbi" ]] || { echo "ERROR: panel SNP VCF was not created"; exit 1; }
    [[ -s "$panel_alleles" && -s "${panel_alleles}.tbi" ]] || { echo "ERROR: panel allele target file was not created"; exit 1; }
}

prepare_panel

bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_MANIFEST")
[[ -n "$bam" ]] || { echo "ERROR: no BAM listed for array task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
[[ -s "$bam" ]] || { echo "ERROR: BAM does not exist or is empty: $bam"; exit 1; }

if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
    echo "ERROR: missing BAM index for $bam"
    exit 1
fi

sample=$(basename "$bam")
sample=${sample%.sorted.bam}
sample=${sample%.bam}
sample=${sample%.asm20}
sample=${sample%.asm10}
sample=${sample%.asm5}

raw_vcf="$RESULTS_DIR/${sample}.${VCF_SUFFIX}.raw.vcf.gz"
named_vcf="$RESULTS_DIR/${sample}.${VCF_SUFFIX}.named.vcf.gz"
filtered_vcf="$RESULTS_DIR/${sample}.${VCF_SUFFIX}.filtered.vcf.gz"
gtcheck_tsv="$RESULTS_DIR/${sample}.gtcheck.tsv"
sample_name_file="$RESULTS_DIR/${sample}.sample_name.txt"
qc_tsv="$RESULTS_DIR/${sample}.${VCF_SUFFIX}.qc.tsv"

echo "Array job:     ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task:    $SLURM_ARRAY_TASK_ID"
echo "Reference:     $REFERENCE_FASTA"
echo "Panel VCF:     $PANEL_VCF"
echo "BAM:           $bam"
echo "Sample:        $sample"
echo "Results dir:   $RESULTS_DIR"
echo "MIN_MQ:        $MIN_MQ"
echo "MIN_GQ:        $MIN_GQ"

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

if [[ ! -s "$gtcheck_tsv" ]]; then
    bcftools gtcheck \
        -u GT,GT \
        -E 0 \
        -O t \
        -o "$gtcheck_tsv" \
        -g "$panel_snps" \
        "$filtered_vcf"
fi

{
    echo -e "metric\tvalue"
    echo -ne "records\t"
    bcftools view -H "$filtered_vcf" | wc -l
    echo -ne "called\t"
    bcftools query -f '[%GT\n]' "$filtered_vcf" | awk '$1!="./." && $1!="."{n++} END{print n+0}'
    echo -ne "missing\t"
    bcftools query -f '[%GT\n]' "$filtered_vcf" | awk '$1=="./." || $1=="."{n++} END{print n+0}'
    echo -ne "heterozygous\t"
    bcftools query -f '[%GT\n]' "$filtered_vcf" | awk '$1=="0/1" || $1=="1/0" || $1=="0|1" || $1=="1|0"{n++} END{print n+0}'
} > "$qc_tsv"

echo "[$sample] Wrote $gtcheck_tsv"
