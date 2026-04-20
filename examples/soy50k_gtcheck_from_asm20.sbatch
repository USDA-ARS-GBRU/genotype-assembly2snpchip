#!/bin/bash
#SBATCH -J soy50k_gtcheck
#SBATCH -p batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH -o 00_log/soy50k_gtcheck-%x_%j.out
#SBATCH -e 00_log/soy50k_gtcheck-%x_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

module purge
module load BCFtools/1.21-GCC-13.3.0

# -----------------------------
# User-configurable paths
# -----------------------------
refG="ref/Gmax_880_v6.0.fa"
bamDir="04b_bams-asm20"
panelVCF="50k/50k.renamed.vcf.gz"
outDir="05c_vcfs_soy50k"
panelDir="05c_vcfs_soy50k/panel"

mkdir -p 00_log "$outDir" "$panelDir"

echo "Job ID:            $SLURM_JOBID"
echo "Node:              ${SLURMD_NODENAME:-NA}"
echo "Submit dir:        $SLURM_SUBMIT_DIR"
echo "Reference genome:  $refG"
echo "BAM directory:     $bamDir"
echo "Panel VCF:         $panelVCF"
echo "Output directory:  $outDir"
echo

# -----------------------------
# Sanity checks
# -----------------------------
[[ -s "$refG" ]] || { echo "ERROR: Reference genome not found: $refG"; exit 1; }
[[ -s "$panelVCF" ]] || { echo "ERROR: Panel VCF not found: $panelVCF"; exit 1; }
[[ -s "${panelVCF}.tbi" || -s "${panelVCF}.csi" ]] || { echo "ERROR: Panel VCF index not found for: $panelVCF"; exit 1; }
[[ -d "$bamDir" ]] || { echo "ERROR: BAM directory not found: $bamDir"; exit 1; }

# -----------------------------
# Prepare targeted SoySNP50K panel
# -----------------------------
panelSNP="${panelDir}/soy50k.biallelic.snps.vcf.gz"
panelPOS="${panelDir}/soy50k.positions.tsv"
panelALLELES="${panelDir}/soy50k.alleles.tsv"

if [[ ! -s "$panelSNP" || ! -s "${panelSNP}.tbi" ]]; then
    echo "[panel] Building biallelic SNP-only SoySNP50K panel"
    bcftools view \
        -m2 -M2 \
        -v snps \
        -Oz \
        -o "$panelSNP" \
        "$panelVCF"

    bcftools index -f -t "$panelSNP"
fi

if [[ ! -s "$panelPOS" ]]; then
    echo "[panel] Writing target positions list"
    bcftools query -f '%CHROM\t%POS\n' "$panelSNP" > "$panelPOS"
fi

if [[ ! -s "$panelALLELES" ]]; then
    echo "[panel] Writing allele target file"
    bcftools query \
        -f '%CHROM\t%POS\t%REF,%ALT\n' \
        "$panelSNP" > "$panelALLELES"
fi

# -----------------------------
# Process each BAM
# -----------------------------
shopt -s nullglob
bamFiles=("$bamDir"/*.sorted.bam)

if [[ ${#bamFiles[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found in $bamDir"
    exit 1
fi

echo "Found ${#bamFiles[@]} BAM files"
printf '  %s\n' "${bamFiles[@]}"
echo

for bam in "${bamFiles[@]}"; do
    base=$(basename "$bam")
    sample=${base%.sorted.bam}
    sample=${sample%-asm20}

    bai1="${bam}.bai"
    bai2="${bam%.bam}.bai"
    if [[ ! -s "$bai1" && ! -s "$bai2" ]]; then
        echo "ERROR: BAM index not found for $bam"
        exit 1
    fi

    rawVCF="${outDir}/${sample}.soy50k.raw.vcf.gz"
    namedVCF="${outDir}/${sample}.soy50k.named.vcf.gz"
    filtVCF="${outDir}/${sample}.soy50k.filtered.vcf.gz"
    diploidVCF="${outDir}/${sample}.soy50k.filtered.diploid.vcf.gz"
    gtcheckTSV="${outDir}/${sample}.soy50k.filtered.gtcheck.tsv"
    sampleNameFile="${outDir}/${sample}.sample_name.txt"
    statsTXT="${outDir}/${sample}.soy50k.qc.txt"

    echo "=================================================="
    echo "Sample: $sample"
    echo "BAM:    $bam"
    echo "=================================================="

    if [[ ! -s "$rawVCF" || ! -s "${rawVCF}.tbi" ]]; then
        echo "[$sample] Calling SoySNP50K target sites"
        bcftools mpileup \
            -Ou \
            --threads "$SLURM_CPUS_PER_TASK" \
            -f "$refG" \
            -T "$panelSNP" \
            --no-BAQ \
            -Q 0 \
            -q 5 \
            --max-depth 100 \
            -a FORMAT/AD,FORMAT/DP \
            "$bam" \
        | bcftools call \
            -Ou \
            -m \
            -C alleles \
            -T "$panelALLELES" \
            -i \
            -V indels \
            -a GQ \
        | bcftools norm \
            -Ou \
            -f "$refG" \
        | bcftools annotate \
            -x INFO,^FORMAT/GT,FORMAT/DP,FORMAT/AD,FORMAT/GQ \
            -Oz \
            -o "$rawVCF"

        bcftools index -f -t "$rawVCF"
    else
        echo "[$sample] Raw targeted VCF already exists, skipping"
    fi

    if [[ ! -s "$namedVCF" || ! -s "${namedVCF}.tbi" ]]; then
        echo "[$sample] Reheadering sample name"
        printf '%s\n' "$sample" > "$sampleNameFile"
        bcftools reheader \
            -s "$sampleNameFile" \
            -o "$namedVCF" \
            "$rawVCF"
        bcftools index -f -t "$namedVCF"
    else
        echo "[$sample] Named VCF already exists, skipping"
    fi

    if [[ ! -s "$filtVCF" || ! -s "${filtVCF}.tbi" ]]; then
        echo "[$sample] Setting weak genotype calls to missing"
        bcftools filter \
            -S . \
            -e 'FMT/DP<1 || FMT/GQ<20' \
            -Oz \
            -o "$filtVCF" \
            "$namedVCF"
        bcftools index -f -t "$filtVCF"
    else
        echo "[$sample] Filtered VCF already exists, skipping"
    fi

    if [[ ! -s "$diploidVCF" || ! -s "${diploidVCF}.tbi" ]]; then
        echo "[$sample] Forcing query GT fields to diploid before gtcheck"
        bcftools +fixploidy \
            "$filtVCF" \
            -Oz \
            -o "$diploidVCF" \
            -- \
            -f 2
        bcftools index -f -t "$diploidVCF"
    else
        echo "[$sample] Diploid-normalized VCF already exists, skipping"
    fi

    if [[ ! -s "$gtcheckTSV" ]]; then
        echo "[$sample] Running gtcheck in GT-vs-GT mode, keeping reference-only sites"
        bcftools gtcheck \
            -u GT,GT \
            -E 0 \
            --keep-refs \
            -O t \
            -o "$gtcheckTSV" \
            -g "$panelSNP" \
            "$diploidVCF"
    else
        echo "[$sample] gtcheck output already exists, skipping"
    fi

    echo "[$sample] Writing simple QC summary"
    {
        echo -e "sample\t$sample"
        echo -ne "records\t"
        bcftools view -H "$diploidVCF" | wc -l
        echo -ne "called\t"
        bcftools query -f '[%GT\n]' "$diploidVCF" | awk '$1!="./."{n++} END{print n+0}'
        echo -ne "missing\t"
        bcftools query -f '[%GT\n]' "$diploidVCF" | awk '$1=="./."{n++} END{print n+0}'
        echo -ne "hom\t"
        bcftools query -f '[%GT\n]' "$diploidVCF" | awk '$1=="0/0" || $1=="1/1" || $1=="2/2" || $1=="0|0" || $1=="1|1" || $1=="2|2"{n++} END{print n+0}'
        echo -ne "het\t"
        bcftools query -f '[%GT\n]' "$diploidVCF" | awk '$1=="0/1" || $1=="1/0" || $1=="0|1" || $1=="1|0" || $1=="1/2" || $1=="2/1" || $1=="1|2" || $1=="2|1"{n++} END{print n+0}'
    } > "$statsTXT"

    echo "[$sample] Done"
    echo

done

echo "All samples finished"
echo "Job ID: $SLURM_JOBID"
echo "-- End of script --"
