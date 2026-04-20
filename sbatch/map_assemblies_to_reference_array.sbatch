#!/bin/bash
#SBATCH --job-name=map_asm_array
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Submit after creating a manifest, for example:
#   find assemblies -maxdepth 1 \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) | sort > assemblies.fofn
#   N=$(wc -l < assemblies.fofn)
#   sbatch --array=1-"$N" sbatch/map_assemblies_to_reference_array.sbatch

if command -v module >/dev/null 2>&1; then
    module purge || true
    module load minimap2 || module load minimap2/2.28-GCCcore-13.2.0 || true
    module load samtools || module load SAMtools/1.21-GCC-13.3.0 || true
fi

REFERENCE_FASTA="${REFERENCE_FASTA:-reference/genome.fa}"
ASSEMBLY_MANIFEST="${ASSEMBLY_MANIFEST:-assemblies.fofn}"
BAM_DIR="${BAM_DIR:-bams}"
PRESET="${PRESET:-asm20}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p logs "$BAM_DIR"

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "ERROR: submit this script as a SLURM array job"; exit 1; }
[[ -s "$REFERENCE_FASTA" ]] || { echo "ERROR: missing reference FASTA: $REFERENCE_FASTA"; exit 1; }
[[ -s "$ASSEMBLY_MANIFEST" ]] || { echo "ERROR: missing assembly manifest: $ASSEMBLY_MANIFEST"; exit 1; }

assembly=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ASSEMBLY_MANIFEST")
[[ -n "$assembly" ]] || { echo "ERROR: no assembly listed for array task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
[[ -s "$assembly" ]] || { echo "ERROR: assembly FASTA does not exist or is empty: $assembly"; exit 1; }

sample=$(basename "$assembly")
sample=${sample%.fa}
sample=${sample%.fasta}
sample=${sample%.fna}
bam="$BAM_DIR/${sample}.${PRESET}.sorted.bam"

echo "Array job:       ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task:      $SLURM_ARRAY_TASK_ID"
echo "Reference:       $REFERENCE_FASTA"
echo "Assembly:        $assembly"
echo "Sample:          $sample"
echo "Output BAM:      $bam"
echo "minimap2 preset: $PRESET"
echo "Threads:         $THREADS"

if [[ -s "$bam" && -s "${bam}.bai" ]]; then
    echo "[$sample] BAM and index already exist, skipping"
    exit 0
fi

minimap2 -ax "$PRESET" -t "$THREADS" "$REFERENCE_FASTA" "$assembly" \
    | samtools sort -@ "$THREADS" -o "$bam"

samtools index "$bam"

echo "[$sample] Wrote $bam"
