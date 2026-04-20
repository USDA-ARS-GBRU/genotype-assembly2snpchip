#!/bin/bash
#SBATCH --job-name=map_assemblies
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Portable module loading. Edit these lines for your HPC, or comment them out
# if minimap2 and samtools are already available in your environment.
if command -v module >/dev/null 2>&1; then
    module purge || true
    module load minimap2 || module load minimap2/2.28-GCCcore-13.2.0 || true
    module load samtools || module load SAMtools/1.21-GCC-13.3.0 || true
fi

REFERENCE_FASTA="${REFERENCE_FASTA:-reference/genome.fa}"
ASSEMBLY_DIR="${ASSEMBLY_DIR:-assemblies}"
BAM_DIR="${BAM_DIR:-bams}"
PRESET="${PRESET:-asm20}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p logs "$BAM_DIR"

echo "Reference:       $REFERENCE_FASTA"
echo "Assembly dir:    $ASSEMBLY_DIR"
echo "Output BAM dir:  $BAM_DIR"
echo "minimap2 preset: $PRESET"
echo "Threads:         $THREADS"

[[ -s "$REFERENCE_FASTA" ]] || { echo "ERROR: missing reference FASTA: $REFERENCE_FASTA"; exit 1; }
[[ -d "$ASSEMBLY_DIR" ]] || { echo "ERROR: missing assembly directory: $ASSEMBLY_DIR"; exit 1; }

shopt -s nullglob
assemblies=("$ASSEMBLY_DIR"/*.fa "$ASSEMBLY_DIR"/*.fasta "$ASSEMBLY_DIR"/*.fna)
[[ ${#assemblies[@]} -gt 0 ]] || { echo "ERROR: no FASTA assemblies found in $ASSEMBLY_DIR"; exit 1; }

for assembly in "${assemblies[@]}"; do
    sample=$(basename "$assembly")
    sample=${sample%.fa}
    sample=${sample%.fasta}
    sample=${sample%.fna}
    bam="$BAM_DIR/${sample}.${PRESET}.sorted.bam"

    echo "[$sample] Mapping $assembly"
    minimap2 -ax "$PRESET" -t "$THREADS" "$REFERENCE_FASTA" "$assembly" \
        | samtools sort -@ "$THREADS" -o "$bam"
    samtools index "$bam"
    echo "[$sample] Wrote $bam"
done

echo "Finished mapping ${#assemblies[@]} assemblies."
