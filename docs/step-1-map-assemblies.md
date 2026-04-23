# Step 1: Map Assemblies to the Reference

## Submit the Example Job

```bash
sbatch sbatch/map_assemblies_to_reference.sbatch
```

## Variables to Review

The loop script is meant to be edited. The most important variables near the top are:

```text
REFERENCE_FASTA=reference/genome.fa
ASSEMBLY_DIR=assemblies
BAM_DIR=bams
PRESET=asm20
THREADS=16
```

Also review the `#SBATCH` resource lines and module-loading block for your cluster.

## What It Produces

For each assembly, the script writes:

```text
bams/<sample>.bam
bams/<sample>.bam.bai
```

## Core Logic

The central command is:

```bash
minimap2 -ax "$PRESET" -t "$THREADS" "$REFERENCE_FASTA" "$assembly" \
  | samtools sort -@ "$THREADS" -o "$bam"
samtools index "$bam"
```

## Choosing the `minimap2` Preset

- `asm5`: closely related assemblies
- `asm10`: intermediate divergence
- `asm20`: more permissive, often a practical within-species starting point

For many crop identity-checking cases, `asm20` is a reasonable default because it helps recover marker positions across assemblies that are not perfectly reference-like.

## Array Version

For larger projects, use:

```text
sbatch/map_assemblies_to_reference_array.sbatch
```

Typical pattern:

```bash
find assemblies -maxdepth 1 \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) | sort > assemblies.fofn
N=$(wc -l < assemblies.fofn)
sbatch --array=1-"$N" sbatch/map_assemblies_to_reference_array.sbatch
```

See also [HPC Notes](hpc_notes.md).
