# HPC Notes

This workflow is written for SLURM-style HPC clusters, but every cluster has its own local habits. Treat the scripts in `sbatch/` as readable templates: keep the workflow logic, then adjust partitions, account names, time limits, memory, and module names for the system you are using.

## What Usually Changes

The most common edits are near the top of each sbatch file:

```bash
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
```

Some clusters also require an account or allocation:

```bash
#SBATCH --account=my_lab_or_project
```

Module names also vary:

```bash
module load minimap2
module load samtools
module load bcftools
module load htslib
```

If a cluster uses conda, mamba, Apptainer/Singularity, or preloaded software paths instead of modules, replace the module block with the local environment setup.

## Sapelo2-Style Notes

Sapelo2 commonly uses EasyBuild-style module names. Depending on what is installed and visible in your module tree, commands may look like:

```bash
module load minimap2/2.28-GCCcore-13.2.0
module load SAMtools/1.21-GCC-13.3.0
module load BCFtools/1.21-GCC-13.3.0
```

The older soybean example scripts in `examples/` use this style. If those exact versions are not available, search:

```bash
module spider minimap2
module spider SAMtools
module spider BCFtools
```

Sapelo2 partition, account, and quality-of-service settings differ by group and allocation. Ask your lab or RC documentation which `#SBATCH` lines are required.

## SciNet-Style Notes

SciNet systems can differ depending on which cluster and software stack you are using. Start by checking available modules:

```bash
module avail minimap2
module avail samtools
module avail bcftools
```

or:

```bash
module spider minimap2
module spider samtools
module spider bcftools
```

If the exact tools are not available as modules, a conda/mamba environment is often the simplest portable option:

```bash
mamba create -n snpchip-id -c bioconda -c conda-forge minimap2 samtools bcftools
conda activate snpchip-id
```

For production work, confirm whether your cluster prefers centrally installed modules over user-managed conda environments.

## HyperGator-Style Notes

On HyperGator, you may need allocation-specific account and partition settings. A script header may require lines like:

```bash
#SBATCH --account=my_allocation
#SBATCH --partition=hpg-default
```

The exact values are lab- and allocation-specific. Check the local documentation or a known working script from your group.

Search for tool modules:

```bash
module spider minimap2
module spider samtools
module spider bcftools
```

If modules are unavailable or version-conflicted, a conda/mamba environment with Bioconda packages is a reasonable fallback.

## Loop Jobs vs Array Jobs

The provided scripts are loop-style jobs because they are easy for a new student to read:

```text
one job starts
  -> loops over all assemblies or BAMs
  -> writes all outputs
```

For a small class exercise or a modest number of samples, that is fine.

For larger projects, SLURM arrays are usually better:

```text
one array job starts
  -> each task handles one sample
  -> each sample gets its own log
  -> failed samples can be rerun without rerunning everything
```

This repository includes optional array versions:

```text
sbatch/map_assemblies_to_reference_array.sbatch
sbatch/call_panel_variants_and_gtcheck_array.sbatch
```

To run mapping as an array job, first make a manifest:

```bash
find assemblies -maxdepth 1 \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) | sort > assemblies.fofn
```

Then each array task can read one line:

```bash
assembly=$(sed -n "${SLURM_ARRAY_TASK_ID}p" assemblies.fofn)
```

Submit with the number of lines in the manifest:

```bash
N=$(wc -l < assemblies.fofn)
sbatch --array=1-"$N" my_array_script.sbatch
```

For targeted SNP-chip genotyping and `gtcheck`, make a BAM manifest:

```bash
find bams -maxdepth 1 -name '*.bam' | sort > bams.fofn
N=$(wc -l < bams.fofn)
sbatch --array=1-"$N" sbatch/call_panel_variants_and_gtcheck_array.sbatch
```

The panel calling array script shares panel-preparation files across tasks and uses a small lock directory to prevent multiple tasks from building those files at the same time.

## Scratch Space

Assembly BAMs and VCF intermediates can be large. On most clusters:

- keep source FASTA/VCF files in project storage,
- write heavy temporary files to scratch,
- copy final summaries back to project storage,
- avoid writing large intermediates to home directories.

The example scripts default to simple local directories because that is easiest to teach. For a real HPC run, set paths explicitly:

```bash
sbatch \
  --export=ALL,REFERENCE_FASTA=/project/ref/genome.fa,ASSEMBLY_DIR=/scratch/$USER/assemblies,BAM_DIR=/scratch/$USER/bams \
  sbatch/map_assemblies_to_reference.sbatch
```

## Quick Preflight Checklist

Before submitting a full run, check:

- `samtools faidx reference/genome.fa` has been run.
- The SNP-chip VCF is bgzipped and indexed.
- Reference contig names match the SNP-chip VCF contig names.
- A single assembly maps successfully.
- A single BAM produces a filtered panel-site VCF.
- `bcftools gtcheck` reports nonzero sites compared.
- The Python summary script works on a tiny fixture or one real `gtcheck.tsv`.

That small preflight saves a lot of queue time.
