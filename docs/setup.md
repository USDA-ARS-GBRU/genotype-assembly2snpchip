# Setup and Environment

## Repository Layout

```text
.
├── scripts/   species-agnostic Python utilities
├── sbatch/    example SLURM scripts for mapping and panel genotyping
├── docs/      workflow documentation
├── examples/  soybean-specific case-study material
└── tests/     tiny regression fixtures and checks
```

## Suggested Working Directories

```text
reference/      reference genome FASTA and indexes
assemblies/     query assembly FASTA files
data/           SNP-chip or marker-panel VCF
bams/           assembly-to-reference BAM files
work/           reusable panel intermediates
results/        query VCFs, gtcheck outputs, and summaries
logs/           SLURM stdout/stderr logs
```

Large FASTA, BAM, and VCF files are intentionally ignored by git.

## Software Requirements

On an HPC cluster, load or install:

- `minimap2`
- `samtools`
- `bcftools`
- `htslib`, including `bgzip` and `tabix`
- Python 3
- `pandas`, `matplotlib`, `seaborn`, `scikit-learn`, and `openpyxl`

## Environment Options

### Mamba

```bash
mamba env create -f environment.yml
mamba activate genotype-assembly2snpchip
```

### Conda

```bash
conda env create -f environment.yml
conda activate genotype-assembly2snpchip
```

### Pixi

```bash
pixi install
pixi shell
```

### pip

```bash
python -m pip install -r requirements.txt
```

## Quick Verification

```bash
minimap2 --version
samtools --version
bcftools --version
python -c "import pandas, matplotlib, seaborn, sklearn, openpyxl; print('Python packages OK')"
```

## SLURM Parameters to Review First

Before submitting jobs, open the `sbatch/` scripts and adjust:

- `#SBATCH --partition`
- `#SBATCH --account` if your cluster requires one
- `#SBATCH --cpus-per-task`
- `#SBATCH --mem`
- `#SBATCH --time`
- module names or environment activation blocks

For more cluster-specific guidance, see [HPC Notes](hpc_notes.md).
