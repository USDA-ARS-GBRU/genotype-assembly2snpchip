# Inputs

## Required Files

You need three main inputs:

1. one or more query assembly FASTA files
2. a reference FASTA in the same coordinate system as the SNP-chip panel
3. a bgzipped and indexed multi-sample SNP-chip / marker-panel VCF

## Reference Genome FASTA

Use the same reference coordinate system as the SNP-chip VCF. In the safest case, this is the exact FASTA used when the panel VCF was generated.

Required files:

```text
reference/genome.fa
reference/genome.fa.fai
```

Create the FASTA index if needed:

```bash
samtools faidx reference/genome.fa
```

## Query Assembly FASTA Files

Put one or more assemblies in:

```text
assemblies/
```

Typical accepted extensions in the example scripts:

```text
.fa
.fasta
.fna
```

Each assembly filename becomes the sample basename in downstream outputs.

## SNP-Chip or Marker-Panel VCF

The panel VCF should be compressed and indexed:

```text
data/snp_chip_panel.vcf.gz
data/snp_chip_panel.vcf.gz.tbi
```

If needed:

```bash
bgzip data/snp_chip_panel.vcf
tabix -p vcf data/snp_chip_panel.vcf.gz
```

The workflow filters this panel to biallelic SNPs before genotyping and comparison.

## Practical Compatibility Checks

Before a full run, verify:

- contig names match between the FASTA and VCF
- panel REF alleles agree with the reference FASTA
- the panel is on the same reference assembly version as the mapping FASTA

If those assumptions fail, downstream discordance can look biologically meaningful while actually being a coordinate-system problem.
