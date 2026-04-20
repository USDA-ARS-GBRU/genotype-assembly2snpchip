# Genotype Assembly to SNP-Chip Panel

This repository teaches a generic workflow for asking whether a genome assembly matches a known genotype in a SNP-chip database.

The same idea can be used for soybean, cotton, maize, wheat, sorghum, peanut, or any other project where you have:

1. a genome assembly FASTA for one or more query samples,
2. a reference genome FASTA for the same species or a close reference,
3. a SNP-chip or marker-panel VCF containing many known samples,
4. access to an HPC cluster with SLURM, `minimap2`, `samtools`, and `bcftools`.

The workflow was motivated by soybean/SoySNP50K work on Sapelo2, but the top-level code and documentation are intentionally species-agnostic. The soybean-specific material is kept in `examples/` as a concrete case study.

## The Question

The biological question is simple:

> Does this assembly look like the accession, line, cultivar, or individual that we think it is?

That question comes up often in plant genomics. Assemblies may be renamed, transferred between projects, derived from seed lots with confusing labels, or compared against public SNP-chip databases. A SNP-chip panel is useful because it contains standardized genotypes for many known samples. If we can genotype the assembly at those same marker positions, we can compare the assembly to every panel sample and identify the closest match.

This workflow is not primarily a variant-discovery workflow. It is an identity-checking workflow.

## Why Not Just Call Whole-Genome Variants?

A common first attempt is:

```text
assembly FASTA
  -> align assembly to reference
  -> call whole-genome variants
  -> compare that VCF to the SNP-chip panel
```

That often gives confusing results, especially when the assembly is very similar to the reference genome.

The reason is subtle but important. A normal variant-calling VCF usually records only sites where the query differs from the reference. If the query assembly is reference-like, the VCF may be very sparse. Many SNP-chip markers where the query has the reference allele are absent from the VCF, even though those reference-state genotypes are highly informative for identity checking.

For example, imagine a SNP-chip has 50,000 markers. At a particular marker, the panel has many samples with `A/A` and many with `G/G`. If your assembly matches the reference allele `A`, a variant-only VCF may omit that site entirely. To an identity comparison tool, omitted is not the same thing as confidently `A/A`; it is often treated as missing or unavailable.

The solution is to reverse the logic:

> Do not ask which variants happened to appear in the assembly. Instead, genotype the assembly at the exact SNP-chip marker positions.

That produces a standardized marker profile for each assembly.

## Workflow Overview

```text
query assembly FASTA
  -> minimap2 assembly-to-reference alignment
  -> sorted and indexed BAM
  -> bcftools mpileup/call restricted to SNP-chip marker sites
  -> single-sample query VCF containing panel-site genotypes
  -> bcftools gtcheck against the panel VCF
  -> Python summary of top panel matches
```

The key idea is that the query VCF should include genotypes at predefined panel positions, including reference matches when they can be called. Then `bcftools gtcheck` can compare the assembly-derived genotype profile to the SNP-chip database.

## Repository Layout

```text
.
├── README.md
├── requirements.txt
├── LICENSE
├── docs/
│   └── hpc_notes.md
├── sbatch/
│   ├── map_assemblies_to_reference.sbatch
│   ├── map_assemblies_to_reference_array.sbatch
│   ├── call_panel_variants_and_gtcheck.sbatch
│   └── call_panel_variants_and_gtcheck_array.sbatch
├── scripts/
│   └── summarize_gtcheck_top_hits.py
├── tests/
│   ├── fixtures/tiny.gtcheck.tsv
│   └── run_tiny_test.sh
└── examples/
    ├── README_with_qc_and_sitecount_notes.md
    ├── minimap-A2G.sbatch
    ├── soy50k_gtcheck_from_asm20.sbatch
    └── results/
        └── *.gtcheck.tsv
```

Suggested working directories for a new analysis:

```text
reference/      reference genome FASTA and indexes
assemblies/     query assembly FASTA files
data/           SNP-chip or marker-panel VCF
bams/           assembly-to-reference BAM files
work/           reusable panel intermediates
results/        query VCFs, gtcheck outputs, and summaries
logs/           SLURM stdout/stderr logs
```

Large FASTA, BAM, and VCF files are ignored by git. Keep the workflow code in the repository and keep large data files in project storage or scratch space.

## Software Requirements

On an HPC cluster, load or install:

- `minimap2`
- `samtools`
- `bcftools`
- `htslib`, including `bgzip` and `tabix`
- Python 3

The included Python script uses only the Python standard library. The `requirements.txt` file is intentionally empty except for a note, so this command is optional:

```bash
python -m pip install -r requirements.txt
```

Most HPC systems use environment modules, but module names differ. The sbatch scripts include generic module-loading lines plus examples from Sapelo2-style names. Edit those lines for your cluster.

## Inputs

### Reference Genome FASTA

Use the same reference coordinate system as the SNP-chip VCF.

This point matters. `bcftools` will only compare records cleanly when chromosome names, positions, and alleles are compatible. If the panel VCF uses `chr01` but your reference uses `1`, you must rename one of them before running the workflow.

Required files:

```text
reference/genome.fa
reference/genome.fa.fai
```

Create the FASTA index if needed:

```bash
samtools faidx reference/genome.fa
```

### Query Assembly FASTA Files

Put one or more assemblies in:

```text
assemblies/
```

Accepted extensions in the example sbatch script:

```text
.fa
.fasta
.fna
```

Each assembly filename becomes the sample name. For example:

```text
assemblies/Cultivar_A.fa
```

becomes sample:

```text
Cultivar_A
```

### SNP-Chip or Marker-Panel VCF

The panel VCF should contain genotypes for the known database samples. It should be compressed and indexed:

```text
data/snp_chip_panel.vcf.gz
data/snp_chip_panel.vcf.gz.tbi
```

If needed:

```bash
bgzip data/snp_chip_panel.vcf
tabix -p vcf data/snp_chip_panel.vcf.gz
```

The workflow filters this panel to biallelic SNPs because simple biallelic SNPs are easier to genotype and compare robustly across many samples.

## Step 1. Map Assemblies to the Reference

Use:

```bash
sbatch sbatch/map_assemblies_to_reference.sbatch
```

By default, the script expects:

```text
REFERENCE_FASTA=reference/genome.fa
ASSEMBLY_DIR=assemblies
BAM_DIR=bams
PRESET=asm20
```

You can override those values at submission time:

```bash
sbatch \
  --export=ALL,REFERENCE_FASTA=/path/to/ref.fa,ASSEMBLY_DIR=/path/to/assemblies,BAM_DIR=/path/to/bams,PRESET=asm20 \
  sbatch/map_assemblies_to_reference.sbatch
```

### What the Mapping Step Does

The core command is:

```bash
minimap2 -ax "$PRESET" -t "$THREADS" "$REFERENCE_FASTA" "$assembly" \
  | samtools sort -@ "$THREADS" -o "$bam"
samtools index "$bam"
```

`minimap2` aligns the query assembly contigs or scaffolds to the reference genome. `samtools sort` produces a coordinate-sorted BAM, and `samtools index` creates the index needed for fast targeted access.

### Choosing `asm5`, `asm10`, or `asm20`

`minimap2` has presets for assembly-to-reference alignment:

- `asm5`: closely related assemblies, roughly up to 5 percent divergence
- `asm10`: intermediate divergence
- `asm20`: more divergent assemblies, roughly up to 20 percent divergence

For many within-species crop genome comparisons, `asm20` is a forgiving default. If your assemblies are very close to the reference and you want stricter placements, try `asm5` or `asm10`.

For identity checking, a slightly permissive alignment can be useful because you want to recover as many panel marker positions as possible. However, overly permissive mapping can create ambiguous placements in repetitive or duplicated regions. The later `bcftools mpileup -q` setting helps filter low-confidence alignments.

### Optional Array-Job Mapping

For larger projects, use the array version so each assembly gets its own SLURM task and log file.

Create a manifest:

```bash
find assemblies -maxdepth 1 \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) | sort > assemblies.fofn
```

Submit one task per assembly:

```bash
N=$(wc -l < assemblies.fofn)
sbatch --array=1-"$N" sbatch/map_assemblies_to_reference_array.sbatch
```

The array script uses the same default variables as the loop script, plus:

```text
ASSEMBLY_MANIFEST=assemblies.fofn
```

You can still override paths at submission time:

```bash
sbatch \
  --array=1-"$N" \
  --export=ALL,REFERENCE_FASTA=/path/to/ref.fa,ASSEMBLY_MANIFEST=/path/to/assemblies.fofn,BAM_DIR=/path/to/bams,PRESET=asm20 \
  sbatch/map_assemblies_to_reference_array.sbatch
```

## Step 2. Prepare the SNP-Chip Panel

The variant-calling sbatch script automatically prepares reusable panel files:

```text
work/panel/panel.biallelic.snps.vcf.gz
work/panel/panel.positions.tsv
work/panel/panel.alleles.tsv.gz
```

Conceptually, this step does:

```bash
bcftools view -m2 -M2 -v snps -Oz \
  -o work/panel/panel.biallelic.snps.vcf.gz \
  data/snp_chip_panel.vcf.gz

bcftools index -f -t work/panel/panel.biallelic.snps.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' \
  work/panel/panel.biallelic.snps.vcf.gz \
  | bgzip -c > work/panel/panel.alleles.tsv.gz

tabix -f -s1 -b2 -e2 work/panel/panel.alleles.tsv.gz
```

The allele file is important for `bcftools call -C alleles`. It tells `bcftools` which alleles are expected at each marker, so the query assembly is genotyped against the SNP-chip allele definitions rather than inventing a different representation.

## Step 3. Genotype Each Assembly at Panel Sites

Use:

```bash
sbatch sbatch/call_panel_variants_and_gtcheck.sbatch
```

Default paths:

```text
REFERENCE_FASTA=reference/genome.fa
PANEL_VCF=data/snp_chip_panel.vcf.gz
BAM_DIR=bams
RESULTS_DIR=results
PANEL_DIR=work/panel
```

Example with explicit paths:

```bash
sbatch \
  --export=ALL,REFERENCE_FASTA=/path/to/ref.fa,PANEL_VCF=/path/to/chip.vcf.gz,BAM_DIR=/path/to/bams,RESULTS_DIR=/path/to/results \
  sbatch/call_panel_variants_and_gtcheck.sbatch
```

### The Central `bcftools` Pipeline

The important command chain is:

```bash
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
```

### Reasoning Through the Parameters

`-f "$REFERENCE_FASTA"` tells `mpileup` which reference bases to use. This must match the coordinate system used for the BAM and the panel.

`-T "$panel_snps"` restricts pileup generation to panel markers. This is the heart of the workflow. We are not discovering all variants; we are asking for genotypes at a predefined set of positions.

`--no-BAQ` disables BAQ adjustment. BAQ is often useful for noisy short-read alignments near indels, but assembly-to-reference BAMs do not behave like ordinary read alignments. For panel-site identity checking, turning BAQ off keeps interpretation simpler.

`-Q 0` keeps bases regardless of base quality. Assembly alignments often have uninformative or artificial base qualities. A normal resequencing threshold such as `-Q 20` can throw away useful assembly-derived bases for no biological reason.

`-q "$MIN_MQ"` filters low mapping-quality alignments. The default `MIN_MQ=5` is mild. It removes the most ambiguous contig placements while retaining many usable sites. If too many sites become missing, try `MIN_MQ=1`. If repetitive regions are causing suspicious calls, try a higher value.

`-C alleles -T "$panel_alleles"` tells `bcftools call` to genotype known alleles from the panel. This helps keep the query VCF compatible with the SNP-chip database.

`-i 1` asks `bcftools call` to emit target sites. Without this idea, you risk drifting back toward a sparse variant-only VCF, which is exactly what this workflow is designed to avoid.

`-V indels` keeps the callset SNP-focused.

`-a GQ` adds genotype quality, which is later used to set weak genotypes to missing.

### Optional Array-Job Panel Calling and `gtcheck`

For larger projects, use the array version after the BAM files have been created.

Create a BAM manifest:

```bash
find bams -maxdepth 1 -name '*.bam' | sort > bams.fofn
```

Submit one task per BAM:

```bash
N=$(wc -l < bams.fofn)
sbatch --array=1-"$N" sbatch/call_panel_variants_and_gtcheck_array.sbatch
```

The array script uses the same default variables as the loop script, plus:

```text
BAM_MANIFEST=bams.fofn
```

The SNP-panel preparation files in `work/panel/` are shared by all array tasks. The array script uses a small lock directory so only one task prepares those files while the others wait. If you prefer a simpler production pattern, run the loop script once on a tiny test BAM to prepare the panel files, then submit the full array.

## Step 4. Set Weak Calls to Missing

The script filters the query VCF like this:

```bash
bcftools filter \
  -S . \
  -e "FMT/DP<1 || FMT/GQ<${MIN_GQ}" \
  -Oz \
  -o "$filtered_vcf" \
  "$named_vcf"
```

This is an important design choice.

For identity checking, it is usually better to keep the marker record but set a weak genotype to missing (`./.`) than to delete the marker entirely. Keeping the record preserves the structure of the panel-site VCF and makes missingness easier to audit.

The default threshold is:

```text
MIN_GQ=20
```

That is a reasonable starting point, not a universal truth. Assemblies, panels, and references differ. If your call rate is extremely low, inspect the BAMs and consider whether `MIN_MQ`, `MIN_GQ`, reference compatibility, or chromosome naming is the issue.

## Step 5. Compare to the Panel with `bcftools gtcheck`

The script runs:

```bash
bcftools gtcheck \
  -u GT,GT \
  -E 0 \
  -O t \
  -o "$gtcheck_tsv" \
  -g "$panel_snps" \
  "$filtered_vcf"
```

`gtcheck` compares the assembly-derived query genotypes to each sample in the panel VCF.

The options matter:

- `-g "$panel_snps"` provides the known panel genotypes.
- `-u GT,GT` forces genotype-vs-genotype comparison.
- `-E 0` makes discordance easier to interpret.
- `-O t` writes tab-delimited output.

With `-u GT,GT -E 0`, the discordance value is interpretable as the count of mismatching genotypes among compared sites. That is much easier to explain than a likelihood-based score.

## Step 6. Summarize the Best Matches

Run the Python summary script:

```bash
python scripts/summarize_gtcheck_top_hits.py \
  -i 'results/*.gtcheck.tsv' \
  -n 10 \
  --min-sites 1000 \
  -o results/gtcheck_top10.tsv
```

For the included soybean example outputs:

```bash
python scripts/summarize_gtcheck_top_hits.py \
  -i 'examples/results/*.gtcheck.tsv' \
  -n 5 \
  --min-sites 1000 \
  -o examples/results/example_gtcheck_top5.tsv
```

The script is species-agnostic. It does not assume soybean, SoySNP50K, a particular filename pattern, or a particular panel size. It reads standard `bcftools gtcheck` `DCv2` rows.

### Useful Summary Columns

`query_sample`: the sample name from the query VCF.

`panel_sample`: the sample in the SNP-chip database being compared to the query.

`rank`: best hit is rank 1.

`discordance`: number of mismatching genotypes in GT-vs-GT mode with `-E 0`.

`sites_compared`: number of marker sites that could actually be compared for that query-panel pair.

`matching_genotypes`: number of compared marker sites that matched.

`discordance_per_site`: `discordance / sites_compared`; lower is better.

`match_fraction`: `matching_genotypes / sites_compared`; higher is better.

`confidence_score`: `match_fraction * log10(sites_compared)`. This is a simple helper score that rewards high concordance while slightly penalizing tiny site counts. Do not treat it as a formal statistic.

`match_fraction_gap_rank1_rank2`: in the sample summary file, the difference between the best and second-best hit. A large gap is reassuring; a tiny gap suggests ambiguity or close relatedness.

## Interpreting Results

First look at `sites_compared`. A beautiful match fraction based on 12 sites is not convincing. A slightly lower match fraction based on 20,000 sites may be much stronger evidence.

Then look at `match_fraction`. For many clean within-species SNP-chip identity checks, a strong match may be above 0.98 or 0.99, but thresholds depend on panel density, heterozygosity, reference quality, imputation, and the biological material.

Finally, look at the gap between rank 1 and rank 2.

### Strong Match Pattern

```text
rank 1   match_fraction 0.995   sites_compared 18000
rank 2   match_fraction 0.962   sites_compared 17950
```

This is usually reassuring. The query assembly has a high-concordance best hit and a clear separation from the next candidate.

### Ambiguous Match Pattern

```text
rank 1   match_fraction 0.991   sites_compared 16000
rank 2   match_fraction 0.989   sites_compared 16010
rank 3   match_fraction 0.988   sites_compared 15980
```

This may mean the panel contains closely related lines, near-duplicates, derivatives, or accessions with shared ancestry. It can still be biologically meaningful, but it is not a clean unique identity call.

### Weak Match Pattern

```text
rank 1   match_fraction 0.780   sites_compared 300
rank 2   match_fraction 0.775   sites_compared 290
```

This is not enough evidence for a confident identity assignment. Investigate missingness, chromosome naming, reference mismatch, panel compatibility, or the possibility that the assembly is not represented well in the panel.

## Why Site Counts Differ

Students often expect `sites_compared` to equal the number of SNP-chip markers. It usually does not.

Sites can be lost because:

- the assembly does not align over that marker,
- the aligned contig has low mapping quality,
- the genotype was set to missing after filtering,
- the panel and reference use incompatible chromosome names,
- the REF/ALT representation does not match cleanly,
- the genotype is not diploid or otherwise not usable by `gtcheck`,
- the query assembly is structurally different from the reference in that region.

This is why the workflow writes both filtered VCFs and small QC files. If a sample has unexpectedly low `sites_compared`, check how many panel sites were called in the query VCF before interpreting the identity result.

## Adapting to Different HPC Systems

The sbatch scripts are examples, not universal cluster law. Edit:

```bash
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
```

and the module lines:

```bash
module load minimap2
module load samtools
module load bcftools
module load htslib
```

On SciNet, Sapelo2, HyperGator, or another university cluster, the exact partition names and module names may differ. The workflow logic stays the same.

For many assemblies, consider converting the loop-style scripts into SLURM array jobs. Arrays are easier to restart because each sample has its own task and log. The current scripts are deliberately loop-style because they are easier for a new student to read from top to bottom.

See [docs/hpc_notes.md](docs/hpc_notes.md) for additional notes on Sapelo2-style modules, SciNet-style environments, HyperGator-style account settings, scratch-space habits, and converting the loop examples into array jobs.

## Testing the Parser

Before using real results, you can run a tiny synthetic `gtcheck` parser test:

```bash
bash tests/run_tiny_test.sh
```

This writes a temporary output file in `tests/tmp/`, compares it to `tests/expected/tiny_top_hits.tsv`, and prints a success message if the parser ranking is working as expected.

## Common Problems and Fixes

### `bcftools` says chromosome names do not match

Your reference FASTA and panel VCF use different contig names. Rename the VCF contigs or use the correct reference version.

### `gtcheck` compares very few sites

Check:

- Is the panel VCF indexed?
- Does the BAM use the same reference coordinates as the panel?
- Are BAM indexes present?
- Are `MIN_MQ` or `MIN_GQ` too strict?
- Did `bcftools call` emit panel sites, or did it produce a sparse VCF?
- Are many genotypes `./.` in the filtered query VCF?

### The top two hits are nearly identical

That may be real. Many crop panels include related breeding lines, derivatives, duplicates, or samples with subtle naming differences. Treat the result as a cluster-level identity clue unless rank 1 separates clearly from rank 2.

### The expected sample is absent from the panel

The workflow can still find the closest available genotype, but it cannot prove identity to a sample that is not represented in the database.

## Minimal End-to-End Example

Prepare directories:

```bash
mkdir -p reference assemblies data bams work results logs
```

Place files:

```text
reference/genome.fa
assemblies/sample1.fa
assemblies/sample2.fa
data/snp_chip_panel.vcf.gz
data/snp_chip_panel.vcf.gz.tbi
```

Index the reference:

```bash
samtools faidx reference/genome.fa
```

Map assemblies:

```bash
sbatch sbatch/map_assemblies_to_reference.sbatch
```

Call panel-site genotypes and run `gtcheck`:

```bash
sbatch sbatch/call_panel_variants_and_gtcheck.sbatch
```

Summarize top hits:

```bash
python scripts/summarize_gtcheck_top_hits.py \
  -i 'results/*.gtcheck.tsv' \
  -n 10 \
  --min-sites 1000 \
  -o results/gtcheck_top10.tsv
```

## Teaching Notes

The central lesson is that file format choices encode biological assumptions.

A variant-only VCF says, “show me where this sample differs from the reference.” That is excellent for variant discovery, but incomplete for identity checking.

A targeted panel-site VCF says, “show me the genotype at each marker used by the database.” That is the right shape of evidence for comparing an assembly to a SNP-chip panel.

Once students understand that distinction, the rest of the workflow becomes logical:

1. Align the assembly so each contig can be interpreted in reference coordinates.
2. Query only the SNP-chip marker positions because those are the coordinates shared with the database.
3. Preserve callable reference genotypes because matches to the reference are still genotype information.
4. Set weak genotypes to missing rather than pretending they are confident.
5. Compare genotype profiles to the panel.
6. Interpret both concordance and the number of sites supporting that concordance.

That reasoning applies whether the crop is soybean, cotton, maize, wheat, or another SNP-chip project.
