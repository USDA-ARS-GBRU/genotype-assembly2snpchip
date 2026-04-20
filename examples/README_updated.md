# SoySNP50K Assembly Identity Validation Pipeline

A reproducible workflow for validating the identity of soybean whole-genome assemblies against the SoySNP50K genotype panel.

This pipeline answers a simple but important question:

> Does this assembly match the accession label we think it does?

It is designed for cases where assemblies have already been aligned to the Wm82a6 reference genome, and the goal is to compare each sample against the large SoySNP50K reference panel with `bcftools gtcheck`.

---

## Why this workflow exists

A common first instinct is to:

1. align each assembly to Wm82a6,
2. call a whole-genome VCF,
3. compare that VCF against the SoySNP50K panel.

That can work for more divergent accessions, but it breaks down for near-reference samples.

Why?

- Whole-genome assembly-derived VCFs are usually **variant-only**.
- Near-reference genomes can therefore produce **very sparse query VCFs**.
- `bcftools gtcheck` then ends up comparing **very few sites**, or even **zero sites**.

That is exactly why a self-mapping control like Wm82a6 mapped back to Wm82a6 can produce a `gtcheck` file with zero informative comparisons.

The fix is to turn the problem around:

> Instead of asking which variants happen to appear in each assembly, genotype every assembly at the exact SoySNP50K marker positions.

That gives every sample a standardized marker profile and makes `gtcheck` much more interpretable.

---

## Workflow logic

The workflow has four conceptual stages:

1. **Align each assembly to Wm82a6 with minimap2**
2. **Restrict genotyping to SoySNP50K marker sites**
3. **Run `bcftools gtcheck` against the SoySNP50K panel**
4. **Summarize the top hits per query sample**

In shorthand:

```text
Assembly FASTA
   -> minimap2 assembly-to-reference alignment
   -> BAM
   -> targeted bcftools mpileup/call at SoySNP50K sites
   -> SoySNP50K-only query VCF
   -> bcftools gtcheck against panel
   -> summarized top hits
```

---

## Project layout

This README assumes a directory structure like this:

```text
00_log/                 # SLURM logs
01_sbatch/              # sbatch scripts
02_scripts/             # helper scripts
03_genomes/             # reference genome and optionally assembly FASTAs
04b_bams-asm20/         # assembly-to-genome BAMs
05_vcfs_soy50k/         # targeted SoySNP50K VCFs and gtcheck output
50/                     # SoySNP50K input VCF
```

Main inputs used below:

- reference genome: `03_genomes/Wm82a6.fa`
- SoySNP50K panel: `50/50k.renamed.vcf.gz`
- BAM directory: `04b_bams-asm20/`

---

## Step 1. Align assemblies to Wm82a6 with minimap2

If your BAMs already exist, you can skip this step.

The alignment mode used here is:

```bash
minimap2 -ax asm20
```

This tells minimap2 to align assembled contigs/scaffolds against a reference genome with settings tuned for assemblies that may be moderately diverged.

### Why `asm20`?

`asm20` is usually a good default for soybean assembly-to-reference comparisons when the samples are not expected to be nearly identical clones but are still within the same species.

Roughly speaking:

- `asm5` is stricter, for more similar assemblies
- `asm10` is intermediate
- `asm20` is more tolerant of divergence

For this use case, `asm20` is a sensible default.

### Existing loop-style alignment script

The current `minimap-A2G.sbatch` follows a one-sample-at-a-time pattern driven by environment variables such as `query`, `sample`, `oDir`, `refG`, and `asm`. The key alignment step is:

```bash
minimap2 -ax $asm -t 16 $refG $query > $oDir/$sample.sam
samtools view -bS $oDir/$sample.sam | samtools sort -o $oDir/$sample.sorted.bam
samtools index $oDir/$sample.sorted.bam
```

That is fine, but for large runs you will usually want either:

- a **loop-style sbatch** that processes every assembly within one job, or
- an **array-style sbatch** that processes one assembly per task.

### Recommended loop-style pattern

```bash
#!/bin/bash
#SBATCH -J minimap_A2G
#SBATCH -p batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=06:00:00
#SBATCH -o 00_log/%x_%j.out
#SBATCH -e 00_log/%x_%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

ml minimap2/2.28-GCCcore-13.2.0
ml SAMtools/1.21-GCC-13.3.0

refG="03_genomes/Wm82a6.fa"
indir="03_genomes/assemblies"
outDir="04b_bams-asm20"
asm="asm20"

mkdir -p "$outDir"
shopt -s nullglob

for query in "$indir"/*.fa "$indir"/*.fasta "$indir"/*.fna; do
    [[ -e "$query" ]] || continue
    sample=$(basename "$query")
    sample=${sample%.fa}
    sample=${sample%.fasta}
    sample=${sample%.fna}

    bam="$outDir/${sample}-asm20.sorted.bam"

    minimap2 -ax "$asm" -t "$SLURM_CPUS_PER_TASK" "$refG" "$query" \
        | samtools sort -@ 8 -o "$bam"

    samtools index "$bam"
done
```

### Recommended array-style pattern

Arrays are often cleaner on SLURM because each sample gets its own task, log, and failure state.

First, make a manifest:

```bash
find 03_genomes/assemblies -maxdepth 1 \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) | sort > 03_genomes/assemblies.fofn
```

Then use an array job:

```bash
#!/bin/bash
#SBATCH -J minimap_A2G
#SBATCH -p batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=06:00:00
#SBATCH --array=1-20
#SBATCH -o 00_log/%x_%A_%a.out
#SBATCH -e 00_log/%x_%A_%a.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

ml minimap2/2.28-GCCcore-13.2.0
ml SAMtools/1.21-GCC-13.3.0

refG="03_genomes/Wm82a6.fa"
outDir="04b_bams-asm20"
manifest="03_genomes/assemblies.fofn"
asm="asm20"

mkdir -p "$outDir"

query=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$manifest")
[[ -n "$query" ]] || { echo "No query found for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

sample=$(basename "$query")
sample=${sample%.fa}
sample=${sample%.fasta}
sample=${sample%.fna}

bam="$outDir/${sample}-asm20.sorted.bam"

minimap2 -ax "$asm" -t "$SLURM_CPUS_PER_TASK" "$refG" "$query" \
    | samtools sort -@ 8 -o "$bam"

samtools index "$bam"
```

Submit with the correct array size for your manifest.

---

## Step 2. Build a targeted SoySNP50K marker panel

The panel VCF may contain records you do not want for this identity-checking workflow. The cleanest target set is:

- SNPs only
- biallelic sites only

Build that panel once:

```bash
mkdir -p 05_vcfs_soy50k/panel

bcftools view \
    -m2 -M2 \
    -v snps \
    -Oz \
    -o 05_vcfs_soy50k/panel/soy50k.biallelic.snps.vcf.gz \
    50/50k.renamed.vcf.gz

bcftools index -f -t 05_vcfs_soy50k/panel/soy50k.biallelic.snps.vcf.gz

bcftools query \
    -f '%CHROM\t%POS\n' \
    05_vcfs_soy50k/panel/soy50k.biallelic.snps.vcf.gz \
    > 05_vcfs_soy50k/panel/soy50k.positions.tsv
```

---

## Step 3. Genotype each BAM only at SoySNP50K sites

This is the key step.

Instead of calling all variants in the BAM, we restrict genotyping to the SoySNP50K marker list.

### Why this is better than a whole-genome VCF

A whole-genome assembly VCF tells you where the sample differs from the reference.

A targeted SoySNP50K VCF tells you:

> what genotype the sample has at each predefined panel site

That means reference matches (`0/0`) are preserved and become informative, which is exactly what identity checking needs.

### Main `mpileup` settings

The recommended targeted call uses:

```bash
bcftools mpileup \
    -f reference.fa \
    -T soy50k.biallelic.snps.vcf.gz \
    --no-BAQ \
    -Q 0 \
    -q 5 \
    -a FORMAT/AD,FORMAT/DP,FORMAT/PL
```

Here is what the key filters mean.

#### `-q, --min-MQ`

```text
-q INT    Skip alignments with mapQ smaller than INT
```

This filters **entire alignments** based on mapping quality.

In resequencing data, that means low-confidence read placements. In this workflow, where BAMs come from assembly contigs aligned back to the reference, it mostly means:

> ignore ambiguous contig placements in repeats, paralogs, or collapsed regions.

That is useful for identity checking. `-q 5` is a sensible default. If too many sites become missing, try `-q 1`.

#### `-Q, --min-BQ`

```text
-Q INT    Skip bases with baseQ/BAQ smaller than INT
```

This filters **individual bases** by base quality.

For normal short-read resequencing, base quality matters a lot. For assembly-alignment BAMs, base qualities are often absent, constant, or not especially meaningful.

So for this workflow:

> `-Q 0` is the right default.

It keeps all aligned bases and avoids discarding usable positions for arbitrary reasons.

#### `--no-BAQ`

BAQ is often useful for noisy read alignments near small indels, but it is less helpful here and can complicate interpretation. For assembly-to-reference genotyping at predefined SNP sites, turning BAQ off is a reasonable choice.

### Main `call` settings

The recommended call stage uses:

```bash
bcftools call \
    -C alleles \
    -T soy50k.biallelic.snps.vcf.gz \
    -i 1 \
    -V indels \
    -f GQ
```

Why:

- `-C alleles` tells `bcftools call` to genotype the known alleles in the target panel.
- `-T panel.vcf.gz` ensures those exact SoySNP50K sites are the calling targets.
- `-i 1` helps preserve target-site output rather than silently collapsing the query back into a sparse variant-only file.
- `-V indels` keeps this workflow SNP-only.
- `-f GQ` adds genotype quality.

### Filtering strategy

A good beginner mistake to avoid is removing low-confidence target sites entirely.

For identity checking, it is usually better to:

- keep the site in the VCF,
- but set the genotype to missing if it is too weak.

Example:

```bash
bcftools filter \
    -S . \
    -e 'FMT/DP<1 || FMT/GQ<20' \
    -Oz \
    -o sample.soy50k.filtered.vcf.gz \
    sample.soy50k.named.vcf.gz
```

For assembly-alignment BAMs, depth is often low and should not be judged like resequencing coverage. A filter like `DP<1` simply requires that something aligned there at all.

---

## Step 4. Compare against the SoySNP50K panel with `bcftools gtcheck`

Recommended settings:

```bash
bcftools gtcheck \
    -u GT,GT \
    -E 0 \
    -g 05_vcfs_soy50k/panel/soy50k.biallelic.snps.vcf.gz \
    sample.soy50k.filtered.vcf.gz
```

### Why `-u GT,GT -E 0` matters

Without care, `gtcheck` can compare query `PL` values against panel `GT` values. That can still rank hits, but the discordance value becomes a score rather than a literal mismatch count.

Using:

- `-u GT,GT`
- `-E 0`

makes the comparison much more interpretable:

> discordance becomes the number of mismatching genotypes among the compared sites.

That is the cleanest mode for PI-facing interpretation.

---

## Step 5. Summarize the top hits

The repository already includes a companion script:

```text
summarize_gtcheck_top_hits_hardened.py
```

Its job is to parse one or more `bcftools gtcheck` output files and report the top N best hits per query sample.

This script assumes bcftools gtcheck was run in GT,GT mode (-u GT,GT -E 0). If not, discordance values are not directly interpretable as mismatch counts.

### What the script does well

The script is sound for its main purpose:

- it parses `DCv2` rows from `gtcheck`
- it calculates `discordance_per_site`
- it calculates `match_fraction`
- it ranks hits per query sample
- it writes a combined TSV output

The ranking order is:

1. `discordance_per_site` ascending
2. `discordance` ascending
3. `avg_neg_log_hwe` descending
4. `matching_genotypes` descending

That is a reasonable ranking strategy.

### Two edge cases to be aware of

The script is basically good, but there are two edge cases worth knowing about:

1. **Zero-site rows**
   - if `sites_compared == 0`, then `discordance_per_site` and `match_fraction` become undefined or infinite.
   - that can happen with old variant-only workflows, especially for self-reference controls.

2. **Filtering everything out with `--min-sites`**
   - if all rows are removed after filtering, the script can end up with no groups to concatenate.

The hardened script automatically removes rows with sites_compared == 0 by default, so these edge cases are handled safely. You can override this behavior with --keep-zero-sites if needed.

### Recommended usage

```bash
python summarize_gtcheck_top_hits_hardened.py \
  -i 05_vcfs_soy50k/*.gtcheck.tsv \
  --min-sites 1000 \
  -n 10 \
  -o 05_vcfs_soy50k/gtcheck_top10.tsv
```

If your targeted callsets are dense and clean, `--min-sites 1000` is a good default safety valve. If needed, you can raise it.

### Main output columns

The summary file includes:

- `query_sample`
- `query_sample_basename`
- `genotyped_sample`
- `rank`
- `discordance`
- `discordance_per_site`
- `avg_neg_log_hwe`
- `sites_compared`
- `matching_genotypes`
- `match_fraction`
- `source_file`

### Most useful columns for interpretation

For practical identity checking, focus first on:

- `match_fraction`
- `sites_compared`
- `rank`
- the gap between rank 1 and rank 2

Note: match_fraction is computed by the companion script and is robust to missing data, but is only meaningful when sites_compared is sufficiently large.

---

## A complete SLURM workflow

### Option A. Run the targeted genotyping script as one loop-style job

This is the simplest approach. One sbatch script walks over all BAMs in `04b_bams-asm20/`.

Submit:

```bash
sbatch 01_sbatch/soy50k_gtcheck_from_asm20.sbatch
```

This style is easy to manage and works well for a modest number of BAMs.

### Option B. Run targeted genotyping as an array job

Array jobs are cleaner if you want:

- one log per sample
- easier restart of failed samples
- better scheduler visibility

First build a BAM manifest:

```bash
find 04b_bams-asm20 -maxdepth 1 -name '*.sorted.bam' | sort > 04b_bams-asm20/bams.fofn
```

Then adapt the script so each task reads one BAM from the manifest:

```bash
#!/bin/bash
#SBATCH -J soy50k_gtcheck
#SBATCH -p batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=1-13
#SBATCH -o 00_log/%x_%A_%a.out
#SBATCH -e 00_log/%x_%A_%a.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load BCFtools/1.21-GCC-13.3.0

refG="03_genomes/Wm82a6.fa"
panelVCF="50/50k.renamed.vcf.gz"
outDir="05_vcfs_soy50k"
panelDir="05_vcfs_soy50k/panel"
manifest="04b_bams-asm20/bams.fofn"

mkdir -p 00_log "$outDir" "$panelDir"

panelSNP="${panelDir}/soy50k.biallelic.snps.vcf.gz"

if [[ ! -s "$panelSNP" || ! -s "${panelSNP}.tbi" ]]; then
    bcftools view -m2 -M2 -v snps -Oz -o "$panelSNP" "$panelVCF"
    bcftools index -f -t "$panelSNP"
fi

bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$manifest")
[[ -n "$bam" ]] || { echo "No BAM found for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

base=$(basename "$bam")
sample=${base%.sorted.bam}
sample=${sample%-asm20}

rawVCF="${outDir}/${sample}.soy50k.raw.vcf.gz"
namedVCF="${outDir}/${sample}.soy50k.named.vcf.gz"
filtVCF="${outDir}/${sample}.soy50k.filtered.vcf.gz"
gtcheckTSV="${outDir}/${sample}.soy50k.filtered.gtcheck.tsv"
sampleNameFile="${outDir}/${sample}.sample_name.txt"

bcftools mpileup \
    -Ou \
    --threads "$SLURM_CPUS_PER_TASK" \
    -f "$refG" \
    -T "$panelSNP" \
    --no-BAQ \
    -Q 0 \
    -q 5 \
    --max-depth 100 \
    -a FORMAT/AD,FORMAT/DP,FORMAT/PL \
    "$bam" \
| bcftools call \
    -Ou \
    -C alleles \
    -T "$panelSNP" \
    -i 1 \
    -V indels \
    -f GQ \
| bcftools norm \
    -Ou \
    -f "$refG" \
| bcftools annotate \
    -x INFO,^FORMAT/GT,FORMAT/DP,FORMAT/AD,FORMAT/PL,FORMAT/GQ \
    -Oz \
    -o "$rawVCF"

bcftools index -f -t "$rawVCF"

printf '%s\n' "$sample" > "$sampleNameFile"
bcftools reheader -s "$sampleNameFile" -o "$namedVCF" "$rawVCF"
bcftools index -f -t "$namedVCF"

bcftools filter \
    -S . \
    -e 'FMT/DP<1 || FMT/GQ<20' \
    -Oz \
    -o "$filtVCF" \
    "$namedVCF"

bcftools index -f -t "$filtVCF"

bcftools gtcheck \
    -u GT,GT \
    -E 0 \
    -O t \
    -o "$gtcheckTSV" \
    -g "$panelSNP" \
    "$filtVCF"
```

Then submit with the array size matching the number of BAMs.

---

## Interpreting `gtcheck` output

A typical `gtcheck` result row contains:

```text
query_sample  panel_sample  discordance  HWE_score  sites_compared  matching_genotypes
```

For this workflow, the three most important interpretation columns are:

### `sites_compared`

How many SoySNP50K markers were actually usable for that pairwise comparison.

- high = strong evidence
- low = weak evidence

As a rough guide:

- thousands of sites: good
- a few hundred: borderline
- single digits: not trustworthy

### `matching_genotypes`

How many of those compared sites had matching genotypes.

### `match_fraction`

This is not a native `bcftools gtcheck` column. It is added by the summary script and calculated as:

```text
matching_genotypes / sites_compared
```

This is the most intuitive concordance measure for day-to-day interpretation.

---

## How to decide whether a sample looks right

### Strong, reassuring result

A convincing identity assignment usually has:

- high `match_fraction`
- large `sites_compared`
- a clear gap between rank 1 and rank 2

Example pattern:

```text
rank 1   match_fraction 0.995   sites_compared 14500
rank 2   match_fraction 0.961   sites_compared 14500
```

That is a strong top hit.

### Ambiguous but still informative result

Example pattern:

```text
rank 1   0.992
rank 2   0.989
rank 3   0.988
```

That suggests a cluster of closely related lines rather than a clean one-to-one identity.

### Weak or suspicious result

Example pattern:

- low `match_fraction`
- few compared sites
- no clear separation between top hits

Possible explanations:

- poor BAM quality
- too many missing sites
- contamination
- mislabeled sample
- a sample that is simply not well represented by the panel

---

## Controls and sanity checks

### Wm82a6 mapped to Wm82a6

In the original variant-only workflow, this often produces zero compared sites and tells you almost nothing.

In the targeted SoySNP50K workflow, it should become informative because reference-state calls are preserved.

### Wm82v4 mapped to Wm82a6

This should also become much more interpretable in the targeted workflow.

If both Wm82 controls behave sensibly, that is a strong sign that the pipeline is doing what it should.

---

## Practical advice for beginners

1. **Do not over-filter.** In this workflow, missing is usually better than deleting a site.
2. **Do not interpret old variant-only Wm82 self-maps as biological evidence.** They are mostly workflow artifacts.
3. **Always look at `sites_compared` before trusting a hit.**
4. **Treat the gap between rank 1 and rank 2 as part of the evidence.**
5. **Use known controls first.** If Wm82a6 and Wm82v4 look sensible, then the rest of the run is easier to trust.
6. **Prefer array jobs when you want robust reruns and sample-level logs.** Prefer loops when you want simplicity.

---

## Suggested command sequence

### 1. Align assemblies

```bash
sbatch 01_sbatch/minimap-A2G.sbatch
```

Or use an array-based alignment script if preferred.

### 2. Run targeted SoySNP50K genotyping and gtcheck

```bash
sbatch 01_sbatch/soy50k_gtcheck_from_asm20.sbatch
```

### 3. Summarize top hits

```bash
python summarize_gtcheck_top_hits_hardened.py \
  -i 05_vcfs_soy50k/*.gtcheck.tsv \
  --min-sites 1000 \
  -n 10 \
  -o 05_vcfs_soy50k/gtcheck_top10.tsv
```

---

## Recommended future improvements

This repository is already quite usable, but the following additions would make it even stronger:

- add per-sample QC summaries (call rate, heterozygosity) into the final table
- generate a PI-facing PDF or HTML report with flagged suspicious samples
- optionally add PCA validation on the targeted marker matrix

---

## Bottom line

This pipeline turns assembly alignments into **panel-compatible genotype profiles** and compares them against the SoySNP50K reference population in a way that is much more interpretable than a generic whole-genome variant workflow.

That is exactly what you want when the real biological question is not “what variants are here?” but:

> “Is this assembly actually the accession we think it is?”
