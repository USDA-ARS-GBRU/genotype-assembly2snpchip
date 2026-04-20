#!/bin/bash
#SBATCH -J Minimap-cntgs
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=150G
#SBATCH -t 6:00:00
#SBATCH -o 00_log/minimap.out
#SBATCH -e 00_log/minimap.err

cd $SLURM_SUBMIT_DIR

echo Input Seqs: $query
echo Sample: $sample
echo Output Directory: $oDir
echo Ref Genome: $refG
echo Assembly Level: $asm
echo Node: $SLURMD_NODENAME

mkdir -p ${oDir}

# load minimap2 and samtools modules
ml minimap2/2.28-GCCcore-13.2.0
ml SAMtools/1.21-GCC-13.3.0

# map genome assembly to reference genome
# choose one of asm5, asm10, or asm20 for the asm param
minimap2 -ax $asm -t 16 $refG $query > $oDir/$sample.sam


echo ' '
echo minimap2 version used: `minimap2 --version`
echo ' '

# sort and compress and index
samtools view -bS $oDir/$sample.sam | samtools sort -o $oDir/$sample.sorted.bam
samtools index $oDir/$sample.sorted.bam

echo ' '
echo samtools version used: `samtools --version`
echo ' '

# attach the job ID to the log
# check job resource usage with: seff <Jobid>
echo ' '
echo Job ID: $SLURM_JOBID
echo ' '