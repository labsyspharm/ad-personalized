#!/bin/bash
#SBATCH -J samtools_bam_to_fastq
#SBATCH -c 8                               # Request one core
#SBATCH -t 0-0:45                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/samtools_bam_to_fastq_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/samtools_bam_to_fastq_%j.err                 # File to which STDERR will be written, including job ID (%j)

# USAGE:
# comm -1 -3 <(ls -1 fastq | grep .done$ | sed s/.done//g | sort | uniq) <(ls -1 raw | grep .bam$ | grep -v sorted | sed s/.bam//g | sort | uniq) |
#  parallel -N 8 sbatch bam_to_fastq.sh {}


module load gcc/14.2.0 samtools/1.21

set -eux

# for i in SRR8571937  SRR8571938 SRR8571939 SRR8571940 SRR8571941 SRR8571942 SRR8571944 SRR8571945 SRR8571947 SRR8571948 SRR8571949 SRR8571950 SRR8571951 SRR8571952
for i in $@
do
echo "$SLURM_JOB_ID" > "fastq/${i}.lock"
if [ ! -f "raw/${i}_namesorted.bam" ]; then
  samtools collate -@ 8 -o "raw/${i}_namesorted.bam" "raw/${i}.bam"
fi
if [ ! -f "fastq/${i}.done" ]; then
  samtools fastq -1 >(pigz -p 3 -c > "fastq/${i}_1.fastq.gz") \
    -2 >(pigz -p 3 -c > "fastq/${i}_2.fastq.gz") \
    -0 >(pigz -p 3 -c > "fastq/${i}_unpaired.fastq.gz") \
    -@ 4 "raw/${i}_namesorted.bam"
  rm "raw/${i}_namesorted.bam"
fi
rm "fastq/${i}.lock"
echo "$SLURM_JOB_ID" > "fastq/${i}.done"
done

echo "done"
