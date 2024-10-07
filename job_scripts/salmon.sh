#!/bin/bash
#SBATCH -J salmon
#SBATCH -c 5                               # Request one core
#SBATCH -t 0-0:40                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

# USAGE EXAMPLE
# sbatch salmon.sh /n/scratch3/users/c/ch305/rna-editing/rna-seq/fastq/Homo_sapiens.GRCh38.gentrome_including_variants_index 962_TCX_2

# use python provided by miniconda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate salmon

# sets options how shell script is executed:
# e = exit immediately on error
# u = treat using undeclared variables as error
# x = print commands to log before running them
set -eux

PREFIX_PATH="$2"
PREFIX_FILE="${PREFIX_PATH##*/}"

# make a temporary lock file so I can see that this sample is
# currently being worked on
touch "quants/${2}.lock"

# $1 contains path to salmon index
# $2 contains sample name
# trimming outputs two fastq files. one contains paired reads,
# the other reads whose pair got lost during trimming.
# I'm using `cat` to concatenate them back together for alignment with salmon

salmon quant -i "$1" -l A \
  --seqBias --gcBias --posBias --recoverOrphans -p 5 \
  -r <(cat "trimmed/${2}_1_unpaired.fastq.gz" "trimmed/${2}_2_unpaired.fastq.gz")
  -1 "trimmed/${2}_1.fastq.gz" \
  -2 "trimmed/${2}_2.fastq.gz" \
  -o "quants/${2}"

# Marking sample as done on the file system
rm "quants/${2}.lock"
touch "quants/${2}.done"
