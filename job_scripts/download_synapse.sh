#!/bin/bash
#SBATCH -J synapse_download
#SBATCH -c 4                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)
#SBATCH -o synapse_download_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e synapse_download_%j.err                 # File to which STDERR will be written, including job ID (%j)
set -eux

# Give a CSV file ($1) with synapse IDs in given column number ($2) and download them all
# EXAMPLE sbatch download_synapse.sh selected_rosmap_samples.csv 2

source ~/miniconda3/etc/profile.d/conda.sh
conda activate synapse

cut -d, -f "$2" "$1" |
  tail -n +2 |
  parallel -j 4 --bar --halt now,fail=1 \
    synapse get --downloadLocation raw {}
