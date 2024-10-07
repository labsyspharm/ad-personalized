
# module load sratoolkit/3.0.2
#
# parallel -j 6 --bar \
#   fastq-dump -O fastq --split-3 --gzip \
#     {} ::: SRR8571937 SRR8571938 SRR8571939 SRR8571940 SRR8571941 SRR8571942 SRR8571943 SRR8571944 SRR8571945 SRR8571946 SRR8571947 SRR8571948 SRR8571949 SRR8571950 SRR8571951 SRR8571952

# Download from ENA instead of from SRA, much faster

module load gcc/9.2.0 R/4.2.1

# First download TSV file with fastq URLs from ENA and turn it into
# a simple list of URLs using R

R -e 'library(tidyverse)
tbl <- read_tsv("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA522295&result=read_run&fields=fastq_ftp,submitted_ftp&format=tsv&download=true&limit=0") %>%
  separate(fastq_ftp, c("fastq_1", "fastq_2"), sep = ";") %>%
  pivot_longer(-c(run_accession, submitted_ftp), names_to = "read", values_to = "url")
write_lines(
  tbl$url, "gse126542_fastq_urls.txt"
)
'

cd fastq

# Download the fastqs. Very fast!

parallel -j 6 --bar \
  wget "ftp://"{} < ../gse126542_fastq_urls.txt

