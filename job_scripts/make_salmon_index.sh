# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

wget https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

wget https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

grep "^>" <(gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) | cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

# Including custom splice variants of STMN2 and UNC13A

cat <(gunzip -c Homo_sapiens.GRCh38.cdna.all.fa.gz) \
  stmn2-short.fa unc13a-ce.fa \
  <(gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) | \
  gzip > \
  Homo_sapiens.GRCh38.gentrome_including_variants.fa.gz

salmon index -t Homo_sapiens.GRCh38.gentrome_including_variants.fa.gz \
  -d decoys.txt -p 12 \
  -i Homo_sapiens.GRCh38.gentrome_including_variants_index

