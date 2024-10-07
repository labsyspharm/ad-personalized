library(tidyverse)
library(synExtra)
library(data.table)
library(tximport)
library(qs)
library(powerjoin)

synapser::synLogin()
syn <- synDownloader("~/data", .cache = TRUE)

rosmap_quants_raw <- syn("syn43841162") %>%
  fread()

rosmap_clinical <- syn("syn3191087") %>%
  read_csv()

specimen_meta <- syn("syn21323366") %>%
  read_csv() %>%
  mutate(
    specimenID = str_replace(specimenID, "Sample_", "")
  )

rosmap_file_meta <- syn("syn44137214") %>%
  fread() %>%
  mutate(
    file_prefix = str_replace(file, r"{(\.fastq\.gz|\.bam)}", "") %>%
      str_replace("_R[12]_001", "")
  )

rna_meta <- syn("syn21088596") %>%
  read_csv() %>%
  mutate(
    across(specimenID, str_replace, fixed("Sample_"), "")
  ) %>%
  extract(
    notes, "library_prep_batch", r"{batch ([0-9]+)}",
    remove = FALSE, convert = TRUE
  ) %>%
  left_join(
    select(rosmap_file_meta, specimen_id, file_prefix),
    by = c("specimenID" = "specimen_id")
  )

rosmap_file_to_clinical <- rosmap_file_meta %>%
  distinct(specimen_id, file_prefix) %>%
  left_join(
    distinct(specimen_meta, individualID, specimenID, tissue),
    by = c("specimen_id" = "specimenID")
  ) %>%
  mutate(
    brain_region = recode(
      tissue,
      `dorsolateral prefrontal cortex` = "DLPFC",
      `posterior cingulate cortex` = "PCC",
      `Head of caudate nucleus` = "HCN"
    )
  ) %>%
  inner_join(
    distinct(
      rosmap_clinical,
      projid, individualID,
      PMI = pmi, AOD = age_death, CDR = cogdx, Braak = braaksc
    ),
    by = "individualID"
  )

genmap <- syn("syn44325472") %>%
  read_csv()

lf <- Sys.glob("quants/*/quant.sf") %>%
  set_names(str_replace(., "quants/", "") %>% str_replace("/quant.sf", ""))

tx2gene <- genmap %>%
  filter(external_gene_source == "HGNC Symbol") %>%
  distinct(TXNAME = ensembl_transcript_id, GENEID = external_gene_name)

rosmap_counts_raw <- tximport(lf, tx2gene = tx2gene, type = "salmon", ignoreTxVersion = TRUE)

qsave(rosmap_counts_raw, "rosmap_tximport.qs")
synStoreMany("rosmap_tximport.qs", "syn43841077")

rosmap_abundance_protein <- rosmap_counts_raw$abundance[
  genmap %>%
    filter(gene_biotype == "protein_coding", external_gene_source == "HGNC Symbol") %>%
    pull(external_gene_name) %>%
    unique() %>%
    intersect(rownames(rosmap_counts_raw$abundance)) %>%
    sort(),
] %>%
  t() %>%
  as_tibble(rownames = "ID")

rosmap_abundance_clinical <- rosmap_abundance_protein %>%
  power_left_join(
    distinct(rosmap_file_to_clinical, file_prefix, brain_region, PMI, AOD, CDR, Braak),
    by = c("ID" = "file_prefix"),
    check = check_specs(
      duplicate_keys_left = "abort",
      duplicate_keys_right = "abort",
      unmatched_keys_left = "abort"
    )
  ) %>%
  select(ID, brain_region, PMI, AOD, CDR, Braak, everything())

fwrite(
  rosmap_abundance_clinical,
  "rosmap_quant_clinical.csv.gz"
)
synStoreMany("rosmap_quant_clinical.csv.gz", "syn43841077")
