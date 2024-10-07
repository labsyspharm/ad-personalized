library(tidyverse)
library(synExtra)
library(synapser)

# Select ROSMAP samples that use rRNA depletion instead of poly-A enrichment
# in library prep. Can't quantify repeat transcripts in samples prepared
# using poly-A enrichment

synLogin()

syn <- synDownloader("~/data", .cache = TRUE)

syns_fastq_is <- synChildren("syn21589959")
syns_fastq_cs <- synChildren("syn8612097")
# All samples from batch 2 are only available as bam
syns_bam_is <- synChildren("syn21188662")

fastq_df <- bind_rows(
  syns_bam_is %>%
    enframe("file", "syn_id") %>%
    bind_cols(
      str_replace(
        .$file,
        fixed(".final"), ""
      ) %>%
        str_match(
          r"{(.*)\.bam}"
        ) %>%
          magrittr::set_colnames(
            c("file", "specimen_id")
          ) %>%
          as_tibble() %>%
          select(-file)
    ),
  syns_fastq_is %>%
    enframe("file", "syn_id") %>%
    bind_cols(
      str_match(
        .$file,
        r"{(.+_[0-9]+(?:_redo|_rerun)?)_?(S[0-9]+)_R([12])_001\.fastq\.gz}"
      ) %>%
        magrittr::set_colnames(
          c("file", "specimen_id", "run_sample_id", "read_number")
        ) %>%
        as_tibble() %>%
        select(-file)
    ),
  syns_fastq_cs %>%
    enframe("file", "syn_id") %>%
    bind_cols(
      # str_replace(
      #   .$file,
      #   fixed("Sample_"), ""
      # ) %>%
      str_match(
        .$file,
        r"{(.*)\.r([12])\.fastq\.gz}"
      ) %>%
        magrittr::set_colnames(
          c("file", "specimen_id", "read_number")
        ) %>%
        as_tibble() %>%
        select(-file)
    )
)

rna_meta <- syn("syn21088596") %>%
  read_csv() %>%
  mutate(
    across(specimenID, str_replace, fixed("Sample_"), "")
  ) %>%
  extract(
    notes, "library_prep_batch", r"{batch ([0-9]+)}",
    remove = FALSE, convert = TRUE
  )

specimen_meta <- syn("syn21323366") %>%
  read_csv()
# Parsing issues in Brodmann are col. Not an issue

# Only use batch 2 and 3. Batch 1 uses polyA enrichment, can't use for
# quantifying repetitive elements
# Or actually just filter directly on libraryPrep column
rna_meta_eligible <- rna_meta %>%
  semi_join(fastq_df, by = c("specimenID" = "specimen_id")) %>%
  filter(libraryPrep == "rRNAdepletion")
  # filter(library_prep_batch %in% c(2, 3))

# set.seed(42)
# rna_meta_selected <- rna_meta_eligible %>%
#   slice_sample(n = 100)

all(rna_meta_eligible$specimenID %in% fastq_df$specimen_id)

fastq_df_selected <- fastq_df %>%
  semi_join(
    rna_meta_eligible,
    by = c("specimen_id" = "specimenID")
  ) %>%
  mutate(file_type = if_else(str_detect(file, "bam"), "bam", "fastq"))

write_csv(
  fastq_df_selected,
  "selected_rosmap_samples.csv"
)
