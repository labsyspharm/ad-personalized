library(tidyverse)
library(synExtra)
library(synapser)
library(powerjoin)

# Select ROSMAP samples that use rRNA depletion instead of poly-A enrichment
# in library prep. Can't quantify repeat transcripts in samples prepared
# using poly-A enrichment

synLogin()

syn <- synDownloader("~/data", .cache = TRUE)

rosmap_clinical <- syn("syn3191087") %>%
  read_csv()

specimen_meta <- syn("syn21323366") %>%
  read_csv() %>%
  mutate(
    specimenID = str_replace(specimenID, "Sample_", "")
  ) %>%
  filter(
    assay == "rnaSeq"
  )

rna_meta <- syn("syn21088596") %>%
  read_csv() %>%
  mutate(
    across(specimenID, \(x) str_remove(x, fixed("Sample_")))
  ) %>%
  extract(
    notes, "library_prep_batch", r"{batch ([0-9]+)}",
    remove = FALSE, convert = TRUE
  )


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


# Only use batch 2 and 3. Batch 1 uses polyA enrichment, can't use for
# quantifying repetitive elements
# Or actually just filter directly on libraryPrep column
# All are paired end
rna_meta_eligible <- rna_meta %>%
  semi_join(fastq_df, by = c("specimenID" = "specimen_id")) %>%
  filter(
    libraryPrep == "rRNAdepletion",
    runType == "pairedEnd"
  )
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
  mutate(
    file_type = if_else(str_detect(file, "bam"), "bam", "fastq"),
    file_prefix = str_replace(file, r"{(\.fastq\.gz|\.bam)}", "") %>%
      str_replace("_R[12]_001", "")
  )


rosmap_file_to_clinical <- fastq_df_selected %>%
  distinct(specimen_id, file_prefix) %>%
  power_left_join(
    distinct(specimen_meta, individualID, specimenID, tissue),
    by = c("specimen_id" = "specimenID"),
    check = check_specs(
      unmatched_keys_left = "warn",
      duplicate_keys_right = "warn"
    )
  ) %>%
  mutate(
    brain_region = recode(
      tissue,
      `dorsolateral prefrontal cortex` = "DLPFC",
      `posterior cingulate cortex` = "PCC",
      `Head of caudate nucleus` = "HCN"
    )
  ) %>%
  power_inner_join(
    distinct(
      rosmap_clinical,
      projid, individualID,
      PMI = pmi, AOD = age_death, CDR = cogdx, Braak = braaksc
    ),
    by = "individualID",
    check = check_specs(
      unmatched_keys_left = "warn",
      duplicate_keys_right = "warn"
    )
  )


fastq_df_selected_brain_region <- fastq_df_selected %>%
  power_inner_join(
    rosmap_file_to_clinical %>%
      select(specimen_id, brain_region),
    by = "specimen_id",
    check = check_specs(
      unmatched_keys_left = "warn",
      duplicate_keys_right = "warn"
    )
  )

write_csv(
  fastq_df_selected_brain_region,
  "selected_rosmap_samples.csv"
)

write_csv(
  rosmap_file_to_clinical,
  "rosmap_file_to_clinical_mapping.csv.gz"
)

synStoreMany(
  c(
    "selected_rosmap_samples.csv",
    "rosmap_file_to_clinical_mapping.csv.gz"
  ),
  "syn43841077",
  forceVersion = FALSE
)
