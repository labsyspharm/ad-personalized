library(tidyverse)
library(synExtra)
library(qs)
library(data.table)
library(powerjoin)
library(ggbeeswarm)

synapser::synLogin()
syn <- synDownloader("~/data", .cache = TRUE)

rosmap_quants <- syn("syn43841162") %>%
  fread()

msbb_quants <- syn("syn50920431") %>%
  fread()

msbb_clinical <- syn("syn6101474") %>%
  read_csv()

msbb_specimen_meta <- syn("syn21893059") %>%
  read_csv()


rosmap_key_mapping <- syn("syn3382527") %>%
  read_csv()

rosmap_clinical <- syn("syn3191087") %>%
  read_csv()

rosmap_specimen_meta <- syn("syn21323366") %>%
  read_csv() %>%
  mutate(
    specimenID = str_replace(specimenID, "Sample_", "")
  )
# Parsing issues in Brodmann are col. Not an issue

rosmap_file_meta <- syn("syn44137214") %>%
  fread() %>%
  mutate(
    file_prefix = str_replace(file, r"{(\.fastq\.gz|\.bam)}", "") %>%
      str_replace("_R[12]_001", "")
  )

rosmap_file_to_clinical <- rosmap_file_meta %>%
  left_join(
    distinct(rosmap_specimen_meta, individualID, specimenID, tissue),
    by = c("specimen_id" = "specimenID")
  ) %>%
  mutate(
    brain_region = recode(
      tissue,
      `dorsolateral prefrontal cortex` = "DLPFC",
      `posterior cingulate cortex` = "PCC",
      `Head of caudate nucleus` = "HCN"
    )
  )

rosmap_rna_meta <- syn("syn21088596") %>%
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


transcripts_of_interest <- tribble(
  ~transcript_id, ~gene, ~variant_association, ~transcript_name,
  "STMN2short", "STMN2", "TDP-43-", "STMN2short",
  "UNC13A-CE1", "UNC13A", "TDP-43-", "UNC13A-CE1",
  "UNC13A-CE2", "UNC13A", "TDP-43-", "UNC13A-CE2",
  "ENST00000220876.12", "STMN2", "TDP-43+", "STMN2\ncanonical",
  "ENST00000519716.7", "UNC13A", "TDP-43+", "UNC13A\ncanonical"
) %>%
  mutate(across(transcript_name, fct_inorder))


marker_counts_msbb <- msbb_quants %>%
  inner_join(
    transcripts_of_interest,
    by = c("Name" = "transcript_id")
  ) %>%
  inner_join(
    distinct(msbb_specimen_meta, individualID, specimenID, brain_region = tissue),
    by = c("file" = "specimenID")
  )


marker_counts_rosmap <- rosmap_quants %>%
  inner_join(
    transcripts_of_interest,
    by = c("Name" = "transcript_id")
  ) %>%
  inner_join(
    distinct(rosmap_rna_meta, file_prefix, library_prep_batch, specimenID),
    by = c("file" = "file_prefix"),
    checks = check_specs(
      unmatched_keys_left = "abort",
      duplicate_keys_right = "abort"
    )
  ) %>%
  inner_join(
    distinct(rosmap_file_to_clinical, specimen_id, brain_region),
    by = c("specimenID" = "specimen_id")
  )

marker_counts <- list(
  rosmap = marker_counts_rosmap,
  msbb = marker_counts_msbb
) %>%
  rbindlist(idcol = "dataset", use.names = TRUE, fill = TRUE)

p <- ggplot(
  marker_counts %>%
    mutate(
      TPM = if_else(TPM == 0, 1e-2, TPM),
      NumReads = if_else(NumReads == 0, 1e-2, NumReads)
    ) %>%
    filter(!(dataset == "msbb" & brain_region == "prefrontal cortex")) %>%
    pivot_longer(
      c(TPM, NumReads),
      names_to = "metric", values_to = "count"
    ),
  aes(metric, count, group = file)
) +
  geom_line(color = "gray", alpha = 0.5) +
  geom_quasirandom() +
  facet_grid(vars(transcript_name), vars(dataset, brain_region)) +
  scale_y_log10() +
  theme_minimal()

ps <- marker_counts %>%
  mutate(
    TPM = if_else(TPM == 0, 1e-2, TPM),
    NumReads = if_el
  pivot_longer(
    c(TPM, NumReads),
    names_to = "metric", values_to = "count"
  ) %>%
  group_by(Name, gene, variant_association, transcript_name) %>%
  summarize(
    plot = list(
      ggplot(
        cur_data(),
        aes(metric, count, group = file)
      ) +
        geom_line(color = "gray", alpha = 0.5) +
        geom_quasirandom() +
        facet_wrap(vars(dataset, brain_region)) +
        scale_y_log10() +
        theme_minimal()
    ),
    .groups = "drop"
  )

pwalk(
  ps,
  function(Name, transcript_name, plot, ...) {
    ggsave(
      paste0("driad_plots/tpm_vs_num_", transcript_name, ".pdf"), plot,
      width = 10, height = 10
    )
  }
)

ps <- map(
  c("TPM", "NumReads") %>%
    set_names(),
  ~ggplot(
    marker_counts %>%
      mutate(
        TPM = if_else(TPM == 0, 0.5 * min(TPM[TPM > 0]), TPM),
        NumReads = if_else(NumReads == 0, 0.5 * min(NumReads[NumReads > 0]), NumReads)
      ) %>%
      filter(
        !(dataset == "rosmap" & brain_region != "PCC")
      ),
    aes(paste(dataset, brain_region, sep = "\n"), !!rlang::sym(.x))
  ) +
    geom_quasirandom() +
    facet_wrap(vars(transcript_name)) +
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = .5))
)



iwalk(
  ps,
  function(p, metric) {
    ggsave(
      paste0("driad_plots/tpm_rosmap_pcc_vs_msbb_", metric, ".pdf"), p,
      width = 10, height = 10
    )
  }
)

marker_counts %>%
  mutate(zero_count = TPM == 0) %>%
  group_by(transcript_name, dataset, brain_region) %>%
  summarize(
    n_zero = sum(zero_count),
    frac_zero = n_zero / n(),
    .groups = "drop"
  )

p <- marker_counts %>%
  mutate(
    zero_count = factor(
      TPM == 0,
      c(TRUE, FALSE), labels = c("0", ">0")
    ),
    dataset = str_to_upper(dataset)
  ) %>%
  filter(
    variant_association == "TDP-43-",
    !(dataset == "rosmap" & brain_region != "PCC"),
    brain_region != "prefrontal cortex"
  ) %>%
  ggplot(
    aes(paste(dataset, brain_region, sep = "\n"), fill = zero_count)
  ) +
  geom_bar(position = position_fill()) +
  geom_text(
    aes(label = after_stat(count), y = after_stat(prop), group = zero_count),
    stat = "count",
    position = position_fill(vjust = 0)
  ) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~transcript_name) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(
    x = NULL,
    y = "Proportion of samples",
    fill = "Samples with\ntranscript count"
  )

ggsave(
  paste0("driad_plots/zero_counts_rosmap_pcc_vs_msbb.pdf"), p,
  width = 9, height = 5
)
