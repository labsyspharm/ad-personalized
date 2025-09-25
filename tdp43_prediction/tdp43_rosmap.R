library(tidyverse)
library(synExtra)
library(data.table)
library(powerjoin)

synapser::synLogin()
syn <- synDownloader("~/data", .cache = TRUE)

rosmap_quant <- syn("syn43841162") %>%
  fread()

rosmap_key_mapping <- syn("syn3382527") %>%
  read_csv()

rosmap_clinical <- syn("syn3191087") %>%
  read_csv()

specimen_meta <- syn("syn21323366") %>%
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



transcripts_of_interest <- tribble(
  ~transcript_id, ~gene, ~variant_association, ~transcript_name,
  "STMN2short", "STMN2", "TDP-43-", "STMN2short",
  "UNC13A-CE1", "UNC13A", "TDP-43-", "UNC13A-CE1",
  "UNC13A-CE2", "UNC13A", "TDP-43-", "UNC13A-CE2",
  "ENST00000220876.12", "STMN2", "TDP-43+", "STMN2\ncanonical",
  "ENST00000519716.7", "UNC13A", "TDP-43+", "UNC13A\ncanonical"
) %>%
  mutate(across(transcript_name, fct_inorder))

marker_counts_rosmap <- rosmap_quant %>%
  inner_join(
    transcripts_of_interest,
    by = c("Name" = "transcript_id")
  ) %>%
  power_inner_join(
    distinct(rna_meta, file_prefix, library_prep_batch, specimenID) %>%
      drop_na(file_prefix),
    by = c("file" = "file_prefix"),
    check = check_specs(
      unmatched_keys_left = "abort",
      duplicate_keys_right = "abort"
    )
  ) %>%
  inner_join(
    distinct(rosmap_file_to_clinical, specimen_id, brain_region),
    by = c("specimenID" = "specimen_id")
  )

toi_thresholds <- tribble(
  ~transcript_name, ~brain_region, ~threshold,
  "STMN2short", "HCN", 1.3,
  "UNC13A-CE1", "HCN", 0.09,
  "UNC13A-CE2", "HCN", 0.08,
  "STMN2short", "DLPFC", 1.3,
  "UNC13A-CE1", "DLPFC", 0.09,
  "UNC13A-CE2", "DLPFC", 0.08,
  "STMN2short", "PCC", 1.5,
  "UNC13A-CE1", "PCC", 0.09,
  "UNC13A-CE2", "PCC", 0.08
)

rosmap_classification_expression <- marker_counts_rosmap %>%
  drop_na(brain_region) %>%
  inner_join(
    toi_thresholds,
    by = c("transcript_name", "brain_region")
  ) %>%
  mutate(
    expressed = TPM > threshold
  ) %>%
  select(file, specimenID, transcript_name, brain_region, expressed, TPM, threshold)

rosmap_classification <- rosmap_classification_expression %>%
  group_by(file, specimenID, brain_region) %>%
  summarize(n_expressed = sum(expressed), .groups = "drop") %>%
  inner_join(
    rosmap_classification_expression %>%
      select(-TPM, -threshold) %>%
      pivot_wider(names_from = transcript_name, values_from = expressed),
    by = c("file", "specimenID", "brain_region")
  ) %>%
  mutate(
    class_low = if_else(n_expressed > 0, "positive", "negative"),
    class_high = if_else(n_expressed > 1, "positive", "negative"),
  )



rosmap_classification %>%
  count(brain_region, n_expressed)

rosmap_classification %>%
  count(brain_region, class_low)

rosmap_classification %>%
  count(brain_region, class_high)

rosmap_classification %>%
  select(STMN2short, `UNC13A-CE1`, `UNC13A-CE2`) %>%
  table() %>%
  fisher.test()

write_csv(
  rosmap_classification,
  "tdp43_prediction/rosmap_classification.csv.gz"
)
synStoreMany(
  "tdp43_prediction/rosmap_classification.csv.gz", "syn43841077", forceVersion = FALSE
)


p <- rosmap_classification %>%
  mutate(
    n_expressed_class = factor(
      n_expressed, levels = c(0, 1, 2, 3), labels = c("none", "1", "2", "3"),
      ordered = TRUE
    )
  ) %>%
  group_by(brain_region, n_expressed_class) %>%
  summarize(n = n(), .groups = "drop") %>%
  ggplot(aes(brain_region, y = n, fill = n_expressed_class)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.5) +
    # shadowtext::geom_shadowtext(aes(label = n), position = position_stack(vjust = 0.5), size = 2) +
    labs(x = "Tissue", y = "Number of samples", fill = "Number of cryptic\ntranscripts expressed\nabove threshold")

ggsave(
  "plots/rosmap_classification_n_samples.pdf", width = 4, height = 3
)ß

p <- rosmap_classification %>%
  mutate(
    n_expressed_class = factor(
      n_expressed, levels = c(0, 1, 2, 3), labels = c("none", "1", "2", "3"),
      ordered = TRUE
    )
  ) %>%
  filter(brain_region == "PCC") %>%
  group_by(brain_region, n_expressed_class) %>%
  summarize(n = n(), .groups = "drop") %>%
  ggplot(aes(n_expressed_class, y = n, fill = n_expressed_class)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  guides(fill = "none") +
  # shadowtext::geom_shadowtext(aes(label = n), position = position_stack(vjust = 0.5), size = 2) +
  labs(x = "Number of cryptic transcripts\nexpressed above threshold", y = "Number of samples", fill = NULL)

ggsave(
  "plots/rosmap_classification_n_samples_pcc_only.pdf", width = 3, height = 2
)


p <- marker_counts_rosmap %>%
  # mutate(
  #   TPM_adj = TPM + .5 * min(salmon_quants_long$TPM[salmon_quants_long$TPM > 0])
  # ) %>%
  mutate(across(transcript_name, fct_rev)) %>%
  ggplot(aes(transcript_name, TPM)) +
  geom_quasirandom(
    # aes(shape = outlier_bottom, color = outlier_bottom),
    aes(shape = outlier_bottom),
    data = ~.x %>%
      mutate(
        outlier_bottom = TPM == 0,
        TPM = if_else(TPM > 0, TPM, 0.5 * min(TPM[TPM > 0]))
      )
  ) +
  scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = 16), guide = "none") +
  # scale_color_manual(values = c(`TRUE` = "gray40", `FALSE` = "black"), guide = "none") +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = pos),
    alpha = 0.2,
    data = function(count_df) {
      toi_thresholds %>%
        mutate(
          transcript_name = factor(transcript_name, levels = levels(count_df$transcript_name)),
          xmin = as.integer(transcript_name) - 0.5,
          xmax = as.integer(transcript_name) + 0.5
        ) %>%
        crossing(pos = c("positive", "negative")) %>%
        mutate(
          ymin = if_else(pos == "positive", threshold, min(count_df$TPM)),
          ymax = if_else(pos == "positive", max(count_df$TPM), threshold)
        )
    },
    inherit.aes = FALSE
  ) +
  # scale_fill_manual(values = c("positive" = "firebrick1", "negative" = "mediumblue")) +
  scale_fill_manual(values = c("positive" = "red", "negative" = "blue")) +
  scale_y_log10() +
  facet_wrap(~brain_region) +
  theme_minimal() +
  coord_flip() +
  labs(x = "Transcript", fill = "TDP-43 pathology\nprediction thresholds")

ggsave(
  "plots/rosmap_toi_counts_colored.pdf", p,
  width = 8, height = 4
)



p <- marker_counts_rosmap %>%
  # mutate(
  #   TPM_adj = TPM + .5 * min(salmon_quants_long$TPM[salmon_quants_long$TPM > 0])
  # ) %>%
  filter(brain_region == "PCC") %>%
  mutate(across(transcript_name, fct_rev)) %>%
  ggplot(aes(transcript_name, TPM)) +
  geom_quasirandom(
    # aes(shape = outlier_bottom, color = outlier_bottom),
    aes(shape = outlier_bottom),
    data = ~.x %>%
      mutate(
        outlier_bottom = TPM == 0,
        TPM = if_else(TPM > 0, TPM, 0.5 * min(TPM[TPM > 0]))
      )
  ) +
  scale_shape_manual(values = c(`TRUE` = 6, `FALSE` = 16), guide = "none") +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = pos),
    alpha = 0.2,
    data = function(count_df) {
      toi_thresholds %>%
        filter(brain_region == "PCC") %>%
        mutate(
          transcript_name = factor(transcript_name, levels = levels(count_df$transcript_name)),
          xmin = as.integer(transcript_name) - 0.5,
          xmax = as.integer(transcript_name) + 0.5
        ) %>%
        crossing(pos = c("positive", "negative")) %>%
        mutate(
          ymin = if_else(pos == "positive", threshold, min(count_df$TPM)),
          ymax = if_else(pos == "positive", max(count_df$TPM), threshold)
        )
    },
    inherit.aes = FALSE
  ) +
  # geom_segment(
  #   aes(x = as.integer(transcript_name) - 0.3, xend = as.integer(transcript_name) + 0.3, y = threshold, yend = threshold),
  #   # alpha = 0.2,
  #   color = "red",
  #   data = function(count_df) {
  #     toi_thresholds %>%
  #       mutate(
  #         transcript_name = factor(
  #           transcript_name, levels = levels(count_df$transcript_name)
  #         )
  #       ) %>%
  #       filter(brain_region == "PCC")
  #   },
  #   inherit.aes = FALSE
  # ) +
  scale_fill_manual(values = c("positive" = "firebrick1", "negative" = "mediumblue")) +
  # scale_fill_manual(values = c("positive" = "red", "negative" = "blue")) +
  scale_y_log10() +
  theme_minimal() +
  labs(x = "Transcript", fill = "TDP-43 pathology\nprediction thresholds")

ggsave(
  "plots/rosmap_toi_counts_pcc_only.pdf", p,
  width = 7, height = 3
)


 toi_all_vs_all <- combn(transcripts_of_interest$transcript_id, 2, simplify = FALSE) %>%
  map(
    function(tids) {
      replacements <- set_names(c("transcript_id_1", "transcript_id_2"), tids)
      marker_counts_rosmap %>%
        filter(Name %in% tids) %>%
        select(file, Name, TPM, brain_region) %>%
        mutate(
          Name = recode(Name, !!!replacements),
          comparison = paste(tids, collapse = " - ")
        )
    }
  ) %>%
  bind_rows() %>%
  pivot_wider(names_from = Name, values_from = TPM)

p <- toi_all_vs_all %>%
  mutate(
    across(
      starts_with("transcript_id_"),
      ~if_else(.x == 0, min(.x[.x > 0]) / 2, .x) %>%
        log10()
    )
  ) %>%
  # mutate(
  #   TPM_adj = TPM + .5 * min(salmon_quants_long$TPM[salmon_quants_long$TPM > 0])
  # ) %>%
  ggplot(aes(transcript_id_1, transcript_id_2)) +
  geom_point() +
  facet_wrap(~comparison, scales = "free")
# geom_quasirandom(
#   aes(shape = outlier_bottom),
#   data = ~.x %>%
#     mutate(
#       outlier_bottom = TPM == 0,
#       TPM = if_else(TPM > 0, TPM, 0.5 * min(TPM[TPM > 0]))
#     )
# ) +
# scale_shape_manual(values = c(`TRUE` = 25, `FALSE` = 16), guide = "none") +
# facet_wrap(~variant_association, scale = "free_x") +
# scale_y_log10()

toi_all_vs_all %>%
  mutate(
    across(
      starts_with("transcript_id_"),
      ~.x > 0
    )
  ) %>%
  group_by(brain_region, comparison) %>%
  summarize(cont_table = table(tibble(transcript_id_1, transcript_id_2)) %>% list(), .groups = "drop") %>%
  mutate(
    obs_exp = map(cont_table, possibly(chisq.test, NULL)) %>%
      map("stdres"),
    test_res = map(cont_table, possibly(fisher.test, NULL)) %>%
      map(possibly(broom::tidy, NULL))
  ) %>%
  unnest(test_res) %>%
  View()

rosmap_classification_braak <- rosmap_classification  %>%
  inner_join(
    rosmap_file_to_clinical %>%
      distinct(file_prefix, individualID),
    by = c("file" = "file_prefix")
  ) %>%
  inner_join(
    rosmap_clinical %>%
      distinct(individualID, braaksc)
  )

library(ggrepel)

p <- rosmap_classification_braak %>%
  count(brain_region, class_high, braaksc) %>%
  ggplot(aes(class_high, y = n, fill = fct_inseq(as.character(braaksc), ordered = TRUE))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text_repel(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = 1,
    direction = "y",
    segment.color = NA
  ) +
  facet_wrap(~brain_region) +
  labs(x = "TDP-43 pathology prediction", y = "N cases", fill = "Braak stage")

ggsave(
  "plots/tdp43_prediction_high_braak_stage_distribution.pdf", width = 8, height = 4
)


p <- rosmap_classification_braak %>%
  filter(brain_region == "PCC") %>%
  count(class_high, braaksc) %>%
  ggplot(aes(class_high, y = n, fill = fct_inseq(as.character(braaksc), ordered = TRUE))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text_repel(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = 1,
    direction = "y",
    segment.color = NA
  ) +
  # facet_wrap(~brain_region) +
  labs(x = "TDP-43 pathology prediction", y = "N cases", fill = "Braak stage")

ggsave(
  "plots/tdp43_prediction_high_pcc_braak_stage_distribution.pdf", width = 6, height = 3
)
