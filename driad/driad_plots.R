library(tidyverse)
library(synExtra)
library(data.table)
library(powerjoin)
library(here)
library(qs)

synapser::synLogin()
syn <- synDownloader("~/data", .cache = TRUE)

rosmap_res <- syn("syn45053249") %>%
  read_csv()

msbb_res <- syn("syn50934769") %>%
  read_csv()

combined_res <- list(
  rosmap = rosmap_res,
  msbb = msbb_res
) %>%
  rbindlist(use.names = TRUE, idcol = "dataset")

plots <- combined_res %>%
  group_nest(
    dataset, classification, comparison
  ) %>%
  mutate(
    plot = map(
      data,
      function(d) {
        d %>%
          ggplot(
            aes(AUC)
          ) +
          geom_density(
            data = ~.x %>%
              filter(str_starts(gs_id, "background"))
          ) +
          geom_vline(
            aes(xintercept = AUC),
            data = ~.x %>%
              filter(str_starts(gs_id, "compound"))
          ) +
          facet_grid(
            rows = vars(name, experiment, class),
            cols = vars(brain_region)
          )
      }
    )
  )

dir.create("driad_plots")
pwalk(
  plots,
  function(dataset, classification, comparison, plot, ...) {
    ggsave(
      paste0("driad_plots/auc_", dataset, "_", classification, "_", comparison, ".pdf"), plot,
      width = 7, height = 15
    )
  }
)

combined_res_pvals <- combined_res %>%
  mutate(
    type = if_else(str_starts(gs_id, "compound"), "compound", "background"),
  ) %>%
  group_by(
    dataset, classification, comparison,
    name, experiment, class, brain_region
  ) %>%
  summarize(
    n_drug = sum(type == "compound"),
    pval = (sum(AUC[type == "compound"] < AUC[type == "background"]) + 1) / n(),
    .groups = "drop"
  ) %>%
  filter(n_drug == 1)

ps <- combined_res_pvals %>%
  group_by(
    classification, comparison,
    dataset, brain_region, class, name
  ) %>%
  mutate(
    name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
  ) %>%
  ungroup() %>%
  group_by(
    classification, comparison
  ) %>%
  summarize(
    plot = list(
      ggplot(
        cur_data(),
        aes(
          class, name_experiment, fill = -log10(pval)
        )
      ) +
        geom_tile() +
        facet_wrap(~dataset + brain_region) +
        scale_fill_viridis_c()
    ),
    .groups = "drop"
  )

pwalk(
  ps,
  function(classification, comparison, plot, ...) {
    ggsave(
      paste0("driad_plots/p_heatmaps_", classification, "_", comparison, ".pdf"),
      plot, width = 10, height = 10
    )
  }
)


ps <- combined_res_pvals %>%
  group_by(
    classification, comparison,
    dataset, brain_region, class, name
  ) %>%
  mutate(
    name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
  ) %>%
  ungroup() %>%
  group_by(
    classification, comparison, dataset
  ) %>%
  summarize(
    plot = list(
      ggplot(
        cur_data(),
        aes(
          class, name_experiment, fill = -log10(pval)
        )
      ) +
        geom_tile() +
        facet_wrap(~brain_region) +
        scale_fill_viridis_c()
    ),
    .groups = "drop"
  )

pwalk(
  ps,
  function(classification, comparison, plot, dataset, ...) {
    ggsave(
      paste0("driad_plots/p_heatmaps_", dataset, "_", classification, "_", comparison, ".pdf"),
      plot, width = 6, height = 4
    )
  }
)


ps <- combined_res_pvals %>%
  group_by(
    classification, comparison,
    dataset, brain_region, class, name
  ) %>%
  mutate(
    name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
  ) %>%
  ungroup() %>%
  group_by(
    classification, comparison
  ) %>%
  summarize(
    plot = list(
      ggplot(
        cur_data(),
        aes(
          class, -log10(pval), group = 1
        )
      ) +
        geom_point() +
        geom_line() +
        facet_grid(vars(dataset, brain_region), vars(name_experiment))
    ),
    .groups = "drop"
  )

pwalk(
  ps,
  function(classification, comparison, plot, ...) {
    ggsave(
      paste0("driad_plots/p_line_", classification, "_", comparison, ".pdf"),
      plot, width = 10, height = 10
    )
  }
)



####
####
####


p <- compound_res_plotting %>%
  filter(classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled") %>%
  group_by(name) %>%
  mutate(
    name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
  ) %>%
  ungroup() %>%
  mutate(
    name_experiment = factor(name_experiment, levels = sort(unique(name_experiment)))
  ) %>%
  ggplot(
    aes(AUC)
  ) +
  geom_density(
    data = ~.x %>%
      filter(str_starts(gs_id, "background"))
  ) +
  geom_vline(
    aes(xintercept = AUC),
    data = ~.x %>%
      filter(str_starts(gs_id, "compound"))
  ) +
  facet_grid(vars(name_experiment), vars(class))

ggsave(
  "driad_plots/pcc_high_pooled_densities.pdf", p, width = 8, height = 12
)

compound_res_p_values <- compound_res_plotting %>%
  mutate(
    type = if_else(str_starts(gs_id, "compound"), "compound", "background"),
  ) %>%
  group_by(brain_region, classification, class, comparison, id, experiment, drug, name, target_name) %>%
  summarize(
    n_drug = sum(type == "compound"),
    pval = (sum(AUC[type == "compound"] < AUC[type == "background"]) + 1) / n(),
    .groups = "drop"
  )

compound_res_p_values_comparison <- compound_res_p_values %>%
  select(brain_region, classification, comparison, class, experiment, drug, name, target_name, pval) %>%
  mutate(across(pval, ~-log10(.x))) %>%
  pivot_wider(names_from = class, values_from = pval)

p <- compound_res_p_values_comparison %>%
  ggplot(aes(negative, positive)) +
  geom_point() +
  facet_grid(vars(brain_region), vars(classification, comparison)) +
  coord_equal()

ggsave(
  "driad_plots/pvalues_positive_vs_negative.pdf", p, width = 12, height = 12
)

library(ggridges)
ggplot(diamonds, aes(x = price, y = cut)) +
  geom_density_ridges(scale = 4) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges()

p <- compound_res_plotting %>%
  filter(classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled") %>%
  group_by(name) %>%
  mutate(
    name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
  ) %>%
  ungroup() %>%
  mutate(
    name_experiment = factor(name_experiment, levels = sort(unique(name_experiment))) %>%
      fct_rev()
  ) %>%
  ggplot(
    aes(AUC, name_experiment)
  ) +
  geom_density_ridges2(
    data = ~.x %>%
      filter(str_starts(gs_id, "background")),
    alpha = 0.8
  ) +
  geom_segment(
    aes(x = AUC, xend = AUC, y = as.numeric(name_experiment), yend = as.numeric(name_experiment) + 0.9),
    data = ~.x %>%
      filter(str_starts(gs_id, "compound")),
    color = "firebrick1"
  ) +
  geom_text(
    aes(x = 0.2, y = as.numeric(name_experiment) + 0.7, label = label),
    data = compound_res_p_values %>%
      filter(classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled") %>%
      group_by(name) %>%
      mutate(
        name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
      ) %>%
      ungroup() %>%
      mutate(
        name_experiment = factor(name_experiment, levels = sort(unique(name_experiment))) %>%
          fct_rev()
      ) %>%
      mutate(
        label = paste0("p=", signif(pval, 2))
      ),
    inherit.aes = FALSE, hjust = 0, vjust = 0.5
  ) +
  facet_wrap(vars(class), ncol = 2) +
  theme_minimal() +
  labs(y = "")

ggsave(
  "driad_plots/pcc_high_pooled_ridges.pdf", p, width = 5, height = 4
)


p <- combined_res %>%
  group_by(name) %>%
  mutate(
    name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
  ) %>%
  ungroup() %>%
  filter(
    name %in% c(
      "Baricitinib", "Ruxolitinib", "Tofacitinib"
    ),
    # !str_detect(name_experiment, "replicate 2"),
    classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled"
  ) %>%
  # mutate(name_experiment = str_replace(name_experiment, " replicate 1", "")) %>%
  mutate(
    name_experiment = factor(name_experiment, levels = sort(unique(name_experiment))) %>%
      fct_rev()
  ) %>%
  ggplot(
    aes(AUC, name_experiment)
  ) +
  geom_density_ridges2(
    data = ~.x %>%
      filter(str_starts(gs_id, "background")),
    alpha = 0.8
  ) +
  geom_segment(
    aes(x = AUC, xend = AUC, y = as.numeric(name_experiment), yend = as.numeric(name_experiment) + 0.9),
    data = ~.x %>%
      filter(str_starts(gs_id, "compound")),
    color = "firebrick1"
  ) +
  geom_text(
    aes(x = 0.2, y = as.numeric(name_experiment) + 0.7, label = label),
    data = combined_res_pvals %>%
      filter(classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled") %>%
      group_by(name) %>%
      mutate(
        name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
      ) %>%
      ungroup() %>%
      filter(
        name %in% c(
          "Baricitinib", "Ruxolitinib", "Tofacitinib"
        ),
        # !str_detect(name_experiment, "replicate 2"),
        classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled"
      ) %>%
      # mutate(name_experiment = str_replace(name_experiment, " replicate 1", "")) %>%
      mutate(
        name_experiment = factor(name_experiment, levels = sort(unique(name_experiment))) %>%
          fct_rev()
      ) %>%
      mutate(
        label = paste0("p=", signif(pval, 2))
      ),
    inherit.aes = FALSE, hjust = 0, vjust = 0.5
  ) +
  facet_wrap(vars(class), ncol = 2) +
  theme_minimal() +
  labs(y = "") +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0))

ggsave(
  "driad_plots/pcc_high_pooled_ridges_small.pdf", p, width = 5, height = 2.2
)


p <- combined_res %>%
  group_by(name) %>%
  mutate(
    name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
  ) %>%
  ungroup() %>%
  filter(
    name %in% c(
      "Baricitinib", "Ruxolitinib", "Tofacitinib"
    ),
    !str_detect(name_experiment, "replicate 2"),
    classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled"
  ) %>%
  # mutate(name_experiment = str_replace(name_experiment, " replicate 1", "")) %>%
  mutate(
    name_experiment = factor(name, levels = sort(unique(name))) %>%
      fct_rev()
  ) %>%
  ggplot(
    aes(AUC, name)
  ) +
  geom_density_ridges2(
    data = ~.x %>%
      filter(str_starts(gs_id, "background")),
    alpha = 0.8
  ) +
  geom_segment(
    aes(x = AUC, xend = AUC, y = as.numeric(name), yend = as.numeric(name) + 0.9),
    data = ~.x %>%
      filter(str_starts(gs_id, "compound")),
    color = "firebrick1"
  ) +
  geom_text(
    aes(x = 0.2, y = as.numeric(name) + 0.7, label = label),
    data = combined_res_pvals %>%
      filter(classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled") %>%
      group_by(name) %>%
      mutate(
        name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
      ) %>%
      ungroup() %>%
      filter(
        name %in% c(
          "Baricitinib", "Ruxolitinib", "Tofacitinib"
        ),
        !str_detect(name_experiment, "replicate 2"),
        classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled"
      ) %>%
      # mutate(name_experiment = str_replace(name_experiment, " replicate 1", "")) %>%
      mutate(
        name = factor(name, levels = sort(unique(name))) %>%
          fct_rev()
      ) %>%
      mutate(
        label = paste0("p=", signif(pval, 2))
      ),
    inherit.aes = FALSE, hjust = 0, vjust = 0.5
  ) +
  facet_wrap(vars(class), ncol = 2) +
  theme_minimal() +
  labs(y = "") +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0))

ggsave(
  "driad_plots/pcc_high_pooled_ridges_small.pdf", p, width = 5, height = 2.2
)
