library(tidyverse)
library(synExtra)
library(data.table)
library(powerjoin)
library(here)
library(qs)
library(ComplexHeatmap)
library(ggbeeswarm)
library(broom)

synapser::synLogin()
syn <- synDownloader("~/data", .cache = TRUE)

rosmap_tdp43_classification <- syn("syn44277761") %>%
  read_csv()

rosmap_quant_clinical <- syn("syn44335073") %>%
  read_csv()

# From Schoggins et al
is_genes <- syn("syn11629935") %>%
  read_lines()

rosmap_quant_clinical_split <- rosmap_quant_clinical %>%
  mutate(
    across(-c(ID, brain_region, PMI, AOD, CDR, Braak), ~log10(.x + 1))
  ) %>%
  power_inner_join(
    distinct(rosmap_tdp43_classification, file, class_low, class_high, n_expressed),
    by = c("ID" = "file"),
    check = check_specs(
      duplicate_keys_left = "abort",
      duplicate_keys_right = "abort",
      unmatched_keys_left = "warn",
      unmatched_keys_right = "abort"
    )
  ) %>% {
    bind_rows(
      select(., -class_high) %>%
        rename(class = class_low) %>%
        mutate(classification = "low_threshold"),
      select(., -class_low) %>%
        rename(class = class_high) %>%
        mutate(classification = "high_threshold"),
    )
  } %>%
  group_nest(brain_region, classification)


rosmap_quant_clinical_split$data[[5]] %>%
  transmute(x = str_sub(ID, 1, 4)) %>%
  count(x) %>%
  arrange(desc(n))

library(seriation)
cluster_fun_eucl <- function(mat, sample_in_col = TRUE) {
  # if (!sample_in_col) {
  #   mat <- t(mat)
  # }
  # mat_imp <- impute.knn(
  #   mat, rng.seed = 42
  # )[["data"]]
  # if (!sample_in_col) {
  #   mat_imp <- t(mat_imp)
  #   mat <- t(mat)
  # }
  # browser()
  dist_mat <- dist(mat)
  # dist_mat <- as.dist(mat)
  clust <- hclust(dist_mat, method = "average")
  reorder(clust, dist_mat, method = "OLO")
}


hms <- rosmap_quant_clinical_split %>%
  rowwise() %>%
  mutate(
    mat = list(
      select(data, -c(PMI, AOD, CDR, Braak, .data[["class"]], n_expressed)) %>%
        column_to_rownames("ID") %>%
        as.matrix() %>% {
          .[, intersect(is_genes, colnames(.))]
        } %>%
        scale() %>%
        t()
    ),
    col_meta = list(
      data %>%
        select(ID, .data[["class"]], n_expressed) %>%
        column_to_rownames("ID")
    ),
    hm = Heatmap(
      mat,
      name = "Scaled log10(TPM + 1)",
      cluster_rows = cluster_fun_eucl,
      cluster_columns = cluster_fun_eucl,
      column_split = data$class,
      top_annotation = HeatmapAnnotation(
        df = col_meta,
        which = "column"
      )
    ) %>%
      list()
  ) %>%
  ungroup()

hms$hm[[5]]

pwalk(
  hms,
  \(hm, brain_region, classification, ...) {
    withr::with_pdf(
      here::here("plots", glue::glue("isg_heatmap_{brain_region}_{classification}.pdf")),
      draw(hm),
      width = 12, height = 10
    )
  }
)

## Scatter plot


ps <- rosmap_quant_clinical_split %>%
  rowwise() %>%
  mutate(
    p_data = list(
      select(data, -c(PMI, AOD, CDR, Braak), any_of(is_genes)) %>%
        pivot_longer(
          any_of(is_genes),
          names_to = "gene",
          values_to = "expression"
        ) %>%
        group_by(gene) %>%
        mutate(
          expression = scale(expression)[, 1]
        ) %>%
        ungroup() %>%
        group_by(ID, `class`, n_expressed) %>%
        summarize(
          across(expression, mean),
          .groups = "drop"
        )
    ),
    p = list(
      p_data %>%
        mutate(
          across(n_expressed, \(x) factor(as.character(x), levels = as.character(0:3), ordered = TRUE))
        ) %>%
        ggplot(
          aes(
            x = `class`,
            y = expression,
            color = n_expressed
          )
        ) +
        geom_quasirandom(shape = 16, alpha = .8) +
        theme_minimal()
    )
  ) %>%
  ungroup()

pwalk(
  ps,
  \(p, brain_region, classification, ...) {
    ggsave(
      here::here("plots", glue::glue("isg_beeswarm_{brain_region}_{classification}.pdf")),
      p, width = 6, height = 4
    )
  }
)

isg_expression_test_res <- ps %>%
  mutate(
    mw_raw = map(
      p_data,
      \(x) wilcox.test(
        expression ~ `class`,
        data = mutate(x, across(`class`, \(y) factor(y, levels = c("positive", "negative")))),
        conf.int = TRUE
      )
    ),
    mw_res = map(mw_raw, tidy)
  )

isg_expression_test_res$mw_raw[[5]]

#
library(ggalluvial)

rosmap_classification <- syn("syn44277761") %>%
  read_csv()

ptx_rescue_type_list %>%
  count(perturbed)

rosmap_classification_alluvials <- rosmap_classification %>%
  pivot_longer(
    starts_with("class_"),
    names_to = "class_threshold",
    values_to = "class"
  ) %>%
  group_nest(brain_region, class_threshold) %>%
  mutate(
    alluvial_data = map(
      data,
      \(x) x %>%
        mutate(
          across(
            c(STMN2short, `UNC13A-CE1`, `UNC13A-CE2`),
            \(y) if_else(y, "detected", "not detected")
          )
        ) %>%
        count(
          STMN2short, `UNC13A-CE1`, `UNC13A-CE2`, class
        )
    ),
    lode_data = map(
      alluvial_data,
      \(x) to_lodes_form(
        x,
        -c(n, class),
        key = "Transcript", value = "Sample"
      )
    ),
    p = map(
      lode_data,
      \(x)
      ggplot(
        x,
        aes(
          x = Transcript, y = n, stratum = `Sample`, alluvium = alluvium
        )
      ) +
        scale_x_discrete(expand = c(0.05, 0.05)) +
        geom_alluvium(
          aes(
            fill = class
          ),
          width = 1/10
        ) +
        ggokabeito::scale_fill_okabe_ito(order = c(2, 1)) +
        labs(x = NULL) +
        geom_stratum(width = 1/6) +
        geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
        geom_text(
          aes(
            x = stage(Transcript, after_scale = x + .12),
            label = n
          ),
          stat = StatAlluvium,
          hjust = 0,
          data = \(x) filter(x, as.integer(Transcript) == 1L)
        ) +
        theme_minimal() +
        # labs(y = "Data points", fill = NA) +
        theme(
            panel.grid.major.x = element_blank()
        ) +
        coord_cartesian(clip = "off")
    )
  )

pwalk(
  rosmap_classification_alluvials,
  \(p, brain_region, class_threshold,...) ggsave(
    file.path("plots", paste0("cryptic_exon_expression_overlap_alluvial_", brain_region, "_", class_threshold, ".pdf")),
    p, width = 6, height = 4
  )
)
