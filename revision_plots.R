library(tidyverse)
library(synExtra)
library(data.table)
library(powerjoin)
library(here)
library(qs)
library(ComplexHeatmap)
library(ggbeeswarm)
library(broom)
library(ggpubr)

synapser::synLogin()
syn <- synDownloader("~/data", .cache = TRUE)

rosmap_tdp43_classification <- syn("syn44277761") %>%
  read_csv()

rosmap_quant_clinical <- syn("syn44335073") %>%
  read_csv()

rosmap_batch_meta <- syn("syn27034471") %>%
  read_csv()

rosmap_file_to_clinical_mapping <- syn("syn66532355") %>%
  read_csv()

# From Schoggins et al
is_genes <- syn("syn11629935") %>%
  read_lines()

rosmap_batch_annotation <- power_inner_join(
  rosmap_quant_clinical %>%
    select(ID),
  rosmap_file_to_clinical_mapping %>%
    distinct(
      ID = file_prefix, specimen_id, individualID
    ),
  by = "ID",
  check = check_specs(
    duplicate_keys_left = "warn",
    duplicate_keys_right = "warn",
    unmatched_keys_left = "warn"
  )
) %>%
  power_inner_join(
    rosmap_batch_meta %>%
      distinct(specimenID, batch = notes) %>%
      mutate(
        across(
          specimenID,
          \(x) str_remove(x, fixed("Sample_"))
        )
      ),
    by = c("specimen_id" = "specimenID"),
    check = check_specs(
      duplicate_keys_left = "warn",
      duplicate_keys_right = "warn",
      unmatched_keys_left = "warn"
    )
  )

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

createAnnotation <- function(df, colors = list(), which = "column") {
  # Load necessary libraries
  color_maps <- list()

  # Iterate over each column in the dataframe
  for (col_name in names(df)) {
    # Check if the column is numeric
    if (is.numeric(df[[col_name]])) {
      # Create a color mapping function for numeric columns
      if (is.null(colors[[col_name]]))
        colors[[col_name]] <- c("white", "red")
      if (is.function(colors[[col_name]])) {
        col_fun <- colors[[col_name]]
      } else if (is.character(colors[[col_name]])) {
        n <- length(colors[[col_name]])
        if (n == 1) {
          col_fun <- circlize::colorRamp2(
            seq(from = min(df[[col_name]]), to = max(df[[col_name]]), length.out = 100),
            viridis(100, option = colors[[col_name]])
          )
        } else {
          col_fun <- circlize::colorRamp2(
            seq(from = min(df[[col_name]]), to = max(df[[col_name]]), length.out = n),
            colors[[col_name]]
          )
        }
      } else
        stop("Dont know how to handle colors for column ", col_name)
      color_maps[[col_name]] <- col_fun
    } else {
      if (is.character(colors[[col_name]])) {
        color_maps[[col_name]] <- colors[[col_name]]
        next
      }
      col_fun <- colors[[col_name]] %||% \(n) ggokabeito::palette_okabe_ito(seq_len(n))
      # Create a named vector of colors for categorical columns
      values <- df[[col_name]]
      unique_values <- if (is.factor(values)) levels(values) else unique(df[[col_name]])
      named_colors <- setNames(col_fun(length(unique_values)), unique_values)
      color_maps[[col_name]] <- named_colors
    }
  }

  # Combine all annotations
  combined_annotations <- HeatmapAnnotation(df = df, col = color_maps, which = which)

  return(combined_annotations)
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
        power_inner_join(
          rosmap_batch_annotation %>%
            select(ID, batch),
          by = "ID",
          check = check_specs(
            duplicate_keys_left = "warn",
            duplicate_keys_right = "warn",
            unmatched_keys_left = "warn"
          )
        ) %>%
        column_to_rownames("ID")
    ),
    hm = Heatmap(
      mat,
      name = "Scaled log10(TPM + 1)",
      cluster_rows = cluster_fun_eucl,
      cluster_columns = cluster_fun_eucl,
      column_split = data$class,
      top_annotation = createAnnotation(
        col_meta %>%
          select(-class),
        which = "column"
      ),
      show_row_names = FALSE,
      show_column_names = FALSE
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


## Compare the expression of ISGs between the two classes

rosmap_quant_isg_data <- rosmap_quant_clinical_split %>%
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
    )
  ) %>%
  ungroup()


isg_expression_test_res <- rosmap_quant_isg_data %>%
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

isg_expression_test_res %>%
  select(brain_region, classification, mw_res) %>%
  unnest(mw_res)


ps <- rosmap_quant_isg_data %>%
  power_inner_join(
    isg_expression_test_res %>%
      select(brain_region, classification, mw_res) %>%
      unnest(mw_res),
    by = c("brain_region", "classification"),
    check = check_specs(
      duplicate_keys_left = "abort",
      duplicate_keys_right = "abort",
      unmatched_keys_left = "warn",
      unmatched_keys_right = "abort"
    )
  ) %>%
  rowwise() %>%
  mutate(
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
        geom_bracket(
          aes(
            xmin = cond_1, xmax = cond_2,
            label = signif(p.value, 2)
          ),
          y.position = 2,
          data = tibble(cond_1 = "positive", cond_2 = "negative", p.value = p.value),
          inherit.aes = FALSE
        ) +
        theme_minimal() +
        labs(
          x = "TDP-43 prediction",
          y = "Scaled log10(TPM + 1)",
          color = "Expressed CEs"
        )
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



# Check enrichment of genes enriched and depleted in the CRISPR screen

crispr_raw <- syn("syn66477979") %>%
  read_csv()

library(enrichR)

enrichr_raw <- enrichr(
  crispr_raw %>%
    arrange(desc(`Average LFC`)) %>%
    head(n = 100) %>%
    pull(`Gene Symbol`) %>%
    unique(),
  databases = "Reactome_Pathways_2024"
)

p <- enrichr_raw[[1]] %>%
  as_tibble() %>%
  mutate(
    Term = fct_reorder(
      Term, Odds.Ratio
    )
  ) %>%
  head(n = 5) %>%
  ggplot(
    aes(
      x = Odds.Ratio,
      y = Term
    )
  ) +
  geom_col() +
  geom_text(
    aes(
      label = cut(Adjusted.P.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns"))
    ),
    nudge_x = 3
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Enrichment Odds Ratio",
    y = NULL
  )

ggsave(
  here::here("plots", "crispr_enrichr_bars.pdf"),
  p, width = 5, height = 1.6
)
