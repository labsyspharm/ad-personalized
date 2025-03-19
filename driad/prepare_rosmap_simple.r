library(tidyverse)
library(synExtra)
library(data.table)
library(powerjoin)
library(here)
library(qs)

synapser::synLogin()
syn <- synDownloader("~/data", .cache = TRUE)

rosmap_quant_clinical <- syn("syn44335073") %>%
  read_csv()

# Use all ROSMAP patient data from PCC region
rosmap_quant_clinical_filtered <- rosmap_quant_clinical %>%
  # Log-transform all counts
  mutate(
    across(
      -c(ID, brain_region, PMI, AOD, CDR, Braak),
      \(x) log10(x + 1)
    )
  ) %>%
  filter(
    brain_region == "PCC"
  ) %>%
  select(
    -c(brain_region)
  )

fwrite(
  rosmap_quant_clinical_filtered,
  "driad/rosmap_quant_simple_clinical_filtered.tsv.gz",
  sep = "\t"
)

# Using `all_pooled` prediction task. Splits
# cases by Braak into early, mid, late
prediction_task <- DRIAD::prepareTask(
  "driad/rosmap_quant_simple_clinical_filtered.tsv.gz",
  task = "all_pooled"
) %>%
  mutate(
    across(
      Label,
      \(x) factor(
        as.character(x),
        levels = c("0", "1", "2"),
        ordered = TRUE
      )
    )
  )

prediction_task_pairs <- DRIAD::preparePairs(
  prediction_task
)

qsave(
  prediction_task_pairs, "driad/rosmap_quant_simple_pairs.qs"
)

dge_gmt <- map(
  c("syn20820056", "syn20820060"),
  ~syn(.x) %>%
    DRIAD::read_gmt()
) %>%
  enframe("experiment", "gene_sets") %>%
  mutate(
    gene_sets = map(gene_sets, ~enframe(.x, "drug", "genes"))
  ) %>%
  unnest(gene_sets)


drug_meta <- syn("syn11801537") %>%
  read_csv()

drugs_of_interest <- drug_meta %>%
  filter(
    name %in% c(
      "Ruxolitinib", "Baricitinib", "Tofacitinib", "Nilotinib", "Vorinostat", "NVP-TAE684"
    )
  )

dge_gmt_of_interest <- dge_gmt %>%
  inner_join(
    drugs_of_interest,
    by = c("drug" = "lincs_id")
  )

gmt_genes <- dge_gmt_of_interest$genes %>%
  reduce(union)
setdiff(gmt_genes, colnames(rosmap_quant_clinical))

gene_set_sizes <- map_int(dge_gmt_of_interest$genes, length) %>%
  unique()

gene_universe <- colnames(rosmap_quant_clinical) %>%
  setdiff(
    c("ID", "brain_region", "PMI", "AOD", "CDR",
      "Braak", "Barcode", "Label")
  )

# For each unique gene set cardinality, prepare
# 1000 background gene sets of the same size drawn
# randomly from the gene universe
set.seed(42)
background_gene_sets <- tibble(
  gene_set_size = gene_set_sizes
) %>%
  mutate(
    genes = map(
      gene_set_size,
      \(n) tibble(
        seq_id = 1:1000,
        genes = map(1:1000, \(y) sample(gene_universe, n, replace = FALSE))
      )
    )
  ) %>%
  unnest(genes)

# Try running DRIAD on a single gene set
sample_res <- DRIAD::evalGeneSets(
  background_gene_sets %>%
    head(1) %>%
    pull(genes),
  prediction_task,
  prediction_task_pairs,
  nBK = 0,
  method = "or"
)

# Try running on a bunch of example gene sets to estimate
# running time based on number of genes
set.seed(42)
running_time_jobs <- background_gene_sets %>%
  group_by(gene_set_size) %>%
  slice_sample(n = 2) %>%
  ungroup()

running_times_raw <- running_time_jobs %>%
  mutate(
    res = map(
      genes,
      \(x) {
        start <- proc.time()["elapsed"]
        res <- DRIAD::evalGeneSets(
          list(x),
          prediction_task,
          prediction_task_pairs,
          nBK = 0,
          method = "or"
        )
        end <- proc.time()["elapsed"]
        list(
          seconds_elapsed = end - start,
          res = res
        )
      }
    )
  )

qsave(
  running_times_raw,
  file.path("driad", "rosmap_simple_running_times_raw.qs")
)

running_times <- running_times_raw %>%
  unnest_wider(res) %>%
  select(gene_set_size, seconds_elapsed)

p <- ggplot(
  running_times,
  aes(gene_set_size, seconds_elapsed)
) +
  geom_point() +
  geom_smooth(method = "lm")

ggsave(
  file.path("driad_plots", "rosmap_simple_running_times.pdf"),
  p, width = 6, height = 4
)


DRIAD::evalGeneSets(
  background_gene_sets %>% filter(id == "prediction_task_11") %>% pull(genes) %>% head(n = 1),
  XY = prediction_tasks_pairs$task[[11]],
  lP = prediction_tasks_pairs$pairs[[11]],
  nBK = 0,
  method = "or"
) %>%
  dplyr::select(-Feats, -BK)
# 10 minutes for a single gene set in prediction task 11 with 512 samples

DRIAD::evalGeneSets(
  background_gene_sets %>% filter(id == "prediction_task_21") %>% pull(genes) %>% head(n = 1),
  XY = prediction_tasks_pairs$task[[21]],
  lP = prediction_tasks_pairs$pairs[[21]],
  nBK = 0,
  method = "or"
) %>%
  dplyr::select(-Feats, -BK)
# 2 minutes for a single gene set in prediction task 21 with 342 samples

DRIAD::evalGeneSets(
  background_gene_sets %>% filter(id == "prediction_task_19") %>% pull(genes) %>% head(n = 1),
  XY = prediction_tasks_pairs$task[[19]],
  lP = prediction_tasks_pairs$pairs[[19]],
  nBK = 0,
  method = "or"
) %>%
  dplyr::select(-Feats, -BK)
# 05:03 - 10:50 347s for a single gene set in prediction task 19 with 461 samples


DRIAD::evalGeneSets(
  background_gene_sets %>% filter(id == "prediction_task_23") %>% pull(genes) %>% head(n = 1),
  XY = prediction_tasks_pairs$task[[23]],
  lP = prediction_tasks_pairs$pairs[[23]],
  nBK = 0,
  method = "or"
) %>%
  dplyr::select(-Feats, -BK)
# 15:20 - 15:43 23s for a single gene set in prediction task 23 with 204 samples

DRIAD::evalGeneSets(
  background_gene_sets %>% filter(id == "prediction_task_3") %>% pull(genes) %>% head(n = 1),
  XY = prediction_tasks_pairs$task[[3]],
  lP = prediction_tasks_pairs$pairs[[3]],
  nBK = 0,
  method = "or"
) %>%
  dplyr::select(-Feats, -BK)
# 21:43 - 22:01 18s for a single gene set in prediction task 3 with 193 samples

DRIAD::evalGeneSets(
  background_gene_sets %>% filter(id == "prediction_task_1") %>% pull(genes) %>% head(n = 1),
  XY = prediction_tasks_pairs$task[[1]],
  lP = prediction_tasks_pairs$pairs[[1]],
  nBK = 0,
  method = "or"
) %>%
  dplyr::select(-Feats, -BK)
# 1s for a single gene set in prediction task 1 with 24 samples

running_times <- tibble(nsamples = c(512, 342, 461, 204, 193, 24), time_elapsed = c(10*60, 2*60, 347, 23, 18, 1)) %>%
  mutate(nsamples2 = nsamples**2, time_elapsed_sqrt = sqrt(time_elapsed), time_elapsed_crt = time_elapsed^(1/3))

qm <- lm(time_elapsed ~ nsamples + nsamples2, data = running_times)
summary(qm)
qm <- lm(time_elapsed_sqrt ~ nsamples, data = running_times)
summary(qm)
qm <- lm(time_elapsed_crt ~ nsamples, data = running_times)
summary(qm)

running_times %>%
  ggplot(aes(nsamples, time_elapsed)) +
    geom_point()

running_times %>%
  ggplot(aes(nsamples, time_elapsed_sqrt)) +
  geom_point() +
  geom_smooth(method = "lm")

running_times %>%
  ggplot(aes(nsamples, time_elapsed_crt)) +
  geom_point() +
  geom_smooth(method = "lm")

# Estimating how much time we need for each task
estimated_running_time <- prediction_tasks_pairs %>%
  transmute(
    id,
    n_samples = map_int(task, nrow),
    est_running_time = (qm$coefficients[["(Intercept)"]] + qm$coefficients[["nsamples"]] * n_samples)**3 %>%
      pmax(1) %>%
      magrittr::multiply_by(1.2)
  ) %>%
  distinct()

set.seed(42)
# Planning jobs so that each one should take ~1h
PLANNED_RUNNING_TIME <- 60*60
background_jobs <- background_gene_sets %>%
  inner_join(estimated_running_time, by = "id") %>%
  slice_sample(prop = 1) %>%
  group_by(id) %>%
  mutate(
    batch_id = floor(cumsum(est_running_time) / PLANNED_RUNNING_TIME)
  ) %>%
  ungroup() %>%
  group_nest(id, batch_id) %>%
  mutate(
    gene_sets = map(data, ~set_names(.x$genes, .x$gs_id)),
    path = paste0("driad/", id, ".qs")
  )

qsave(
  background_jobs, "driad/background_jobs.qs"
)
# background_jobs <- qread("driad/background_jobs.qs")


run_pred_job <- function(task_path, gene_sets, ...) {
  library(tidyverse)
  library(furrr)
  library(DRIAD)
  library(qs)
  library(here)

  plan(sequential)

  message("Reading task...", task_path)
  task <- qread(task_path)
  message("Evaluating gene sets...", length(gene_sets))
  out <- tryCatch({
    DRIAD::evalGeneSets(
      gene_sets,
      XY = task[["task"]],
      lP = task[["pairs"]],
      nBK = 0,
      method = "or"
    ) %>%
      dplyr::select(-Feats, -BK)
  },
  error = function(e) {
    warning(e)
    NULL
  }
  )
  message("Done...")
  out
}

library(batchtools)

reg <- makeRegistry(
  file.dir = here(paste0("registry_", gsub(" ", "_", Sys.time()))),
  seed = 1
)
# reg <- loadRegistry("reg  ")

batchMap(
  fun = run_pred_job,
  task_path = background_jobs$path,
  gene_sets = background_jobs$gene_sets
)

# run_bk_job(
#   task_id = prediction_tasks_gene_sets[1:5,][["id"]][[1]],
#   background_sets = prediction_tasks_gene_sets[1:5,][["background_sets"]][[1]],
#   background_task_id = prediction_tasks_gene_sets[1:5,][["background_task_id"]][[1]]
# )

job_table <- findJobs(reg = reg) %>%
  # Chunk jobs into a single array job
  mutate(chunk = 1)

submitJobs(
  job_table[findExpired()],
  resources = list(
    memory = "4gb",
    ncpus = 1L,
    partition = "short",
    walltime = 150*60,
    chunks.as.arrayjobs = TRUE,
    # For some reason these nodes fail to execute R because of an "illegal instruction"
    exclude = "compute-f-17-[09-25]"
  )
)

background_res <- reduceResultsList(
  job_table[findDone(reg = reg)], reg = reg
)

qsave(
  background_res, "driad/background_results_raw.qs"
)
# background_res <- qread("driad/background_results_raw.qs")
synStoreMany(
  "driad/background_results_raw.qs", "syn45053117", forceVersion = FALSE
)


compound_gene_sets <- prediction_tasks_pairs %>%
  select(-task, -pairs) %>%
  crossing(
    dge_gmt_of_interest
  ) %>%
  mutate(
    gs_id = paste0("compound_", name, "_", experiment, "_", id)
  )

set.seed(42)
compound_jobs <- compound_gene_sets %>%
  inner_join(estimated_running_time, by = "id") %>%
  slice_sample(prop = 1) %>%
  group_by(id) %>%
  mutate(
    batch_id = floor(cumsum(est_running_time) / PLANNED_RUNNING_TIME)
  ) %>%
  ungroup() %>%
  group_nest(id, batch_id) %>%
  mutate(
    gene_sets = map(data, ~set_names(.x$genes, .x$gs_id)),
    path = paste0("driad/", id, ".qs")
  )

reg_drugs <- makeRegistry(
  file.dir = here(paste0("registry_drugs_", gsub(" ", "_", Sys.time()))),
  seed = 1
)
reg_drugs <- loadRegistry("")

batchMap(
  fun = run_pred_job,
  task_path = compound_jobs$path,
  gene_sets = compound_jobs$gene_sets,
  reg = reg_drugs
)

job_table_drugs <- findJobs(reg = reg_drugs) %>%
  # Chunk jobs into a single array job
  mutate(chunk = 1)

submitJobs(
  job_table_drugs,
  resources = list(
    memory = "4gb",
    ncpus = 1L,
    partition = "short",
    walltime = 80*60,
    chunks.as.arrayjobs = TRUE,
    # For some reason these nodes fail to execute R because of an "illegal instruction"
    exclude = "compute-f-17-[09-25]"
  ),
  reg = reg_drugs
)

drug_res <- reduceResultsList(
  job_table_drugs[findDone(reg = reg_drugs)], reg = reg_drugs
)

qsave(
  drug_res, "driad/drug_results_raw.qs"
)

# drug_res <- qread("driad/drug_results_raw.qs")
synStoreMany(
  "driad/drug_results_raw.qs", "syn45053117", forceVersion = FALSE
)

# predictions <- background_gene_sets %>%
  # power_inner_join(
  #   rbindlist(background_res),
  #   by = c("gs_id" = "Set"),
  #   check = check_specs(
  #     duplicate_keys_left = "abort",
  #     duplicate_keys_right = "abort",
  #     unmatched_keys_left = "abort",
  #     unmatched_keys_right = "abort"
  #   )
  # ) %>%
#   bind_rows(
#     compound_gene_sets %>%
#       power_inner_join(
#         rbindlist(drug_res),
#         by = c("gs_id" = "Set"),
#         check = check_specs(
#           # duplicate_keys_left = "abort",
#           duplicate_keys_right = "abort",
#           unmatched_keys_left = "abort",
#           unmatched_keys_right = "abort"
#         )
#       )
#   )

compound_res_joined <- compound_gene_sets %>%
  power_inner_join(
    rbindlist(drug_res),
    by = c("gs_id" = "Set"),
    check = check_specs(
      # duplicate_keys_left = "abort",
      duplicate_keys_right = "abort",
      unmatched_keys_left = "abort",
      unmatched_keys_right = "abort"
    )
  )
fwrite(
  select(compound_res_joined, where(negate(is.list))),
  "driad/compound_res_raw.csv.gz"
)

background_res_joined <- background_gene_sets %>%
  power_inner_join(
    rbindlist(background_res),
    by = c("gs_id" = "Set"),
    check = check_specs(
      duplicate_keys_left = "abort",
      duplicate_keys_right = "abort",
      unmatched_keys_left = "abort",
      unmatched_keys_right = "abort"
    )
  )
fwrite(
  select(background_res_joined, where(negate(is.list))),
  "driad/background_res_raw.csv.gz"
)

compound_res_plotting <- compound_res_joined %>%
  bind_rows(
    background_res_joined %>%
      power_inner_join(
        compound_gene_sets %>%
          select(-gs_id, -link) %>%
          mutate(gene_set_size = map_int(genes, length)) %>%
          select(-genes),
        by = c("brain_region", "classification", "class", "comparison", "id", "gene_set_size"),
        check = check_specs(
          unmatched_keys_left = "abort",
          unmatched_keys_right = "abort"
        )
      )
  )

fwrite(
  select(compound_res_plotting, where(negate(is.list))),
  "driad/compound_and_background_res.csv.gz"
)

synStoreMany(
  c(
    "driad/compound_res_raw.csv.gz",
    "driad/background_res_raw.csv.gz",
    "driad/compound_and_background_res.csv.gz"
  ),
  "syn45053117", forceVersion = FALSE
)


plots <- compound_res_plotting %>%
  group_nest(
    classification, comparison
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
  function(classification, comparison, plot, ...) {
    ggsave(
      paste0("driad_plots/auc_", classification, "_", comparison, ".pdf"), plot,
      width = 7, height = 15
    )
  }
)


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


p <- compound_res_plotting %>%
  group_by(name) %>%
  mutate(
    name_experiment = if (length(unique(experiment)) > 1) paste0(name, " replicate ", experiment) else name
  ) %>%
  ungroup() %>%
  filter(
    name %in% c(
      "Baricitinib", "Ruxolitinib", "NVP-TAE684", "Tofacitinib"
    ),
    !str_detect(name_experiment, "replicate 2"),
    classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled"
  ) %>%
  mutate(name_experiment = str_replace(name_experiment, " replicate 1", "")) %>%
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
      filter(
        name %in% c(
          "Baricitinib", "Ruxolitinib", "NVP-TAE684", "Tofacitinib"
        ),
        !str_detect(name_experiment, "replicate 2"),
        classification == "high_threshold", brain_region == "PCC", comparison == "all_pooled"
      ) %>%
      mutate(name_experiment = str_replace(name_experiment, " replicate 1", "")) %>%
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
  "driad_plots/pcc_high_pooled_ridges_small.pdf", p, width = 5, height = 2
)
