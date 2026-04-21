
## Importing packages

library(tidyverse)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(beemStatic) # remotes::install_github('lch14forever/beem-static')
library(seqtime) # remotes::install_github("hallucigenia-sparsa/seqtime")
library(phyloseq)


## Load the data

mae <- readRDS("MAE_processed_20260329.rds")

## Import functions

lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source)

## Get the abundance matrix for each dataset

### IBD datasets

IBD_long <- mae[["IBD"]] |>
  assay() |>
  t() |>
  as.data.frame() |>
  rownames_to_column("uid") |>
  pivot_longer(-uid, names_to = "taxon", values_to = "rel_ab") |>
  left_join(
    mae[["IBD"]] |> colData() |> as_tibble(rownames = "uid") |> select(uid, pid, visit_num, diagnosis),
    by = join_by(uid)
  )

IBD_healthy_long <- IBD_long |> filter(diagnosis == "nonIBD")
IBD_diseased_long <- IBD_long |> filter(diagnosis == "UC" | diagnosis == "CD")

IBD_healthy <- IBD_healthy_long |>
  select(uid, taxon, rel_ab) |>
  pivot_wider(names_from = taxon, values_from = rel_ab) |>
  column_to_rownames("uid") |>
  as.matrix() |>
  t()

IBD_diseased <- IBD_diseased_long |>
  select(uid, taxon, rel_ab) |>
  pivot_wider(names_from = taxon, values_from = rel_ab) |>
  column_to_rownames("uid") |>
  as.matrix() |>
  t()

rm(IBD_long, IBD_healthy_long, IBD_diseased_long)

### IBS dataset

IBS_long <- mae[["IBS"]] |> 
  assay() |> 
  t() |> 
  as.data.frame() |> 
  rownames_to_column("uid") |> 
  pivot_longer(-uid, names_to = "taxon", values_to = "count") |>
  group_by(uid) |>
  mutate(rel_ab = count/sum(count)) |>
  left_join(
    mae[["IBS"]] |> colData() |> as_tibble(rownames = "uid") |> select(uid, pid, time, cohort), 
    by = join_by(uid)
  )

IBS_healthy_long <- IBS_long |> filter(cohort == "H")
IBS_diseased_long <- IBS_long |> filter(cohort == "C" | cohort == "D")


IBS_healthy <- IBS_healthy_long |>
  select(uid, taxon, rel_ab) |>
  pivot_wider(names_from = taxon, values_from = rel_ab) |>
  column_to_rownames("uid") |>
  as.matrix() |>
  t()

IBS_diseased <- IBS_diseased_long |>
  select(uid, taxon, rel_ab) |>
  pivot_wider(names_from = taxon, values_from = rel_ab) |>
  column_to_rownames("uid") |>
  as.matrix() |>
  t()

rm(IBS_long, IBS_healthy_long, IBS_diseased_long)

### CRC dataset

CRC_long <- mae[["CRC"]] |>
  assay() |>
  t() |>
  as.data.frame() |>
  rownames_to_column("uid") |>
  pivot_longer(-uid, names_to = "taxon", values_to = "count") |>
  group_by(uid) |>
  mutate(rel_ab = count/sum(count)) |>
  left_join(
    mae[["CRC"]] |> colData() |> as_tibble(rownames = "uid") |> select(uid, cohort),
    by = join_by(uid)
  )

CRC_healthy_long <- CRC_long |> filter(cohort == "H")
CRC_diseased_long <- CRC_long |> filter(cohort == "CRC")

CRC_healthy <- CRC_healthy_long |>
  select(uid, taxon, rel_ab) |>
  pivot_wider(names_from = taxon, values_from = rel_ab) |>
  column_to_rownames("uid") |>
  as.matrix() |>
  t()

CRC_diseased <- CRC_diseased_long |>
  select(uid, taxon, rel_ab) |>
  pivot_wider(names_from = taxon, values_from = rel_ab) |>
  column_to_rownames("uid") |>
  as.matrix() |>
  t()

rm(CRC_long, CRC_healthy_long, CRC_diseased_long)

### CDI dataset

CDI_long <- mae[["CDI"]] |>
  assay() |>
  t() |>
  as.data.frame() |>
  rownames_to_column("uid") |>
  pivot_longer(-uid, names_to = "taxon", values_to = "count") |>
  group_by(uid) |>
  mutate(rel_ab = count/sum(count)) |>
  left_join(
    mae[["CDI"]] |> colData() |> as_tibble(rownames = "uid") |> select(uid, status, subject_disease_status),
    by = join_by(uid)
  ) |>
  filter(subject_disease_status == "CDI" | subject_disease_status == "CTR")

CDI_healthy_long <- CDI_long |> filter(subject_disease_status == "CTR")
CDI_diseased_long <- CDI_long |> filter(subject_disease_status == "CDI")

CDI_healthy <- CDI_healthy_long |>
  select(uid, taxon, rel_ab) |>
  pivot_wider(names_from = taxon, values_from = rel_ab) |>
  column_to_rownames("uid") |>
  as.matrix() |>
  t()

CDI_diseased <- CDI_diseased_long |>
  select(uid, taxon, rel_ab) |>
  pivot_wider(names_from = taxon, values_from = rel_ab) |>
  column_to_rownames("uid") |>
  as.matrix() |>
  t()

rm(CDI_long, CDI_healthy_long, CDI_diseased_long)

## Run LIMITS

### IBD datasets
# n_bagging = 200 as per the article's implementation

limits_IBD_healthy <- run_network_bootstrap(IBD_healthy, run_fn = run_limits, n_bootstrap = 100, prevalence_threshold = 1, max_attempts = 100, bagging.iter = 200)
limits_IBD_diseased <- run_network_bootstrap(IBD_diseased, run_fn = run_limits, n_bootstrap = 100, prevalence_threshold = 1, max_attempts = 100, bagging.iter = 200)

## Run BEEM-static

### IBD datasets
beem_IBD_healthy <- run_network_bootstrap(IBD_healthy, run_fn = run_beem_static, n_bootstrap = 100, prevalence_threshold = 0.3, max_attempts = 100)
beem_IBD_diseased <- run_network_bootstrap(IBD_diseased, run_fn = run_beem_static, n_bootstrap = 100, prevalence_threshold = 0.3, max_attempts = 100)

### IBS datasets
beem_IBS_healthy <- run_network_bootstrap(IBS_healthy, run_fn = run_beem_static, n_bootstrap = 100, prevalence_threshold = 0.3, max_attempts = 100)
beem_IBS_diseased <- run_network_bootstrap(IBS_diseased, run_fn = run_beem_static, n_bootstrap = 100, prevalence_threshold = 0.3, max_attempts = 100)

### CRC datasets
beem_CRC_healthy <- run_network_bootstrap(CRC_healthy, run_fn = run_beem_static, n_bootstrap = 100, prevalence_threshold = 0.3, max_attempts = 100)
beem_CRC_diseased <- run_network_bootstrap(CRC_diseased, run_fn = run_beem_static, n_bootstrap = 100, prevalence_threshold = 0.3, max_attempts = 100)

### CDI datasets
beem_CDI_healthy <- run_network_bootstrap(CDI_healthy, run_fn = run_beem_static, n_bootstrap = 100, prevalence_threshold = 0.3, max_attempts = 100)
beem_CDI_diseased <- run_network_bootstrap(CDI_diseased, run_fn = run_beem_static, n_bootstrap = 100, prevalence_threshold = 0.3, max_attempts = 100)

## Save results

saveRDS(list(
  limits_IBD_healthy = limits_IBD_healthy,
  limits_IBD_diseased = limits_IBD_diseased,
  beem_IBD_healthy = beem_IBD_healthy,
  beem_IBD_diseased = beem_IBD_diseased,
  beem_IBS_healthy = beem_IBS_healthy,
  beem_IBS_diseased = beem_IBS_diseased,
  beem_CRC_healthy = beem_CRC_healthy,
  beem_CRC_diseased = beem_CRC_diseased,
  beem_CDI_healthy = beem_CDI_healthy,
  beem_CDI_diseased = beem_CDI_diseased
), "network_bootstrap_results.rds")

