
# We first import the necessary libraries and the data.

library(tidyverse)
library(beemStatic) # remotes::install_github('lch14forever/beem-static')

# Load functions

source("functions.R")

# Load the abundance data

## IBD datasets
IBD_healthy_BEEM <- readRDS("IBD_healthy_BEEM_20260318.rds")
IBD_diseased_BEEM <- readRDS("IBD_diseased_BEEM_20260318.rds")

## IBS datasets
IBS_healthy_BEEM <- readRDS("IBS_healthy_BEEM_20260318.rds")
IBS_diseased_BEEM <- readRDS("IBS_diseased_BEEM_20260318.rds")

## CRC datasets
CRC_healthy_BEEM <- readRDS("CRC_healthy_BEEM_20260318.rds")
CRC_healthy_BEEM <- readRDS("CRC_diseased_BEEM_20260318.rds")

## CDI datasets (for now skip)

# Run BEEM-Static on each dataset

run_beem_static_bootstrap(IBD_healthy_BEEM, n_bootstrap = 2, subsample_pcts = c(0.4), ncpu = 4, max_iter = 2, output_file = "IBS_healthy_beem_bootstrap.csv")

# Load the results

res <- read.csv("IBS_healthy_beem_bootstrap.csv")

# Get the results

df_res <- get_beem_static_results(res)
