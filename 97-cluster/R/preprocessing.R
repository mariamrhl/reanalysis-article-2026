
phyloseq_to_SE <- function(ps) {
  counts <- as.matrix(otu_table(ps))
  
  if (!taxa_are_rows(ps)) {
    counts <- t(counts)
  }
  
  col_data <- sample_data(ps) |> as("data.frame") |> as.data.frame()
  
  row_data <- tax_table(ps) |> as("matrix") |> as.data.frame()
  
  SummarizedExperiment(
    assays  = list(counts = counts),
    colData = col_data,
    rowData = row_data
  )
}

filter_prevalence <- function(df, prevalence_threshold = 0.2) {
  
  # Get taxa passing prevalence threshold
  taxa_to_keep <- df |>
    group_by(taxon) |>
    summarise(prevalence = mean(rel_ab > 0), .groups = "drop") |>
    filter(prevalence >= prevalence_threshold) |>
    pull(taxon)
  
  cat("Keeping", length(taxa_to_keep), "out of", 
      df$taxon |> unique() |> length(), 
      "taxa with prevalence >=", prevalence_threshold * 100, "%\n")
  
  df |>
    mutate(taxon = ifelse(taxon %in% taxa_to_keep, taxon, "Other")) |>
    group_by(uid, taxon) |>
    summarise(rel_ab = sum(rel_ab), .groups = "drop")
  
}

filter_prevalence_matrix <- function(abundance_mat, prevalence_threshold = 0.2, other_label = "Other") {
  
  # prevalence per taxon
  prevalence <- rowMeans(abundance_mat > 0)
  
  taxa_to_keep <- names(prevalence[prevalence >= prevalence_threshold])
  
  cat(
    "Keeping", length(taxa_to_keep), "out of",
    nrow(abundance_mat),
    "taxa with prevalence >=",
    prevalence_threshold * 100, "%\n"
  )
  
  # taxa failing threshold
  taxa_drop <- setdiff(rownames(abundance_mat), taxa_to_keep)
  
  if (length(taxa_drop) == 0)
    return(abundance_mat)
  
  # collapse rare taxa into "Other"
  other_row <- colSums(abundance_mat[taxa_drop, , drop = FALSE])
  
  abundance_filtered <- abundance_mat[taxa_to_keep, , drop = FALSE]
  
  abundance_filtered <- rbind(
    abundance_filtered,
    Other = other_row
  )
  
  return(abundance_filtered)
}


filter_abundance <- function(df, abundance_threshold = 0.001) {
  
  taxa_to_keep <- df |>
    group_by(taxon) |>
    summarise(mean_ab = mean(rel_ab), .groups = "drop") |>
    filter(mean_ab >= abundance_threshold) |>
    pull(taxon)
  
  cat("Keeping", length(taxa_to_keep), "out of", 
      df$taxon |> unique() |> length(), 
      "taxa with mean abundance >=", abundance_threshold, "\n")
  
  df |>
    mutate(taxon = ifelse(taxon %in% taxa_to_keep, taxon, "Other")) |>
    group_by(uid, taxon) |>
    summarise(rel_ab = sum(rel_ab), .groups = "drop")
}