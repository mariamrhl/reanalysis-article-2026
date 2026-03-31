
subsample_matrix <- function(matrix, subsample_pct) {
  
  n_samples <- ncol(matrix)
  subsample_idx <- sample(seq_len(n_samples), size = floor(n_samples * subsample_pct), replace = FALSE)
  matrix_sub <- matrix[, subsample_idx, drop = FALSE]
  
  return(matrix_sub)
  
}

subsample_unique_matrix <- function(matrix, n_samples, pct, seen_hashes) {
  
  matrix_sub <- subsample_matrix(matrix, pct)
  
  subset_hash <- digest::digest(colnames(matrix_sub))
  
  if (!exists(subset_hash, envir = seen_hashes)) {
    assign(subset_hash, TRUE, envir = seen_hashes)
    return(matrix_sub)
  }
  
  return(NULL)
}

filter_prevalence_matrix <- function(matrix, prevalence_threshold = 0.2, label = "Other") {
  
  tryCatch({
    # prevalence per taxon
    prevalence <- rowMeans(matrix > 0)
    taxa_to_keep <- names(prevalence[prevalence >= prevalence_threshold])
    
    cat("Keeping", length(taxa_to_keep), "out of", nrow(matrix), "taxa with prevalence >=", prevalence_threshold * 100, "%\n")
    
    # taxa failing threshold
    taxa_drop <- setdiff(rownames(matrix), taxa_to_keep)
    if (length(taxa_drop) == 0) return(matrix)
    
    # collapse rare taxa into "Other"
    other_row <- colSums(matrix[taxa_drop, , drop = FALSE])
    
    matrix_filtered <- matrix[taxa_to_keep, , drop = FALSE]
    matrix_filtered <- rbind(matrix, label = other_row)
    
    return(matrix_filtered)}, 
    
    error = function(e) {cat("Filtering failed:", conditionMessage(e), "\n")
      return(NULL)})
}

get_network_df <- function(network, pct, iter, n_taxa, metrics = NULL) {
  network |>
    as.data.frame() |>
    rownames_to_column("taxon_from") |>
    pivot_longer(-taxon_from, names_to = "taxon_to", values_to = "weight") |>
    mutate(fraction = pct, iteration = iter, n_taxa = n_taxa)
}

get_run_metadata <- function(tool, package_name = tool) {
  tibble(
    tool = tool, tool_version = as.character(packageVersion(package_name)), 
    r_version = R.version$version.string, date = Sys.Date()
  )
}

compute_diversity_metrics <- function(matrix) {
  # vegan expects samples as rows, taxa as columns
  matrix_t <- t(matrix)
  
  shannon  <- vegan::diversity(matrix_t, index = "shannon")
  richness <- rowSums(matrix_t > 0)  # samples as rows now
  sheldon  <- exp(shannon) / richness
  pielou   <- shannon / log(richness)
  
  tibble(sample = rownames(matrix_t), shannon = shannon, richness = richness, sheldon = sheldon, pielou = pielou)
}


compute_rho <- function(weights) {
  
  pos <- weights[weights > 0]
  neg <- weights[weights < 0]
  
  sum_pos <- sum(pos)
  sum_neg <- sum(neg)
  
  (sum_pos + sum_neg) / (sum_pos - sum_neg)
}


get_stats_df <- function(bootstrap_results, symmetric = FALSE) {
  
  results_split <- bootstrap_results |>
    group_by(fraction, iteration) |>
    group_split()
  
  metrics <- lapply(results_split, function(df) {
    
    b_est <- df |>
      select(taxon_from, taxon_to, weight) |>
      pivot_wider(names_from = taxon_to, values_from = weight) |>
      column_to_rownames("taxon_from") |>
      as.matrix()
    
    # For directed tools (BEEM-Static, LIMITS): use both triangles
    # For undirected tools (FlashWeave, SPIEC-EASI): use upper triangle only
    # rho ratio cancels double-counting so results are comparable
    values <- if (symmetric) {
      b_est[upper.tri(b_est)]
    } else {
      c(b_est[upper.tri(b_est)], b_est[lower.tri(b_est)])
    }
    
    pep <- seqtime::getPep(b_est)
    edge_num <- sum(values != 0)
    rho <- compute_rho(values)
    
    tibble(fraction = unique(df$fraction), iteration = unique(df$iteration), pep = pep, edge_num = edge_num, rho = rho, 
           shannon = unique(df$shannon), richness = unique(df$richness), sheldon = unique(df$sheldon), pielou = unique(df$pielou))
  })
  
  bind_rows(metrics)
}


