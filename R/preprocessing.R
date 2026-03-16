
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