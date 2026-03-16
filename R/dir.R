

data_dir <- function(){
  if (str_detect(getwd(), "maria"))
    data_dir <- "/Users/maria/OneDrive - UCL/Bureau/Mariam/2025-2026/03_others/FWO/04_projects/00_eLetter_science/00_data"
  else
    stop("You need to specify the path to the data directory in `reanalysis-article-2026/R/dir.R`")
  data_dir
}


result_dir <- function(){
  if (str_detect(getwd(), "maria"))
    data_dir <- "/Users/maria/OneDrive - UCL/Bureau/Mariam/2025-2026/03_others/FWO/04_projects/00_eLetter_science/01_results"
  else
    stop("You need to specify the path to the data directory in `reanalysis-article-2026/R/dir.R`")
  data_dir
}