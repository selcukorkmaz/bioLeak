## Install missing dependencies for simulation study
if (!requireNamespace("PRROC", quietly = TRUE))
  install.packages("PRROC", repos = "https://cran.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cran.r-project.org")

if (!requireNamespace("curatedOvarianData", quietly = TRUE))
  BiocManager::install("curatedOvarianData", ask = FALSE, update = FALSE)

cat("All dependencies installed.\n")
