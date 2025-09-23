# install_packages.R
# This script installs all required CRAN and Bioconductor packages for static_threebody_sim.R
# Written using CoPilot
options(repos = c(CRAN = "https://ftp.osuosl.org/pub/cran/")) #Change to whatever local mirror makes sense


# Install pak and BiocManager if not already installed
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# List of CRAN packages
cran_packages <- c(
  "nosoi",
  "tidyverse",
  "viridis",
  "igraph",
  "ggraph",
  "ape",
  "tidytree",
  "argparse"
)

# List of Bioconductor packages
bioc_packages <- c(
  "treeio",
  "ggtree"
)

# Install CRAN packages in parallel using pak
pak::pkg_install(cran_packages)

# Install Bioconductor packages
BiocManager::install(bioc_packages)
