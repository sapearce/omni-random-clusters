#!/usr/bin/env Rscript
# Random clustering baseline - assigns every cell to a random cluster

library(argparse)
library(jsonlite)
library(zellkonverter)

# command-line arguments in same style as the OSCA module
parser <- ArgumentParser(description = "Random clusters baseline")
parser$add_argument("--output_dir", "-o", dest = "output_dir", type = "character", required = TRUE)
parser$add_argument("--name", "-n", dest = "name", type = "character", required = TRUE)
parser$add_argument("--data.h5ad", dest = "data_path", type = "character", required = TRUE)
parser$add_argument("--data.clusters_truth_num", dest = "clusters_truth_num_path", type = "character", required = TRUE)
parser$add_argument("--seed", dest = "seed", type = "integer", default = 42, required = FALSE)
args <- parser$parse_args()

# output files are the same four files the metrics stage expects
timings_path         <- file.path(args$output_dir, paste0(args$name, ".timings.json"))
clustering_info_path <- file.path(args$output_dir, paste0(args$name, ".clustering_info.json"))
clusters_path        <- file.path(args$output_dir, paste0(args$name, ".clusters.tsv"))
pca_path             <- file.path(args$output_dir, paste0(args$name, ".pca.tsv"))

# read the data in -  cell IDs and the true number of clusters
sce <- readH5AD(args$data_path, reader = "python")
cell_ids  <- colnames(sce)
n_cells   <- length(cell_ids)
n_cluster <- as.integer(readLines(args$clusters_truth_num_path, n = 1))

# random clusters - each cell gets a cluster from 1 to n_cluster (reproducible via the seed)
set.seed(args$seed)
start_time <- Sys.time()
random_clusters <- sample(1:n_cluster, size = n_cells, replace = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

# clusters file - the same random labels in the louvain and leiden columns
write.table(
  data.frame(cell_id = cell_ids, louvain = random_clusters, leiden = random_clusters),
  clusters_path, sep = "\t", quote = FALSE, row.names = FALSE
)

# PCA file = no PCA is computed, so we write two random columns only so the metrics stage can run
write.table(
  data.frame(cell_id = cell_ids, PC1 = rnorm(n_cells), PC2 = rnorm(n_cells)),
  pca_path, sep = "\t", quote = FALSE, row.names = FALSE
)

# timings and clustering info (minimal, but required by the metrics stage)
write_json(list(louvain = elapsed, leiden = elapsed), timings_path, auto_unbox = TRUE, pretty = TRUE)
write_json(
  list(resolutions = list(louvain = NA, leiden = NA), num_runs = list(louvain = 1, leiden = 1)),
  clustering_info_path, auto_unbox = TRUE, pretty = TRUE
)
