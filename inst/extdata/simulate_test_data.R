library(DropletUtils)

sim_empty <- function(empty_prof, n_empty, empty_rate, n_genes) {
  empty_prof <- empty_prof / sum(empty_prof)
  total_count <- rexp(n_empty, rate = empty_rate)
  empty_counts <- matrix(rpois(n_genes * n_empty, lambda = outer(empty_prof, total_count)), ncol = n_empty, nrow = n_genes)
  empty_counts <- as(empty_counts, "CsparseMatrix")
  rownames(empty_counts) <- paste0("Gene_", seq_len(n_genes))
  colnames(empty_counts) <- paste0("Empty_", seq_len(n_empty))
  return(empty_counts)
}

sim_small <- function(small_prof, n_small, small_rate, small_shape, n_genes) {
  small_prof <- small_prof / sum(small_prof)
  total_count <- rgamma(n_small, shape = small_shape, rate = small_rate)
  small_counts <- matrix(rpois(n_genes * n_small, lambda = outer(small_prof, total_count)), ncol = n_small, nrow = n_genes)
  small_counts <- as(small_counts, "CsparseMatrix")
  rownames(small_counts) <- paste0("Gene_", seq_len(n_genes))
  colnames(small_counts) <- paste0("Small_", seq_len(n_small))
  return(small_counts)
}

sim_large <- function(large_prof = large_prof, n_large = n_large, large_rate = large_rate, large_shape = large_shape, n_genes = n_genes) {
  large_prof <- large_prof / sum(large_prof)
  total_count <- rgamma(n_large, shape = large_shape, rate = large_rate)
  large_counts <- matrix(rpois(n_genes * n_large, lambda = outer(large_prof, total_count)), ncol = n_large, nrow = n_genes)
  large_counts <- as(large_counts, "CsparseMatrix")
  rownames(large_counts) <- paste0("Gene_", seq_len(n_genes))
  colnames(large_counts) <- paste0("Large_", seq_len(n_large))
  return(large_counts)
}

sim_counts <- function(
    n_genes = 100,
    n_empty = 10000, empty_prof = seq_len(n_genes), empty_rate = 0.04,
    n_small = 100, small_prof = runif(n_genes), small_shape = 20, small_rate = 0.1,
    n_large = 1000, large_prof = empty_prof, large_shape = 10, large_rate = 0.01,
    out_dir = NULL
)
# Used to simulate counts for testing. Shamelessly
# 'adapted' from Aaron Lun's simCounts in DropletUtils.
{
  empty <- sim_empty(empty_prof = empty_prof, n_empty = n_empty, empty_rate = empty_rate, n_genes = n_genes)
  small <- sim_small(small_prof = small_prof, n_small = n_small, small_rate = small_rate, small_shape = small_shape, n_genes = n_genes)
  large <- sim_large(large_prof = large_prof, n_large = n_large, large_rate = large_rate, large_shape = large_shape, n_genes = n_genes)

  filtered <- cbind(small, large)

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    write10xCounts(path = file.path(out_dir, "filtered_feature_bc_matrix"), x = filtered)
    raw <- cbind(filtered, empty)
    write10xCounts(path = file.path(out_dir, "raw_feature_bc_matrix"), x = raw)
    return(filtered)
  } else {
    return(filtered)
  }
}

sim_10x_counts <- function(out_dir) {
  mat <- sim_counts(out_dir = out_dir)
  dir.create(file.path(out_dir, "analysis", "clustering", "gene_expression_graphclust"), recursive = TRUE)
  snn <- scran::buildSNNGraph(mat)
  leiden <- igraph::cluster_leiden(snn)$membership
  write.csv(data.frame(list(Cluster=leiden,Barcode=colnames(mat))), file = file.path(out_dir, "analysis", "clustering", "gene_expression_graphclust", "clusters.csv"), quote = FALSE, row.names = FALSE)
  for (k in 2:10) {
    dir.create(file.path(out_dir, "analysis", "clustering", paste0("gene_expression_kmeans_", k, "_clusters")), recursive = TRUE)
    km <- unname(kmeans(Matrix::t(mat), centers = k)$cluster)
    write.csv(data.frame(list(Cluster=km,Barcode=colnames(mat))), file = file.path(out_dir, "analysis", "clustering", paste0("gene_expression_kmeans_", k, "_clusters"), "clusters.csv"))
  }
}

sim_10x_counts(out_dir = "inst/extdata/simdata")
