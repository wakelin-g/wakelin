quick_cluster <- function(sce, method = c("walktrap", "leiden")) {
  method <- match.arg(method)
  dimred <- ifelse("corrected" %in% SingleCellExperiment::reducedDimNames(sce), "corrected", "PCA")
  snn <- scran::buildSNNGraph(sce, use.dimred = dimred)
}
