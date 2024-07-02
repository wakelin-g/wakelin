#' Quickly cluster \linkS4class{SingleCellExperiment} object
#'
#' @param sce A \linkS4class{SingleCellExperiment} object to be clustered.
#' @param method A character vector specifying either \code{walktrap} or \code{leiden} graph clustering algorithms.
#' @param label A character vector specifying the column in \code{colData} the cluster results should be stored in.
#' Defaults to \code{method}.
#'
#' @return Clustered \linkS4class{SingleCellExperiment} with clusters saved in \code{label} column of \code{colData}.
#' @export
#'
#' @examples
#' sce <- scater::mockSCE()
#' sce <- scater::runPCA(sce)
#' sce <- quick_cluster(sce, method = "walktrap", label = "cluster")
#' str(sce$cluster)
quick_cluster <- function(sce, method = c("walktrap", "leiden"), label = NULL) {
  method <- match.arg(method)
  label <- ifelse(is.null(label), method, label)
  dimred <- ifelse("corrected" %in% SingleCellExperiment::reducedDimNames(sce), "corrected", "PCA")
  snn <- scran::buildSNNGraph(sce, use.dimred = dimred)

  if (method == "walktrap") {
    sce@colData[[label]] <- factor(igraph::cluster_walktrap(snn)$membership)
  } else if (method == "leiden") {
    sce@colData[[label]] <- factor(igraph::cluster_leiden(snn)$membership)
  } else {
    stop("Clustering method not found. Please choose either `walktrap` or `leiden`.")
  }

  sce
}
