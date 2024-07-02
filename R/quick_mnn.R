#' MNN Integrate \linkS4class{SingleCellExperiment} Objects
#'
#' @param ... \linkS4class{SingleCellExperiment} classes to be MNN-integrated.
#'
#' @return \linkS4class{SingleCellExperiment} object with \code{corrected} object in \code{reducedDims}.
#' @export
#'
#' @examples
#' sce1 <- scater::mockSCE()
#' sce2 <- scater::mockSCE()
#' sce.combined <- quick_mnn(sce1, sce2)
quick_mnn <- function(..., cluster = FALSE) {
  combined <- batchelor::correctExperiments(...)
  combined
}
