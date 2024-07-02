#' Quickly filter cells based on QC metrics
#'
#' @param sce A \linkS4class{SingleCellExperiment} to be filtered.
#' @param features An optional list of features in colData to be filtered on.
#' @param discard Logical scalar indicating whether cells not passing the filters should be discarded.
#' @param nmads Numeric scalar indicating the of median absolute deviations to filter on.
#' @param min_counts Numeric scalar indicating the lower threshold for UMI number.
#' @param max_counts Numeric scalar indicating the upper threshold for UMI number.
#'
#' @return A filtered \linkS4class{SingleCellExperiment} (if \code{discard = TRUE}) or a \linkS4class{SingleCellExperiment}
#' with logical \code{discard} column added to its \code{colData}, indicating which cells did not pass the filters.
#'
#' @export
#'
#' @examples
#' sce <- scater::mockSCE()
#' sce.qc <- quick_qc(sce)
#' sce.filt <- quick_filter(sce.qc)
quick_filter <- function(sce, features = NULL, discard = TRUE, nmads = 3, min_counts = NULL, max_counts = NULL) {
  qc_list <- list()
  qc_metrics <- c("percent_mt", "percent_rb", "percent_hb", "cell_doublet_scores")

  if (is.null(features)) {
    if (is.null(sce@colData$UMIs)) {
      stop("Did not find `UMIs` in colData. Did you run `quick_qc()` beforehand?", call. = FALSE)
    }
    for (metric in qc_metrics) {
      if (!is.null(sce@colData[[metric]])) {
        qc_list[[metric]] <- as.logical(isOutlier(sce@colData[[metric]], nmads = nmads))
      }
    }
  } else {
    qc_metrics <- c(qc_metrics, features)
    for (metric in qc_metrics) {
      if (!is.null(sce@colData[[metric]])) {
        qc_list[[metric]] <- as.logical(isOutlier(sce@colData[[metric]], nmads = nmads))
      }
    }
  }

  if (!is.null(min_counts)) {
    qc_list[["min_counts"]] <- sce@colData$UMIs < min_counts
  }

  if (!is.null(max_counts)) {
    qc_list[["max_counts"]] <- sce@colData$UMIs > max_counts
  }

  qc_list <- do.call(cbind.data.frame, qc_list)
  cells_discard <- apply(qc_list, 1, any)
  if (!discard) {
    sce@colData$discard <- cells_discard
    return(sce)
  } else {
    sce <- sce[, !cells_discard]
  }
}
