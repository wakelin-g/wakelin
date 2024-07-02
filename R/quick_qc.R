#' Perform quick one-step QC
#'
#' @param sce A \linkS4class{SingleCellExperiment} for QC metrics to be calculated.
#' @param mito Logical scalar indicating if percentage of mitochondrial reads should be added to colData.
#' @param ribo Logical scalar indicating if percentage of ribosomal reads should be added to colData.
#' @param hemo Logical scalar indicating if percentage of hemoglobin reads should be added to colData.
#' @param doublets Logical scalar indicating if per-cell droplet calculations should be performed.
#' @param verbose Logical scalar indicating if outputs should be verbose.
#'
#' @return A \linkS4class{SingleCellExperiment} with QC metrics added to colData.
#' @export
#'
#' @examples
#' sce <- scuttle::mockSCE()
#' sce.qc <- quick_qc(sce, verbose = TRUE)
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom sparseMatrixStats colSums2
#' @importFrom scDblFinder computeDoubletDensity
#' @importFrom scuttle isOutlier
quick_qc <- function(sce, mito = TRUE, ribo = TRUE, hemo = TRUE, doublets = TRUE, verbose = TRUE) {
  stopifnot(class(sce) == "SingleCellExperiment")

  mat <- counts(sce)
  genes <- rownames(sce)

  sce@colData$UMIs <- colSums2(mat[, 1:ncol(sce)])
  sce@colData$unique_genes <- colSums2((\(x)(x != 0))(mat))

  if (mito) {
    msg_verbose("COMPUTING MITO GENES")
    mito_genes <- genes[grepl("^mt-", genes)]
    sce@colData$total_mt <- colSums2(mat[mito_genes, ])
    sce@colData$percent_mt <- sce@colData$total_mt / sce@colData$UMIs
  }
  if (ribo) {
    msg_verbose("COMPUTING RIBO GENES")
    ribo_genes <- genes[grepl("^Rp[sl]", genes)]
    sce@colData$total_rb <- colSums2(mat[ribo_genes, ])
    sce@colData$percent_rb <- sce@colData$total_rb / sce@colData$UMIs
  }
  if (hemo) {
    msg_verbose("COMPUTING HEMO GENES")
    hemo_genes <- genes[grepl("^Hb[^(p|e|s)]", genes)]
    sce@colData$total_hb <- colSums2(mat[hemo_genes, ])
    sce@colData$percent_hb <- sce@colData$total_hb / sce@colData$UMIs
  }
  if (doublets) {
    msg_verbose("COMPUTING DOUBLETS")
    sce@colData$cell_doublet_scores <- computeDoubletDensity(sce)
  }
  sce
}
