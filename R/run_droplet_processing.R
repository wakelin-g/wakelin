#' Load and Perform Droplet Processing on Data From 10x Genomics Dataset
#'
#' Creates a \linkS4class{SingleCellExperiment} from cellranger output directories, and performs simple barcode calling using \code{\link{emptyDrops}}.
#'
#' @param fpath A character vector indicating the directory containing the \code{filtered_feature_bc_matrix/} subdirectory.
#' @param out_save Logical scalar indicating whether output \linkS4class{SingleCellExperiment} should be written to disk.
#' @param out_dir A character vector indicating the directory to save outputs in.
#' @param plot_save Logical scalar indicating whether barcode rank plot should be saved to disk.
#' @param force Logical scalar indicating that, if the \linkS4class{SingleCellExperiment} output is already found, it should be overwritten.
#' @param verbose Logical scalar indicating whether function should intermittently print status.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object passed to relevant Bioconductor functions.
#'
#' @return A \linkS4class{SingleCellExperiment} object containing count data for each gene (row) and cell (column), with empty droplets removed.
#' @export
#'
#' @examples
#' run_droplet_processing(fpath = "/Users/griffen/Documents/thesis-code/scrnaseq/outs/tw2_uninjured/")
#' @importFrom BiocParallel MulticoreParam SerialParam
#' @importFrom graphics plot abline legend
#' @importFrom grDevices svg dev.off
#' @importFrom DropletUtils read10xCounts emptyDrops barcodeRanks
#' @importFrom SingleCellExperiment counts
#' @importFrom S4Vectors metadata
#' @importFrom readr write_rds
run_droplet_processing <- function(fpath, out_save = FALSE, out_dir = "scrnaseq/", plot_save = FALSE, force = FALSE, verbose = FALSE, BPPARAM = SerialParam()) {
  fname <- get_fname(fpath)
  outfile <- paste0(out_dir, "sces/", fname, "_droplet.rds")
  plotting_dir <- paste0(out_dir, "plots/", fname)

  msg_verbose("--- Droplet processing ---")

  if (file.exists(outfile)) {
    if (!force) {
      msg_verbose("File: '", outfile, "' already exists. Skipping.")
      return(NULL)
    }
  }

  msg_verbose("Opening dataset: '", fname, "'...")

  sce <- read10xCounts(paste0(fpath, "/raw_feature_bc_matrix"), col.names = TRUE, row.names = "symbol", BPPARAM = BPPARAM)

  if (plot_save) {
    msg_verbose("Creating barcode rank plots.")
    if (!dir.exists(plotting_dir)) {
      dir.create(plotting_dir)
    }
    bc_ranks <- barcodeRanks(counts(sce), BPPARAM = BPPARAM)
    bc_ranks_unique <- !duplicated(bc_ranks$rank)
    svg(filename = paste0(plotting_dir, "/barcode_ranks.svg"))
    plot(bc_ranks$rank[bc_ranks_unique], bc_ranks$total[bc_ranks_unique], log = "xy", xlab = "Rank", ylab = "Total UMI Counts", cex.lab = 1.2, main = fname)
    abline(h = metadata(bc_ranks)$inflection, col = "darkgreen", lty = 2)
    abline(h = metadata(bc_ranks)$knee, col = "skyblue", lty = 2)
    legend("bottomleft", legend = c("Inflection", "Knee"), col = c("darkgreen", "skyblue"), lty = 2, cex = 1.2)
    dev.off()
  }

  msg_verbose("Calculating empty droplets.")

  empty_out <- emptyDrops(counts(sce), lower = 50, test.ambient = TRUE, BPPARAM = BPPARAM)
  sce <- sce[,which(empty_out$FDR <= 0.001)]

  if (out_save) {
    msg_verbose("Saving file.")
    write_rds(sce, outfile)
  }
}
