run_droplet_processing <- function(fpath, out_save = FALSE, out_dir = "scrnaseq/", plot_save = FALSE, force = FALSE, verbose = FALSE) {
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

  sce <- DropletUtils::read10xCounts(paste0(fpath, "/raw_feature_bc_matrix"), col.names = TRUE, row.names = "symbol")

  if (plot_save) {
    msg_verbose("Creating barcode rank plots.")
    if (!dir.exists(plotting_dir)) {
      dir.create(plotting_dir)
    }
    bc_ranks <- DropletUtils::barcodeRanks(SingleCellExperiment::counts(sce))
    bc_ranks_unique <- !duplicated(bc_ranks$rank)
    grDevices::svg(filename = paste0(plotting_dir, "/barcode_ranks.svg"))
    plot(bc_ranks$rank[bc_ranks_unique], bc_ranks$total[bc_ranks_unique], log = "xy", xlab = "Rank", ylab = "Total UMI Counts", cex.lab = 1.2, main = fname)
    abline(h = S4Vectors::metadata(bc_ranks)$inflection, col = "darkgreen", lty = 2)
    abline(h = S4Vectors::metadata(bc_ranks)$knee, col = "skyblue", lty = 2)
    legend("bottomleft", legend = c("Inflection", "Knee"), col = c("darkgreen", "skyblue"), lty = 2, cex = 1.2)
    grDevices::dev.off()
  }

  msg_verbose("Calculating empty droplets.")

  empty_out <- DropletUtils::emptyDrops(SingleCellExperiment::counts(sce), lower = 50, test.ambient = TRUE, BPPARAM = BiocParallel::MulticoreParam(8))
  sce <- sce[,which(empty_out$FDR <= 0.001)]

  if (out_save) {
    msg_verbose("Saving file.")
    readr::write_rds(sce, outfile)
  }
}
