run_soupx <- function(fpath, out_save = FALSE, out_dir = "scrnaseq/outs_soupx/", plot_save = FALSE, force = FALSE, verbose = TRUE) {
  fname <- basename(tools::file_path_sans_ext(fpath))
  outfile <- paste0(out_dir, fname, "_soupx.mtx")
  if (verbose) {
    message("--- Ambient RNA estimation ---")
  }
  if (file.exists(outfile)) {
    if (!force) {
      message("File: '", outfile, "' already exists. Skipping.")
      return(NULL)
    }
  }

  if (verbose) {
    message("Opening dataset '", fname, "'...")
  }
  soup <- SoupX::load10X(fpath, verbose = verbose)
  if (plot_save) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    if (verbose) {
      message("Estimating contamination; saving soupX plot.")
    }
    grDevices::svg(filename = paste0(out_dir, fname, "_soupx_est_cont.svg"))
    soup <- SoupX::autoEstCont(soup, doPlot = TRUE, verbose = verbose)
    grDevices::dev.off()
  } else {
    if (verbose) {
      message("Estimating contamination; not saving soupX plot.")
    }
    soup <- SoupX::autoEstCont(soup, doPlot = FALSE, verbose = verbose)
  }
  if (verbose) {
    message("Adjusting counts based on estimated contamination.")
  }
  soup_out <- SoupX::adjustCounts(soup, verbose = 1)
  if (out_save) {
    if (verbose) {
      message("Saving soupX outs to '", out_dir, fname, "_soupx.mtx'...")
      if (!dir.exists(out_dir)) {
        dir.create(out_dir)
      }
      Matrix::writeMM(soup_out, file = outfile)
    }
  }
}
