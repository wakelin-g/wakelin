#' Run SoupX decontamination
#'
#' `run_soupx()` wraps the normal SoupX workflow, with default parameters, into a single
#' function for convenience, with options as to whether the plots and outputs should be saved to disk.
#'
#' @param fpath Path to cellranger `outs` directory, which needs to contain `raw_feature_bc_matrix`, `filtered_feature_bc_matrix`, and `analysis` subdirectories.
#' @param out_save Save SoupX output matrix?
#' @param out_dir Path to output directory where results (and possibly plots) will be saved.
#' @param plot_save Save SoupX plot?
#' @param force If SoupX output matrix already found, should it be overwritten?
#' @param verbose Be verbose?
#'
#' @return SoupX-corrected gene expression matrix.
#' @export
#'
#' @examples
#' fpath <- "/Users/griffen/Documents/thesis-code/scrnaseq/outs/tw2_uninjured/"
#' soup <- wakelin::run_soupx(fpath = fpath, out_save = FALSE, plot_save = FALSE)
#'
#' @importFrom Matrix writeMM
#' @importFrom SoupX load10X autoEstCont adjustCounts
#' @importFrom grDevices svg dev.off
run_soupx <- function(fpath, out_save = FALSE, out_dir = "scrnaseq/outs_soupx/", plot_save = FALSE, force = FALSE, verbose = TRUE) {
  fname <- get_fname(fpath)
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
  soup <- load10X(fpath, verbose = verbose)
  if (plot_save) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    if (verbose) {
      message("Estimating contamination; saving soupX plot.")
    }
    svg(filename = paste0(out_dir, fname, "_soupx_est_cont.svg"))
    soup <- autoEstCont(soup, doPlot = TRUE, verbose = verbose)
    dev.off()
  } else {
    if (verbose) {
      message("Estimating contamination; not saving soupX plot.")
    }
    soup <- autoEstCont(soup, doPlot = FALSE, verbose = verbose)
  }
  if (verbose) {
    message("Adjusting counts based on estimated contamination.")
  }
  soup_out <- adjustCounts(soup, verbose = 1)
  if (out_save) {
    if (verbose) {
      message("Saving soupX outs to '", out_dir, fname, "_soupx.mtx'...")
      if (!dir.exists(out_dir)) {
        dir.create(out_dir)
      }
      writeMM(soup_out, file = outfile)
    }
  }
  soup_out
}
