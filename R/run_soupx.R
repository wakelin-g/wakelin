#' Run SoupX decontamination
#'
#' Wraps the normal SoupX workflow, with default parameters, into a single
#' function for convenience, with options as to whether the plots and outputs should be saved to disk.
#'
#' @param fpath A character vector indicating cellranger \code{outs} directory, which itself which must contain \code{raw_feature_bc_matrix}, \code{filtered_feature_bc_matrix}, and \code{analysis} subdirectories.
#' @param out_save Logical scalar indicating whether SoupX-corrected matrix should be written to disk.
#' @param out_dir A character vector indicating the directory to save outputs in.
#' @param plot_save Logical scalar indicating whether SoupX plot should be saved to disk.
#' @param force Logical scalar indicating that, if the SoupX output is already found, it should be overwritten.
#' @param verbose Logical scalar indicating whether function should intermittently print status.
#'
#' @return \linkS4class{sparseMatrix} corresponding to SoupX-corrected counts.
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

  msg_verbose("--- Ambient RNA estimation ---")

  if (file.exists(outfile)) {
    if (!force) {
      msg_verbose("File: '", outfile, "' already exists. Skipping.")
      return(NULL)
    }
  }

  msg_verbose("Opening dataset '", fname, "'...")
  soup <- load10X(fpath, verbose = verbose)
  if (plot_save) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    msg_verbose("Estimating contamination; saving soupX plot.")
    svg(filename = paste0(out_dir, fname, "_soupx_est_cont.svg"))
    soup <- autoEstCont(soup, doPlot = TRUE, verbose = verbose)
    dev.off()
  } else {
    msg_verbose("Estimating contamination; not saving soupX plot.")
    soup <- autoEstCont(soup, doPlot = FALSE, verbose = verbose)
  }
  msg_verbose("Adjusting counts based on estimated contamination.")
  soup_out <- adjustCounts(soup, verbose = 1)
  if (out_save) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    msg_verbose("Saving SoupX outs to'", out_dir, fname, "_soupx.mtx'.")
    writeMM(soup_out, file = outfile)
  }
  soup_out
}
