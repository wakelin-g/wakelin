basepath <- "inst/extdata/pbmc5k"
dir.create(basepath)

for (f in c("filtered_feature_bc_matrix.tar.gz", "raw_feature_bc_matrix.tar.gz", "analysis.tar.gz")) {
  url_current <- paste0("https://cf.10xgenomics.com/samples/cell-exp/7.0.1/SC3pv3_GEX_Human_PBMC/SC3pv3_GEX_Human_PBMC_", f)
  path_current <- file.path(basepath, f)
  download.file(url_current, path_current)
  untar(path_current, exdir = basepath)
  unlink(path_current)
}
