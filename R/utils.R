get_fname <- function(fpath) {
  basename(tools::file_path_sans_ext(fpath))
}

msg_verbose <- function(...) {
  if (parent.frame()$verbose == TRUE) {
    message(...)
  }
}
