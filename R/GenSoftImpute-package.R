## usethis namespace: start
#' @useDynLib GenSoftImpute
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("GenSoftImpute", libpath)
}
