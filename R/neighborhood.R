#' Title
#'
#' @param meta_data
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
neighborhood_list_function <- function(meta_data){
  i<- NULL
  uid_split <- split(meta_data, meta_data$uID)
  cores = parallel::detectCores()
  clust <- parallel::makeSOCKcluster(cores)
  doSNOW::registerDoSNOW(clust)

  nbd_mat <- foreach::foreach(i = uid_split, .combine = rbind,
                     .packages = c('dplyr', 'data.table', 'rdist')) %dopar%
    (data.table::data.table(nbd_maker(i)))
  parallel::stopCluster(clust)

  nbd_mat2 = lapply(1:nrow(nbd_mat), function(i)(nbd_mat$V1[i][[1]]))

  return(nbd_mat2)
}
