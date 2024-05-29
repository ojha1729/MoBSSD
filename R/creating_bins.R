
#' Title
#'
#' @param data
#'
#' @return
#' @export
#' @importFrom foreach %dopar%
#' @examples
bining_function <- function(data){
  gene_names_all = unique(data$gene)
  data$inNucleus = ifelse(data$inNucleus >=0, 0,1)

  uid_split <- split(data, data$uID)

  cores = parallel::detectCores()
  clust <- parallel::makeSOCKcluster(cores)
  doSNOW::registerDoSNOW(clust)
  bin_data <- foreach::foreach(i = uid_split, .combine = rbind,
                      .packages = c('dplyr', 'data.table')) %dopar%
    (data.table::data.table(cell_bins(i)))
  parallel::stopCluster(clust)

}

