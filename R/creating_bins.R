#' Binning Each Cell
#'
#' The function that takes transcripts coordinates as inputs and divide each cell into square bins. User needs to specify the binsize.
#' @param data
#' @param binsize
#'
#' @return bin_data: The binned data where rows are bins and columns are all genes. Each entry (i,j) is the number of transcripts of Gene(j) in Bin(i).
#' @return meta_data: The meta data for each bin. This contains uID, centers of each bin (node_center_x, node_center_y), nucleus state (0: inside the nucleus, 1: outside the nucleus, -1: boundary of the nucleus ), and index column
#' @export
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @examples
bining_function <- function(data, binsize){
  i <- NULL
  gene_names_all = unique(data$gene)
  data$inNucleus = ifelse(data$inNucleus >=0, 0,1)
  all_genes = unique(data$genes)

  uid_split <- split(data, data$uID)

  cores = parallel::detectCores()
  clust <- parallel::makeSOCKcluster(cores)
  doSNOW::registerDoSNOW(clust)
  bin_data <- foreach::foreach(i = uid_split, .combine = rbind,
                      .packages = c('dplyr', 'data.table')) %dopar%
    (data.table::data.table(cell_bins(i, binsize = binsize, all_genes = all_genes)))
  parallel::stopCluster(clust)


  meta_data = bin_data %>%
    dplyr::select(c("uID", "node_center_x", "node_center_y","nucleus")) %>%
    dplyr::rename("nuc_state" = "nucleus")
  meta_data = meta_data %>% dplyr::mutate(index = seq(1,nrow(meta_data)))

  bin_data = bin_data %>%
    dplyr::select(!c("uID", "node_center_x", "node_center_y","nucleus"))

  return(list(bin_data = bin_data,
              meta_data = meta_data))
}

