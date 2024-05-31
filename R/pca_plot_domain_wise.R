mean_calc <- function(group){
  state <- NULL
  group = group %>% dplyr::select(!state)
  return(t(colMeans(group)))
}

#' Title
#'
#' @param uid
#' @param data
#' @param meta_data
#'
#' @return
#' @export
#'
#' @examples
mean_dat_generator <- function(uid,data, meta_data){
  domain <- uID <- . <- NULL
  x_vec = meta_data$domain

  cell = meta_data %>% dplyr::filter(uID == uid)
  uid_index = which(meta_data$uID ==uid)
  x_vec_cell = x_vec[uid_index]

  cell_data = data[uid_index,] %>% dplyr::mutate("domain" = factor(x_vec_cell))

  mean_vectors <- data.table::data.table(cell_data %>%
                                           dplyr::group_by(domain) %>% dplyr::do(data.table::data.table(mean_calc(.))) %>%
                                           dplyr::ungroup())
  mean_vectors = mean_vectors %>% dplyr::mutate("uID" = uid)

  return(mean_vectors)
}



#' PCA of mean expression of each domain from all cells
#'
#' @param meta_data
#' @param data
#' @param color_values
#'
#' @return
#' @export
#'
#' @examples
mean_expression_cell_wise_pca <- function(meta_data, data, color_values = NULL, mean_data = NULL){

  i <- domain <- PC1 <- PC2 <- NULL

  if (is.null(mean_data)) {
    uid_list <- unique(meta_data$uID)
    cores = parallel::detectCores()
    clust <- parallel::makeSOCKcluster(cores)
    doSNOW::registerDoSNOW(clust)
    mean_data <- foreach::foreach(i = uid_list, .combine = rbind,
                                  .packages = c('dplyr', 'data.table')) %dopar%
      (data.table::data.table(mean_dat_generator(uid = i,data = as.data.frame(data), meta_data = meta_data)))

    parallel::stopCluster(clust)
  }

  pc_data = scale(mean_data %>% dplyr::select(!c("domain","uID")))
  pc <- stats::prcomp(pc_data)
  pc_df = pc$x
  pc_df = data.frame("PC1" = pc_df[,1], "PC2" = pc_df[,2]) %>% dplyr::mutate("domain" = factor(mean_data$domain))
  if (is.null(color_values)) {
    p1 = ggplot2::ggplot(data = pc_df,ggplot2::aes(x=PC1,y=PC2,color = domain)) + ggplot2::geom_point(size=1)  +
      ggplot2::theme(axis.text = ggplot2::element_blank(),axis.ticks = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = 'white', color = 'purple'), legend.position = "bottom")+
      ggplot2::ggtitle("PCA of Mean Expression of different domains in each cell") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 10)))
  }else{
    p1 = ggplot2::ggplot(data = pc_df,ggplot2::aes(x=PC1,y=PC2,color = domain)) + ggplot2::geom_point(size=1)  +
      ggplot2::scale_color_manual(values = color_values)+
      ggplot2::theme(axis.text = ggplot2::element_blank(),axis.ticks = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = 'white', color = 'purple'), legend.position = "bottom")+
      ggplot2::ggtitle("PCA of Mean Expression of different domains in each cell") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 10)))
  }

  print(p1)

}

