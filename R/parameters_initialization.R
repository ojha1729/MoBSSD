#' Parameters Initialization Function for EM algorithm
#'
#' This function uses UMAP or PCA based clustering approach to provide intilization values for the model parameters.
#'
#' @param threshold_data
#' @param meta_data
#' @param projection_type
#' @param num_domains
#'
#' @return domain_initial
#' @return pmat_initial
#' @return omega_initial
#'
#' @export
#'
#' @examples
parameters_init_function <- function(threshold_data, meta_data, projection_type,num_domains){
L <- num_domains
cluster <- x <- NULL
if (projection_type %in% c("u","umap","UMAP","Umap")) {

  data_sparse = Matrix::Matrix(data = threshold_data, sparse = T)
  umap_out = umap::umap(data_sparse)

  umap_df = umap_out$layout
  kmeans = stats::kmeans(x = umap_df,centers = L)

  data_clustered = as.data.frame(apply(threshold_data,2,as.numeric)) %>% dplyr::mutate(cluster = kmeans$cluster)

  p_mat_kmeans = data_clustered %>% dplyr::group_by(cluster) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean)) %>%
    dplyr::select(!cluster)

  p_mat_kmeans = as.matrix(p_mat_kmeans)

  x_vec_kmeans = data_clustered$cluster

  omega_vec_kmeans = as.vector(table(x_vec_kmeans)/length(x_vec_kmeans))

}

if (projection_type %in% c("pca","pc","PC","PCA","p","Pca")){

  data_sparse = Matrix::Matrix(data = threshold_data,sparse = T)
  pca = stats::prcomp(x = threshold_data,center = F,scale. = F,rank. = 10)
  pc_df = pca$x

  kmeans = stats::kmeans(x=x, centers = L)
  data_clustered = as.data.frame(apply(threshold_data,2,as.numeric)) %>% dplyr::mutate(cluster = kmeans$cluster)

  p_mat_kmeans = data_clustered %>% dplyr::group_by(cluster) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean)) %>%
    dplyr::select(!cluster)

  p_mat_kmeans = as.matrix(p_mat_kmeans)

  x_vec_kmeans = data_clustered$cluster

  omega_vec_kmeans = as.vector(table(x_vec_kmeans)/length(x_vec_kmeans))
}

else{
  print("Provide appropriate value for projection_type (umap or pca)")
}

return(list(pmat_initial = p_mat_kmeans, domain_initial = x_vec_kmeans, omega_initial = omega_vec_kmeans))

}
