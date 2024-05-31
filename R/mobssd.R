#' MoBSSD function
#'
#' This function fits the latent states based model to the threshold data.
#'
#' @param threshold_data
#' @param meta_data
#' @param L
#' @param neighborhood_list
#' @param parameters_initial_list
#' @param niter
#' @param gamma_fixed
#' @param gamma
#' @param projection_type
#'
#' @return domains
#' @return omega
#' @return probability_matrix
#' @return meta_data_with_domains
#' @export
#'
#' @examples
mobssd <- function(threshold_data, meta_data, num_domains, neighborhood_list = NULL, parameters_initial_list = NULL, niter = 200, gamma_fixed = FALSE, gamma = 0, projection_type = NULL){

  L <- num_domains

  if (is.null(parameters_initial_list)) {
    if (is.null(projection_type)){
      stop("Provide a projection type (pca or umap) for parameters initialization or a list of initial values (parameters_initial_list)")
    }else{
    parameters_initial_list = parameters_init_function(threshold_data = threshold_data,meta_data = meta_data,projection_type = projection_type,num_domains = L)
    }
  }
  if (is.null(neighborhood_list)) {
    neighborhood_list = neighborhood_list_function(meta_data = meta_data)
  }

  p_mat_kmeans = parameters_initial_list$pmat_initial
  x_vec_kmeans = parameters_initial_list$domain_initial
  omega_vec_kmeans = parameters_initial_list$omega_initial

  result_list = potts_model(L = L,niter = niter, threshold_data = threshold_data,omega_vec = omega_vec_kmeans,
                            p_mat = p_mat_kmeans,x_vec = x_vec_kmeans,meta_data = meta_data,
                            nbd_mat = neighborhood_list,gamma_fixed = gamma_fixed,gamma=gamma)

  meta_data = meta_data %>% dplyr::mutate(domain = result_list$assigned_domains)

return(domains = result_list$assigned_domains,
       omega = result_list$omega,
       probability_matrix = result_list$pmat, meta_data_with_domains = meta_data)

}
