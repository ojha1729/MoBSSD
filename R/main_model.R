#' Model for Subcellular Spatial Domain Discovery
#'
#' @param L
#' @param niter
#' @param threshold_data
#' @param omega_vec
#' @param p_mat
#' @param x_vec
#' @param meta_data
#' @param nbd_mat
#' @param gamma_fixed
#'
#' @return assigned_domains
#' @return omega
#' @return gamma
#' @return pmat
#'
#' @export
#'
#' @examples
#' @importFrom magrittr %>%
potts_model <- function(L,niter=200,threshold_data,omega_vec,p_mat,x_vec, meta_data, nbd_mat, gamma_fixed = FALSE, gamma = 0){
  uID <- NULL
  y_mat = threshold_data

  G = ncol(y_mat)

  diff = 100
  k=2

  ncores_available = parallel::detectCores()


  while((diff> 0.0001)&(k<niter)) {

    if (ncores_available < 24) {
      ncells = length(unique(meta_data$uID))

      cell_samples = sample((unique(meta_data$uID)),ceiling(ncells/5))
      sampled_index = meta_data %>% dplyr::filter(uID %in% cell_samples)
      sampled_index = as.vector(sampled_index$sampling_index)

      x_vec_sample = x_vec[sampled_index]
      y_mat_sample = y_mat[sampled_index,]

      if (gamma_fixed) {
        iteration = em_iteration_fixed_gamma(L = L,omega_vec = omega_vec,
                                 x_vec = x_vec_sample,p_mat = p_mat,y_mat = y_mat_sample,
                                 nbd_mat = nbd_mat,gamma = gamma)

        gamma_new = gamma
      }
      else{
        iteration = em_iteration_variable_gamma(L = L,omega_vec = omega_vec,
                                 x_vec = x_vec_sample,p_mat = p_mat,y_mat = y_mat_sample,
                                 nbd_mat = nbd_mat,gamma = gamma)

        gamma_new = iteration$gamma
      }


      omega_vec_new = iteration$omega_vec
      p_mat_new = iteration$p_mat

      x_vec_sample = iteration$x_vec
      x_vec[sampled_index] = x_vec_sample

      diff = max(max(abs((omega_vec - omega_vec_new)/omega_vec)),
                 max(abs((p_mat_new - p_mat)/p_mat)),abs((gamma - gamma_new)/gamma))

      omega_vec = omega_vec_new
      p_mat = p_mat_new
      gamma = gamma_new

    }

    else{
      if (gamma_fixed) {
        iteration = em_iteration_fixed_gamma(L = L,omega_vec = omega_vec,
                                             x_vec = x_vec,p_mat = p_mat,y_mat = y_mat,
                                             nbd_mat = nbd_mat,gamma = gamma)

        gamma_new = gamma
      }
      else{
        iteration = em_iteration_variable_gamma(L = L,omega_vec = omega_vec,
                                                x_vec = x_vec,p_mat = p_mat,y_mat = y_mat,
                                                nbd_mat = nbd_mat,gamma = gamma)

        gamma_new = iteration$gamma
      }


      omega_vec_new = iteration$omega_vec
      p_mat_new = iteration$p_mat
      x_vec = iteration$x_vec

      diff = max(max(abs((omega_vec - omega_vec_new)/omega_vec)),
                 max(abs((p_mat_new - p_mat)/p_mat)),abs((gamma - gamma_new)/gamma))

      omega_vec = omega_vec_new
      p_mat = p_mat_new
      gamma = gamma_new



    }

    k=k+1
  }


  return(list(assigned_domains = x_vec, omega = omega_vec,gamma = gamma, pmat = p_mat))

}
