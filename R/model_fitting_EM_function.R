em_iteration_variable_gamma <- function(L,x_vec,y_mat,p_mat,nbd_mat, omega_vec,gamma){
  state_list <- seq(1,L)

  ncores = parallel::detectCores()

  map_list <- do.call(rbind,foreach::mclapply(1:length(x_vec), function(i)map_calc_i2(i = i,p_mat = p_mat,
                                                                             y_mat = y_mat,
                                                                             state_list = state_list,
                                                                             nbd_mat = nbd_mat,
                                                                             x_vec = x_vec,
                                                                             omega_vec = omega_vec, gamma=gamma),mc.cores = ncores))
  x_vec_new = as.numeric(unlist(map_list[,1]))
  log_l_mat = apply(do.call(rbind,map_list[,2]),2,as.numeric)
  interaction_mat1 = apply(do.call(rbind,map_list[,3]),2,as.numeric)

  posterior_coeffs <- foreach::mclapply(1:length(x_vec), function(i)posterior_coeffs_i2(i = i,x_vec = x_vec_new,
                                                                               nbd_mat = nbd_mat,
                                                                               state_list = state_list,
                                                                               log_l_mat = log_l_mat,
                                                                               omega_vec = omega_vec, gamma=gamma),mc.cores = ncores)
  posterior_coeff_mat = do.call(rbind, posterior_coeffs)
  posterior_mat <- do.call(rbind, posterior_coeff_mat[,1])
  interaction_mat <- do.call(rbind, posterior_coeff_mat[,2])

  posterior_mat_p <- sweep(t(posterior_mat),1,FUN = "/",STATS = colSums(posterior_mat))

  p_mat_new <- (posterior_mat_p%*%y_mat)

  confun <- function(beta){
    len_b = length(beta)
    omega = beta[1:(len_b-1)]
    alpha = exp(-omega)
    f = sum(alpha) -1
    return(f)
  }

  objfun <- function(beta){
    len_b = length(beta)
    omega = beta[1:(len_b-1)]
    gamma = beta[len_b]
    alpha=  exp(-omega)
    val1 = sum(posterior_mat%*%omega) + gamma*sum(posterior_mat*interaction_mat)
    val2 = sum(log(exp(-gamma*interaction_mat)%*%alpha))
    total  =val1+val2
    return(total)
  }
  beta = c(omega_vec, gamma)

  non_linear_optim = Rsolnp::solnp(pars = beta, fun = objfun, eqfun = confun, LB = c(rep(0,L),-10))
  beta = non_linear_optim$pars

  omega_vec_new = beta[1:L]
  gamma_new = beta[L+1]

  return(list(x_vec = x_vec_new, p_mat = p_mat_new, omega_vec = omega_vec_new, gamma = gamma_new))
}



em_iteration_fixed_gamma <- function(L,x_vec,y_mat,p_mat,nbd_mat, omega_vec,gamma){
  state_list <- seq(1,L)

  ncores = parallel::detectCores()

  map_list <- do.call(rbind,foreach::mclapply(1:length(x_vec), function(i)map_calc_i2(i = i,p_mat = p_mat,
                                                                             y_mat = y_mat,
                                                                             state_list = state_list,
                                                                             nbd_mat = nbd_mat,
                                                                             x_vec = x_vec,
                                                                             omega_vec = omega_vec, gamma=gamma),mc.cores = ncores))
  x_vec_new = as.numeric(unlist(map_list[,1]))
  log_l_mat = apply(do.call(rbind,map_list[,2]),2,as.numeric)
  interaction_mat1 = apply(do.call(rbind,map_list[,3]),2,as.numeric)


  posterior_coeffs <- foreach::mclapply(1:length(x_vec), function(i)posterior_coeffs_i2(i = i,x_vec = x_vec_new,
                                                                               nbd_mat = nbd_mat,
                                                                               state_list = state_list,
                                                                               log_l_mat = log_l_mat,
                                                                               omega_vec = omega_vec, gamma=gamma),mc.cores = ncores)
  posterior_coeff_mat = do.call(rbind, posterior_coeffs)
  posterior_mat <- do.call(rbind, posterior_coeff_mat[,1])
  interaction_mat <- do.call(rbind, posterior_coeff_mat[,2])

  posterior_mat_p <- sweep(t(posterior_mat),1,FUN = "/",STATS = colSums(posterior_mat))

  p_mat_new <- (posterior_mat_p%*%y_mat)

  confun <- function(beta){
    alpha = exp(-beta)
    f = sum(alpha) -1
    return(f)
  }


  objfun <- function(beta){
    alpha=  exp(-beta)
    val1 = sum(posterior_mat%*%beta) + gamma*sum(posterior_mat*interaction_mat)
    val2 = sum(log(exp(-gamma*interaction_mat)%*%alpha))
    total = val1+val2
    return(total)
  }
  beta = omega_vec

  non_linear_optim = Rsolnp::solnp(pars = beta, fun = objfun, eqfun = confun, LB = c(rep(0,L)))
  omega_vec_new = non_linear_optim$pars

  return(list(x_vec = x_vec_new, p_mat = p_mat_new, omega_vec = omega_vec_new))
}
