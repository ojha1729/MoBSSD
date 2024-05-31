log_l_i <- function(p_mat, y_vec,state){
  p_vec <- p_mat[state,]
  p_vec = unlist(lapply(1:length(p_vec), function(k)min(0.9999, p_vec[k])))
  p_vec = unlist(lapply(1:length(p_vec), function(k)max(0.0001, p_vec[k])))
  value = sum(y_vec*log(p_vec) + (1-y_vec)*log(1-p_vec))
  return(value)
}

vc_function2 <- function(i,state_l,x_vec,nbd_states,omega_l,gamma){
  potential_values <- sum(state_l == nbd_states)
  return(list(total=(gamma*potential_values + omega_l), interaction = potential_values))
}


map_calc_i2 <- function(i,p_mat,y_mat,state_list,nbd_mat,x_vec, omega_vec,gamma){
  y_vec = y_mat[i,]
  nbd_i = nbd_mat[[i]]
  nbd_states <- x_vec[nbd_i]
  l_list <- lapply(1:length(state_list), function(j)log_l_i(p_mat = p_mat,state = state_list[j],
                                                            y_vec = y_vec))
  l_vec = unlist(l_list)
  vc_list = do.call(rbind,lapply(1:length(state_list), function(j)vc_function2(i= i,state_l = state_list[j],
                                                                               x_vec = x_vec,nbd_states = nbd_states,
                                                                               omega_l = omega_vec[j],gamma = gamma)))
  vc_vec  = as.numeric(unlist(vc_list[,1]))
  interaction_vec1 = as.numeric(unlist(vc_list[,2]))
  map_val_vec = -l_vec + vc_vec
  map_val = which.min(map_val_vec)

  return(list(map = map_val,log_l_vec = l_vec, interaction_vec1 = interaction_vec1))
}

posterior_coeffs_i2 <- function(i,x_vec,nbd_mat,state_list,log_l_mat, omega_vec,gamma){
  nbd_i = nbd_mat[[i]]
  log_l_vec = exp(log_l_mat[i,])
  nbd_states <- x_vec[nbd_i]
  vc_list = do.call(rbind,lapply(1:length(state_list), function(j)vc_function2(i= i,state_l = state_list[j],
                                                                               x_vec = x_vec,nbd_states = nbd_states,
                                                                               omega_l = omega_vec[j], gamma=gamma)))
  vc_vec1 = exp(-as.numeric(unlist(vc_list[,1])))
  interaction_vec = as.numeric(unlist(vc_list[,2]))
  vc_vec = vc_vec1/sum(vc_vec1)


  if (min(log_l_mat[i,]) > 740) {
    min_index = which.min(log_l_mat[i,])
    coeff_vec = rep(0,length(omega_vec))
    coeff_vec[min_index] = 1
  }else{
    l_vec = exp(log_l_mat[i,])
    coeff_vec1 <- vc_vec*log_l_vec
    coeff_vec <- coeff_vec1/sum(coeff_vec1)}

  return(list(coeff_vec = coeff_vec, interaction_vec = interaction_vec ))
}
