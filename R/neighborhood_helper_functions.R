nbd_def <- function(dist_vec, nbd_size = nbd_size){
  return(which((dist_vec < nbd_size)&(dist_vec>0)))
}

nearest_def <- function(dist_vec, i){
  new_vec = dist_vec[-i]
  return(which.min(new_vec))
}


nbd_maker <- function(group){
  coordinates_temp = as.matrix(group %>% dplyr::select("node_center_x","node_center_y"))

  dist_mat = rdist::pdist(coordinates_temp)
  nbd_mat_temp = lapply(1:nrow(dist_mat), function(i)nbd_def(dist_vec = dist_mat[i,]))
  nbd_nearest = lapply(1:nrow(dist_mat), function(i)nearest_def(dist_vec = dist_mat[i,],i = i))

  problem_index = which(unlist(lapply(1:length(nbd_mat_temp), function(i)(is.na(nbd_mat_temp[[i]][1])))))

  for (i in problem_index) {
    nbd_mat_temp[[i]] = nbd_nearest[[i]]
  }

  nbd_mat_correct = lapply(1:length(nbd_mat_temp), function(i)group$index[nbd_mat_temp[[i]]])
  return(nbd_mat_correct)
}
