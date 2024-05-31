contiguity_score <- function(assigned_domains = NULL, meta_data, neighborhood_list){

  domain <- uID <- state <- score_vec <- median_contiguity <- NULL

  if (is.null(assigned_domains)) {
    state_data = data.table::data.table(state = meta_data$domain)
  }else{
    state_data = data.table::data.table(state = assigned_domains)
  }

  state_data = state_data %>% dplyr::mutate(uID = meta_data$uID)
  state_data = state_data %>% dplyr::mutate(index= seq(1,nrow(state_data)))

  state_data = state_data %>% dplyr::mutate(score_vec =  unlist(lapply(1:nrow(state_data),
                                                                    function(i)mean(state_data$state[i] == state_data$state[neighborhood_list[[(state_data$index)[i]]]]))))

  temp = state_data %>% dplyr::group_by(uID, state) %>% dplyr::summarise(median_contiguity = stats::median(score_vec))

  temp2 = temp %>% dplyr::group_by(state) %>% dplyr::summarise(mean = mean(median_contiguity,na.rm = T))

  return(temp2$mean)

}





