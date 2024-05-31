#' Title
#'
#' @param x
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
threshold_non_zero <- function(x,threshold){
  q_val = stats::quantile(x[which(x>0)],probs = threshold)
  return(ifelse(test = (x >= q_val),yes = 1,no = 0))
}

#' Threshold Binned Data
#'
#' This function is used to threshold binned data.
#' This is an approach to reduce the impact of highly-expressed genes on the model outcomes.
#' @param bin_data
#' @param threshold
#' @param meta_data
#'
#' @return thresholded_data: A matrix of threshold data. Each gene is a binary vector where 1 represents those bins that have counts above the chosen threshold of that gene.
#' @return meta_data: Updated meta data
#'
#' @export
#'
#' @examples+
threshold_data_function <- function(bin_data, threshold, meta_data){
  temp = apply(bin_data,2,as.numeric)
  temp = do.call(cbind,lapply(1:ncol(temp), function(i)threshold_non_zero(x = temp[,i],threshold = threshold)))
  temp = as.matrix(apply(temp,2,as.numeric))
  colnames(temp) = colnames(bin_data)

  zero_index = which(rowSums(temp)== 0)
  temp = as.matrix(temp[-zero_index,])
  meta_data = meta_data[-zero_index,]

    return(list(thresholded_data = temp,
                meta_data = meta_data))
}
