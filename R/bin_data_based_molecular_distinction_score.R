gene_function_wilcox <- function(gene,meta_data, bin_data, num_domains){
  triplets = list()
  k=1
  for (state1 in seq(1,num_domains)) {
    for (state2 in seq(1,num_domains)) {

      if (state1 != state2) {

        gene_data = bin_data[,gene]
        state_indices1 = which(meta_data$domain == state1)
        state_indices2 = which(meta_data$domain == state2)

        state_gene_vec1 = gene_data[state_indices1]
        state_gene_vec2 = gene_data[state_indices2]

        test =  stats::wilcox.test(x=state_gene_vec1, y = state_gene_vec2,alternative = "greater")
        triplets[[k]] = list(gene_name = gene, state1 = state1, state2 = state2, p_val = test$p.value)

      }
      k=k+1
    }
  }
  return(triplets)
}

#' Distinction p values for each gene for each pair of domains
#'
#' @param bin_data
#' @param meta_data
#'
#' @return
#' @export
#'
#' @examples
bin_data_based_distinction <- function(bin_data, meta_data){

  domain <- i <- state1 <- state2 <- NULL

  num_domains = length(unique(meta_data$domain))
  selected_index = meta_data[,"index"]
  bin_data = bin_data[selected_index,]

  gene_names = colnames(bin_data)

  cores = parallel::detectCores()
  clust <- parallel::makeSOCKcluster(cores)
  doSNOW::registerDoSNOW(clust)
  result_triplet <- foreach::foreach(i = 1:length(gene_names),
                          .packages = c('dplyr', 'data.table')) %dopar%
    gene_function_wilcox(gene = gene_names[i], bin_data = bin_data, num_domains = num_domains,meta_data = meta_data)
  parallel::stopCluster(clust)

  result_triplet = data.table::data.table(result_triplet)

  score_mat = array(data = NA, dim = c(num_domains, num_domains))
  for (s1 in seq(1,num_domains)) {
    for (s2 in seq(1,num_domains)) {
      k=1
      if (state1 != state2) {
        p_vec = as.vector(result_triplet %>% dplyr::filter(state1==s1, state2==s2) %>% dplyr::select("p_val"))
        p_vec_adjusted = stat::p.adjust(p=p_vec, method = "BH")
        score_mat[state1,state2] = sum(p_vec_adjusted <= 0.05)
      }
      k=k+1
    }
  }

  average_score = mean(score_mat,na.rm=T)

  return(list(p_value_df = result_triplet, score_mat = score_mat, average_score = average_score))


}






