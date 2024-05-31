
gene_function <- function(gene,df,num_domains){
  domain <- NULL
  triplets = list()
  k=1
  for (state1 in seq(1,num_domains)) {
    for (state2 in seq(1,num_domains)) {

      if (state1 != state2) {

        bins_state1 = df %>% dplyr::filter(domain == state1)
        success1 = sum(bins_state1$gene == gene)
        failures1 = nrow(bins_state1) - success1

        bins_state2 = df %>% dplyr::filter(domain == state2)
        success2 = sum(bins_state2$gene == gene)
        failures2 = nrow(bins_state2) - success2
        x = rbind(c(success1, failures1),
                  c(success2, failures2))

        test = stats::prop.test(x=x,alternative = "greater")
        triplets[[k]] = list(gene_name = gene, state1 = state1, state2 = state2, p_val = test$p.value)

      }
      k=k+1
    }
  }
  return(triplets)
}



transcript_based_molecular_distinction <- function(transcript_data_with_domains, bin_data, meta_data){

  domain <- i <-  state1 <- state2 <-NULL

  num_domains = length(unique(meta_data$domain))
  selected_index = meta_data[,"index"]
  bin_data = bin_data[selected_index,]

  gene_names = colnames(bin_data)

  cores = parallel::detectCores()
  clust <- parallel::makeSOCKcluster(cores)
  doSNOW::registerDoSNOW(clust)
  result_triplet <- foreach::foreach(i = 1:length(gene_names), .combine = rbind,
                                     .packages = c('dplyr', 'data.table')) %dopar%
    gene_function(gene = gene_names[i],num_domains = num_domains, df = transcript_data_with_domains)
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


