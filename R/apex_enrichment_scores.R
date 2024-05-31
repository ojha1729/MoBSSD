#' APEX-seq Enrichment Scores
#'
#' This function calculates enrichment scores using decoupleR package
#' @param weights_file
#' @param exp_data
#' @param meta_data
#'
#' @return
#' @export
#'
#' @examples
enrichment_score <- function(weights_file, exp_data, meta_data){

source <- condition <- score <- domain <- statistic <- NULL
enrichment_df = decoupleR::run_wsum(mat = t(exp_data),network = weights_file,.source = "source",.target = "target",.mor = "weight")

normalised_scores = enrichment_df %>% dplyr::filter(statistic == "norm_wsum") %>% dplyr::select(source, condition, score)
temp = reshape2::dcast(normalised_scores, condition~source)
temp$condition = as.integer(temp$condition)
temp = temp %>% dplyr::arrange(condition)
temp = temp %>% dplyr::mutate(domain = meta_data$domain)

num_domains = length(unique(meta_data$domain))

enrich_vec = list()
for (val in seq(1,num_domains)) {

  df = temp %>% dplyr::filter(domain == val) %>% dplyr::select(!domain)
  df = df>0
  enrich_vec[[val]] = colMeans(df)
}

res = data.frame(do.call(rbind, enrich_vec)) %>% dplyr::mutate(domain = seq(1,num_domains))

return(res)
}




#' Enrichment Score Heatmap
#'
#' @param enrichment_score_df
#'
#' @return
#' @export
#'
#' @examples
enrichment_score_plotter <- function(enrichment_score_df){
domain <- locale <- score <- NULL
  dat <- reshape2::melt(enrichment_score_df,id.vars = "domain")
  colnames(dat) = c("domain","locale","score")
  # Heatmap with annotations
  p1 <- ggplot2::ggplot(dat, ggplot2::aes(domain, locale)) +
    ggplot2::geom_tile(ggplot2::aes(fill = score)) +
    ggplot2::scale_fill_gradient(low = "white", high = "red") +
    ggplot2::geom_text(ggplot2::aes(label = round(score, 2)), size = 3)

  # Display the plot
  print(p1)

}

