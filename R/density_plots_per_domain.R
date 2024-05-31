#' Total number of transcripts in each domain
#'
#' This function returns a boxplot for total number of transcript in each domain. Each point of boxplot is a cell.
#' @param transcript_data_with_domains
#'
#' @return
#' @export
#'
#' @examples
density_plots_per_domain <- function(transcript_data_with_domains){
  count <- domain <- gene <- uID <- NULL
  temp = transcript_data_with_domains %>% dplyr::group_by(domain, uID) %>% dplyr::summarise(count = length(gene))
  p1 = ggplot2::ggplot(data = temp, ggplot2::aes(x = factor(domain), y = count)) + ggplot2::geom_boxplot()+
    ggplot2::xlab("domains") + ggplot2::ylab("number of transcripts per bin")+
    ggplot2::ggtitle("# transcripts in each bin")

  print(p1)
  return(p1)
}
