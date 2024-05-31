#' Assign domain to each trasncript for a cell
#'
#' @param cell
#' @param binned
#' @param binsize
#'
#' @return
#' @export
#'
#' @examples
merge_function <- function(cell, binned, binsize){

  x <- y <- cut_x <- cut_y <- node_center_x <-  node_center_y <- center_x <- center_y <- uID <- NULL

  xmin=  min(cell$x)
  xmax = max(cell$x)
  ymin = min(cell$y)
  ymax = max(cell$y)
  cell = cell %>% dplyr::mutate(cut_x = cut(x, breaks = seq(xmin-binsize,xmax+binsize,by = binsize),dig.lab = 9),
                         cut_y = cut(y, breaks = seq(ymin-binsize,ymax+binsize,by = binsize),dig.lab = 9))

  cell2 = cell %>% dplyr::group_by(uID,cut_x,cut_y) %>%  dplyr::mutate("center_x" = center_calc_v(cut_x)) %>%
    dplyr::mutate("center_y" = center_calc_v(cut_y)) %>% dplyr::ungroup() %>%
    dplyr::select(!c("cut_x","cut_y","uID")) %>%
    dplyr::rename(node_center_x = center_x, node_center_y = center_y)

  cell2$node_center_y = ceiling(cell2$node_center_y)
  cell2$node_center_x = ceiling(cell2$node_center_x)
  cell2 = cell2 %>% dplyr::mutate(coords = paste(node_center_x,"-",node_center_y))

  binned$node_center_x = ceiling(binned$node_center_x)
  binned$node_center_y = ceiling(binned$node_center_y)

  binned = binned %>% dplyr::mutate(coords = paste(node_center_x,"-",node_center_y))

  temp = merge(cell2, binned, by = c("coords")) %>% dplyr::select(c("x","y","gene","uID","domain"))

  return(temp)
}




#' Assign domain to each transcript
#'
#' This function merges the transcript data and the discovered domain data. All the transcripts in a bin are assigned the domain of that bin.
#'
#' @param meta_data
#' @param transcript_data
#' @param binsize
#'
#' @return
#' @export
#'
#' @examples
transcripts_with_assigned_domains <- function(meta_data, transcript_data, binsize){
  i <- NULL

  uid_split1 <- split(transcript_data, transcript_data$uID)
  uid_split2 <- split(meta_data, meta_data$uID)

  cores = parallel::detectCores()
  clust <- parallel::makeSOCKcluster(cores)
  doSNOW::registerDoSNOW(clust)

  merged_data <- foreach::foreach(i = 1:length(uid_split1), .combine = rbind,
                         .packages = c('dplyr', 'data.table')) %dopar%
    ((merge_function(cell = uid_split1[[i]], binned = uid_split2[[i]],binsize = binsize)))
  parallel::stopCluster(clust)
  colnames(merged_data) = c("absX","absY","gene","uID","res")



  return(merged_data)

}
