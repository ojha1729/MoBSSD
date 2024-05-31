#' Title
#'
#' @param group A dataframe of coordinates of a cell with distance from nuclear membrane column named as "inNucleus"
#'
#' @return Returns 0 if bins are inside
#' @export
#'
#' @examples
nucleus_marker <- function(group){
  if (length(unique(group$inNucleus))>1){
    return(-1)
  }
  else{
    return(unique(group$inNucleus))
  }
}

#' Title
#'
#' @param interval
#'
#' @return
#' @export
#'
#' @examples
center_calc <- function(interval){
  interval = as.character(interval)
  vals = as.numeric(unlist(regmatches(interval,gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",interval, perl=TRUE))))
  return(mean(vals))
}

center_calc_v <- function(column){
  return(unlist(lapply(1:length(column), function(i)center_calc(column[i]))))
}



#' Title
#'
#' @param group
#' @param gene_names_all
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom magrittr %>%
gene_dummy <- function(group, gene_names_all){
  gene <- NULL
  count_df = group %>% dplyr::group_by(gene) %>% dplyr::summarise(count = dplyr::n()) %>% dplyr::arrange(by = gene)
  dummy_vec = data.table::data.table(gene = sort(gene_names_all),count = rep(0,length(gene_names_all)))
  gene_names = group$gene
  dummy_vec[which(dummy_vec$gene %in% gene_names),"count"] = count_df$count
  vec = t(dummy_vec$count)
  colnames(vec) = sort(gene_names_all)

  return(vec)
}


#' Title
#'
#' @param cell
#' @param binsize
#' @param all_genes
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom magrittr %>%
cell_bins <- function(cell,binsize,all_genes){
  x <- y <- uID <- cut_x <- cut_y <- . <- NULL
  xmin=  min(cell$x)
  xmax = max(cell$x)
  ymin = min(cell$y)
  ymax = max(cell$y)

  cell = cell %>% dplyr::mutate(cut_x = cut(x, breaks = seq(xmin-binsize,xmax+binsize,by = binsize),dig.lab = 9),
                         cut_y = cut(y, breaks = seq(ymin-binsize,ymax+binsize,by = binsize),dig.lab = 9))

  cell2 = cell %>% dplyr::group_by(uID,cut_x,cut_y) %>% foreach::do(data.table::data.table(nucleus = nucleus_marker(.),
                                                               gene = gene_dummy(., gene_names_all = all_genes))) %>% dplyr::ungroup()
  cell3 = cell2 %>% dplyr::mutate(node_center_x = center_calc_v(cut_x)) %>%
    dplyr::mutate(node_center_y = center_calc_v(cut_y)) %>%
    dplyr::select(!c("cut_x","cut_y"))
  return(cell3)
}
