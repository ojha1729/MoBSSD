#' Title
#'
#' @param path_to_file
#' @param file_name
#'
#' @return
#' @export
#'
#' @examples
load_dataset <- function(path_to_file = NULL,file_name){
  if (is.null(path_to_file)){
    return(data.table::fread(file_name))
  }
  else{
    return(data.table::fread(paste0(path_to_file, file_name)))
  }

}
