#' Read Dataset
#'
#' A function that reads dataset and modifies the 'inNucleus' column.
#' It is expected that coordinates are named as 'x' and 'y'.
#' Column name for genes must be 'gene'.
#' Column that contains distance from the nuclear membrane must be named 'inNucleus'.
#' @param path_to_file path to the folder where transcripts coordinates file is saved
#' @param file_name name of the file (passed as string)
#'
#' @return A csv file with coordinates of each transcript.
#' The 'inNucleus' column is updated.
#' It is changed to 0 for all the transcripts outside the nucleus (inNucleus = 0 in original file) and 1 for the trasncripts inside nucleus (inNucleus =1 in original file).
#' @export
#'
#' @examples
load_dataset <- function(path_to_file = NULL,file_name){
  inNucleus <- NULL
  if (is.null(path_to_file)){
    data_csv = data.table::fread(file_name)
    data_csv$inNucleus = ifelse(data_csv$inNucleus >=0, 0,1)
    return(data_csv)
  }
  else{
    data_csv = data.table::fread(paste0(path_to_file, file_name))
    data_csv$inNucleus = ifelse(data_csv$inNucleus >=0, 0,1)
    return(data_csv)
  }

}
