#' Write Parquet database from processed methylation data
#' @description write_parquet_meth writes Parquet database from methyltion data without loading it into
#' local machine. For a seamless MethParquet workflow, this function requires that the methylation data
#' contains a column representing chromosome number.
#' Currently supports comma-separated values (CSV), delimited text file (txt) and tab-separated values (TSV) files.
#' For more details please refer to \code{\link[arrow]{write_dataset}}.
#'
#' @param data_path Path to methylation data.
#' @param format A string identifier of the file format, can be one of the following: ('csv','txt','tsv').
#' @param chunk_size Maximum number of rows per file in Parquet database. A value larger than 0 indicates
#' the number of rows of each single Parquet file. Default: 0L.
#' @param chr_col Character string indicating the column name for the chromosome number. Default: 'CHR'
#' @param group_by The variable (column name) of the methylation data to partition the Parquet data (as path segments).
#' @param parquet_path Path to write the Parquet directory
#'
#' @return A Parquet directory for the methylation data
#' @export
#' @import arrow
#'
#' @examples
#' library(arrow)
#' library(tidyverse)
#' wdir <- getwd()
#' methpath <- paste0(wdir,'/inst/extdata/MethData.csv')
#' path <- paste0(wdir,'/Parquet_Directory')
#' write_parquet_meth(data_path=methpath,format='csv',group_by='CHR',parquet_path = path)
#' list.files(path)
#' unlink(path,recursive=TRUE)

write_parquet_meth <- function(data_path, format=c('csv','txt','tsv'), chunk_size=0,
                               chr_col='CHR', group_by=NULL, parquet_path){

  if(format=='csv') {
    data <- read_csv_arrow(data_path, as_data_frame = FALSE)
    stopifnot('Please include the chromosome column in methylation data' <- isTRUE(chr_col%in%data$schema$names))
  }
  if(format=='txt') {
    data <- read_delim_arrow(data_path, as_data_frame = FALSE)
    stopifnot('Please include the chromosome column in methylation data' <- isTRUE(chr_col%in%data$schema$names))
  }
  if(format=='tsv') {
    data <- read_tsv_arrow(data_path, as_data_frame = FALSE)
    stopifnot('Please include the chromosome column in methylation data' <- isTRUE(chr_col%in%data$schema$names))
  }

  write_dataset(data, path=parquet_path, format = 'parquet',max_rows_per_group = bitwShiftL(1, 20),
                partitioning = group_by,max_rows_per_file = chunk_size)
}
