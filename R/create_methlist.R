
#' Create MethList
#' @description create_methlist creates the MethList object containing connection to parquet database,
#' and the annotation for CpGs and samples in data.frame.
#'
#' @param db_path Filepath to a methylation parquet database. Should contain a column named 'CpG'
#' for CpG site names
#' @param cpg_col_db Column for CpG site names in parquet data, which will be named as `CpG` in
#' CpG data inside meth list
#' @param subject_annot Subject annotation file
#' @param subject_col_keep Numeric vector containing columns to include in the subject annotation
#' file. Default value: `all`
#' @param subject_id Column for subject IDs in `subject_annot`, which will be named as `sample_id`
#' in `subject_annot` inside meth list
#' @param cpg_annot CpG annotation file
#' @param cpg_col_annot Column for CpG site names in `cpg_annot`
#' @param cpgAnnot_col_keep Numeric vector containing columns to include in the CpG annotation
#' file. Default value: `all`
#' @param gene_col_name Name of column containing gene names in `cpg_annot`. If povided, this column
#' will be named as `Gene` in `cpg_annot` of the meth list. Default value: `FALSE`
#' @param chr_col_name Name of column containing chromosome number in `cpg_annot`. If povided,
#' this column will be named as `CHR` in `cpg_annot` of the meth list. Default value: `FALSE`
#'
#' @return
#' A `MethList` object which contains:
#' `db`: open connection to a methylation parquet databse
#' `subject_annot`: subject annotation file
#' `cpg_annot`: CpG annotation file
#' @export
#' @import dplyr
#' @import rlang
#' @importFrom arrow open_dataset
#' @examples
#' library(tidyverse)
#' library(arrow)
#'
#' wdir <- getwd()
#' methpath <- system.file('extdata','MethData.csv',package='MethParquet')
#' data(phenoData)
#' data(chrAnnotation)
#'
#' # Create Parquet data in 'path'
#' path <- paste0(wdir,'/Parquet_Directory')
#' write_parquet_meth(data_path=methpath,format='csv',group_by='CHR',parquet_path = path)
#'
#' # Create MethList object
#' mlist <- create_methlist(db_path = path,cpg_col_db='CpG',subject_annot = phenoData,
#' subject_col_keep='all',cpgAnnot_col_keep=c(1:2,12:13,16),cpg_annot = chrAnnotation,
#' subject_id='sample_id',cpg_col_annot='Name', gene_col_name = 'UCSC_RefGene_Name')
#' names(mlist)
#' unlink(path,recursive=TRUE)


create_methlist <- function(db_path, cpg_col_db, subject_annot, subject_col_keep='all',
                            subject_id, cpg_annot, cpg_col_annot, cpgAnnot_col_keep='all',
                            gene_col_name = FALSE, chr_col_name=FALSE) {

  stopifnot('Subject_id is not in subject annotation'=subject_id %in% colnames(subject_annot))
  stopifnot('CpG column is not in CpG annotation'=cpg_col_annot %in% colnames(cpg_annot))

  subject_annot <- subject_annot %>% dplyr::rename(sample_id = as.character(subject_id))
  row.names(subject_annot) <- subject_annot$sample_id
  cpg_annot <- cpg_annot %>% dplyr::rename(CpG = as.character(cpg_col_annot))

  if(gene_col_name != FALSE) {
    cpg_annot <- cpg_annot %>% dplyr::rename(Gene = as.character(gene_col_name))
  }

  if(all(subject_col_keep!='all')==TRUE) {
    subject_annot <- subject_annot[,subject_col_keep]
  }

  if(chr_col_name != FALSE) {
    cpg_annot <- cpg_annot %>% dplyr::rename(CHR = as.character(chr_col_name))
  }

  if(all(cpgAnnot_col_keep!='all')==TRUE) {
    cpg_annot <- cpg_annot[,cpgAnnot_col_keep]
  }

  db <- open_dataset(db_path)
  db <- db %>% dplyr::rename(CpG=as.character(cpg_col_db))

  meth_list <- list(db=db, subject_annot=subject_annot, cpg_annot=cpg_annot)

  return(meth_list)
}

