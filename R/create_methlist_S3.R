#' Create MethList from Amazon Web Services S3 bucket
#' @description create_methlist_S3 creates the MethList object linking parquet database stored in
#' Amazon Web Services(AWS) Simple Storage Service (S3) bucket.
#'
#' @param db_s3_path Filepath to the methylation parquet database in AWS_S3 bucket. Should contain a column
#' named 'CpG' for CpG site names
#' @param subject_annot Subject annotation file
#' @param subject_path Filepath to the subject annotation file in AWS_S3 bucket if not provided as
#' `subject_annot`. Default value: `FALSE`
#' @param cpg_annot CpG annotation file
#' @param cpg_path Filepath to the CpG annotation file in AWS_S3 bucket if not provided as
#' `cpg_annot`. Default value: `FALSE`
#' @param subject_readFun Function to read the subject annotation file in AWS_S3 bucket.
#' @param cpg_readFun Function to read the cpg annotation file in AWS_S3 bucket.
#' @param access_key access key for the AWS_S3 bucket, if applicable.
#' @param secret_access_key secret access key for the AWS_S3 bucket, if applicable.
#' @param default_region AWS default region, if applicable.
#' @param cpg_col_db Column for CpG site names in parquet data, which will be named as `CpG` in
#' CpG data inside meth list
#' @param subject_col_keep Numeric vector containing columns to include in the subject annotation
#' file. Default value: `all`
#' @param subject_id Column for subject IDs in `subject_annot`, which will be named as `sample_id`
#' in `subject_annot` inside meth list
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
#' `db`: open connection to a methylation parquet databse in AWS_S3
#' `subject_annot`: subject annotation file
#' `cpg_annot`: CpG annotation file
#' @export
#' @import tidyverse aws.s3
#' @importFrom arrow open_dataset
#' @examples
#' # Example code only
#'
#' # mlist=create_methlist_S3(db_s3_path = "path to methylation Parquet data in AWS_S3",
#' # subject_path="path to subject annotation file in AWS_S3",
#' # cpg_path="path to CpG annotation file in AWS_S3",
#' # access_key = "access key to S3 bucket", secret_access_key = "secret access key to S3 bucket",
#' # default_region = "us-east-1", cpgAnnot_col_keep=c(1:2,12:13,16),
#' # subject_readFun = load, cpg_readFun = read.csv, cpg_col_db='CpG',
#' # subject_id='sample_id', cpg_col_annot='Name', gene_col_name = 'UCSC_RefGene_Name')

create_methlist_S3 <- function(db_s3_path, subject_annot=FALSE, subject_path=FALSE, cpg_annot=FALSE,
                               cpg_path=FALSE, subject_readFun=FALSE, cpg_readFun=FALSE, access_key=FALSE,
                               secret_access_key=FALSE, default_region=FALSE, cpg_col_db,
                               subject_col_keep='all', subject_id, cpg_col_annot, cpgAnnot_col_keep='all',
                               gene_col_name = FALSE, chr_col_name=FALSE) {

  if(access_key != FALSE) {
    Sys.setenv("AWS_ACCESS_KEY_ID" = access_key)
  }
  if(secret_access_key != FALSE) {
    Sys.setenv("AWS_SECRET_ACCESS_KEY" = secret_access_key)
  }
  if(default_region != FALSE) {
    Sys.setenv("AWS_DEFAULT_REGION" = default_region)
  }

  db <- open_dataset(db_s3_path)
  db <- db %>% dplyr::rename(CpG=as.character(cpg_col_db))

  loading <- function(dataset){
    merchants <- load(dataset)
    return(get(merchants))
  }

  if(isFALSE(subject_annot)) {
    subject_annot <- s3read_using(FUN = subject_readFun, object = subject_path)
    if(class(subject_annot)=='character') {
      subject_annot <- s3read_using(FUN = loading, object = subject_path)
    }
  }
  if(isFALSE(cpg_annot)) {
    cpg_annot <- s3read_using(FUN = cpg_readFun, object = cpg_path)
    if(class(cpg_annot) == 'character') {
      cpg_annot <- s3read_using(FUN = loading, object = cpg_path)
    }
  }

  stopifnot('Subject_id is not in subject annotation'=subject_id %in% colnames(subject_annot))
  stopifnot('CpG column is not in CpG annotation'=cpg_col_annot %in% colnames(cpg_annot))

  subject_annot <- subject_annot %>% dplyr::rename(sample_id = as.character(subject_id))
  row.names(subject_annot) <- subject_annot$sample_id
  cpg_annot <- cpg_annot %>% dplyr::rename(CpG = as.character(cpg_col_annot))

  if(gene_col_name != FALSE) {
    cpg_annot <- cpg_annot %>% dplyr::rename(Gene = as.character(gene_col_name))
  }
  if(chr_col_name != FALSE) {
    cpg_annot <- cpg_annot %>% dplyr::rename(CHR = as.character(chr_col_name))
  }

  if(all(subject_col_keep!='all')==TRUE) {
    subject_annot <- subject_annot[,subject_col_keep]
  }
  if(all(cpgAnnot_col_keep!='all')==TRUE) {
    cpg_annot <- cpg_annot[,cpgAnnot_col_keep]
  }

  meth_list <- list(db=db, subject_annot=subject_annot, cpg_annot=cpg_annot)
  return(meth_list)
}
