#' Extract subset of CpG sites
#' @description cpg_extract extracts subset of CpG methylation data based on argument provided
#'
#' @param db_obj MethList created using \code{\link{create_methlist}}
#' @param gene_list Character string of gene names to select CpGs from. Default: `FALSE`.
#' @param gene_col Column for gene names in the CpG annotation data from `Methlist`.
#' @param select_chr Character string of chromosome number to select CpGs from. Default: `FALSE`.
#' @param chr_col Column for chromosome in `cpg_annot` of the Methlist. Default: `FALSE`
#' @param site_names Character string for CpG names
#' @param row_indices Numeric vector indicating the row indices to extract CpG data
#'
#' @return Subset of CpG data along with their annotation based on selection criteria
#' @export
#' @import dplyr
#' @import rlang
#' @examples
#' library(tidyverse)
#' library(arrow)
#' data(phenoData)
#' data(chrAnnotation)
#'
#' wdir <- getwd()
#' methpath <- system.file('extdata','MethData.csv',package='MethParquet')
#' data(phenoData)
#' data(chrAnnotation)
#'
#' # Create Parquet data and MethList
#' path <- paste0(wdir,'/Parquet_Directory')
#' write_parquet_meth(data_path=methpath,format='csv',group_by='CHR',parquet_path = path)
#'
#' mlist <- create_methlist(db_path = path,cpg_col_db='CpG',subject_annot = phenoData,
#' subject_col_keep='all',cpgAnnot_col_keep=c(1:2,12:13,16),cpg_annot = chrAnnotation,
#' subject_id='sample_id',cpg_col_annot='Name', gene_col_name = 'UCSC_RefGene_Name')
#'
#' # Extract CpG data based on chromosome
#' chr <- as.character(seq(1:5))
#' extr_chr <- cpg_extract(mlist,select_chr=chr,chr_col='CHR')
#' head(extr_chr[[1]])[1:5,1:5]
#' head(extr_chr[[2]])
#'
#' # Based on row index
#' extr_row <- cpg_extract(mlist,row_indices = 1:100)
#' unlink(path,recursive=TRUE)

cpg_extract <- function(db_obj, gene_list=FALSE,gene_col=FALSE, select_chr=FALSE, chr_col=FALSE, site_names=FALSE, row_indices=FALSE) {
  db <- db_obj$db
  annot <- db_obj$cpg_annot

  # check how many inputs there are
  input <- c(paste(gene_list,collapse=','),paste(select_chr,collapse=','),paste(site_names,collapse=','),paste(row_indices,collapse=','))
  variable <- input[which(input!='FALSE')]
  stopifnot('Please specify one condition at a time' <- isTRUE(length(variable)==1))

  if (isTRUE(all(gene_list!=FALSE))) {
    stopifnot('Please enter gene column name in CpG annotation data' <- isFALSE(gene_col==FALSE))
    annot <- annot %>% rename(Gene = as.character(gene_col))
    annot$Gene <- paste0(";", annot$Gene, ";")
    genes <- paste0(";", genes, ";")
    gene_annot <- filter(annot, grepl(paste(genes, collapse="|"), annot$Gene))
    CpGs <- db %>% filter(CpG %in% gene_annot$CpG) %>% as.data.frame()
    res<-list(cpg_data=CpGs,cpg_annotation=gene_annot)
  } else if (isTRUE(all(select_chr!=FALSE))) {
    stopifnot('Please enter chromosome column name in CpG annotation data'=isFALSE(chr_col==FALSE))
    annot <- annot %>% rename(CHR = as.character(chr_col))
    gene_annot <- annot %>% filter(CHR%in%select_chr)
    CpGs <- db %>% filter(CpG %in% gene_annot$CpG) %>% as.data.frame()
    res<-list(cpg_data=CpGs,cpg_annotation=gene_annot)
  } else if (isTRUE(all(site_names!=FALSE))) {
    gene_annot <- annot %>% filter(CpG%in%site_names)
    CpGs <- db %>% filter(CpG %in% gene_annot$CpG) %>% as.data.frame()
    res<-list(cpg_data=CpGs,cpg_annotation=gene_annot)
  } else if (all(row_indices)!=FALSE) {
    CpGs <- db %>% as.data.frame() %>% rownames_to_column(var = "row") %>% filter(row %in% row_indices)
    gene_annot <- annot %>% filter(CpG%in%CpGs$CpG)
    res<-list(cpg_data=CpGs,cpg_annotation=gene_annot)
  }

  return(res)

}
