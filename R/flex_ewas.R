#' Flexible implementation of association functions on methylation data
#' @description `flex_ewas` takes another regression function as input to test association between
#' each CpG site and the phenotype defined in the regression model
#'
#' @param db_obj MethList created using \code{\link{create_methlist}}.
#' @param fun A regression function applied to each CpG site
#' @param phen Phenotype data based on which the function (`fun`) is built. Default: `NULL` is the subject annotation from the `db_obj`
#' @param select_sites Character string consisting of CpG(s) to be tested. Default: `full` (all CpGs).
#' @param select_chr Character string of chromosome number to select CpGs from. Default: `FALSE`.
#' @param gene_list Character string of gene names to select CpGs from. Default: `FALSE`.
#' @param gene_col Column for gene names in the CpG annotation data from `Methlist`.
#' @param out_position Column names for CpG position in chromosome from `cpg_annot` of `MethList`.
#'  If given, output will contain these columns and chromosome number.
#' @param block_size Integer to specify number of CpGs in each iteration block. Default: 50,000.
#' @param NAs_to_zero Convert missing values to 0 in methylation data. Default: `FALSE`.
#'
#' @return A data frame with results obtained using `fun`.
#' @export
#' @import tidyverse
#' @examples
#' library(tidyverse)
#' library(arrow)
#' data(phenoData)
#' data(chrAnnotation)
#' wdir <- getwd()
#' path <- paste0(wdir,'/data/Parquet_Directory')
#'
#' # Create MethList
#' mlist <- create_methlist(db_path = path,cpg_col_db='CpG',subject_annot = phenoData,
#' subject_col_keep='all',cpgAnnot_col_keep=c(1:2,12:13,16),cpg_annot = chrAnnotation,
#' subject_id='sample_id',cpg_col_annot='Name', gene_col_name = 'UCSC_RefGene_Name')
#'
#' # Create a robust linear regression function with methylation as exposure
#' covariates_string = c('age','sex')
#' rlm_function <- function(x) {
#' mod <- base::suppressWarnings(MASS::rlm(as.formula(paste('body_fat', "~x+", covariates_string),data = phenoData)))
#' return(summary(mod)$coefficients[2,])
#' }
#'
#' # Run flex_ewas
#' flex_rlm <- flex_ewas(mlist,rlm_function,out_position = c('MAPINFO','Gene'))
#' head(flex_rlm)

flex_ewas <- function(db_obj,fun,phen=NULL,select_sites='full',select_chr=FALSE,gene_list=FALSE,
                      gene_col=FALSE,out_position=FALSE,block_size=50000,NAs_to_zero=FALSE){
  db <- db_obj$db
  if (is.null(phen)==TRUE) {
    phen <- db_obj$subject_annot
  }
  # Which CpGs to test
  if (sum(select_chr!=FALSE)==length(select_chr) & sum(select_sites==FALSE)==1 & sum(gene_list==FALSE)==1) {
    CpG_test <- chr_dat %>% filter(CHR %in% select_chr) %>% select(CpG,CHR) # Based on chromosome number
    CpG_test <- db %>% filter(CpG %in% CpG_test$CpG) %>% dplyr::select(CpG) %>% as.data.frame()
  } else if (sum(select_chr==FALSE)==1 & sum(select_sites!=FALSE)==1 & sum(gene_list==FALSE)==1) {
    CpG_test <- as.data.frame(db['CpG']) # Full data
  } else if (sum(select_sites!="full")==length(select_sites) & sum(select_chr==FALSE)==1 & sum(gene_list==FALSE)==1) {
    CpG_test <- as.data.frame(db['CpG']) %>% filter(CpG%in%select_sites) # Sites
  } else if (sum(select_chr==FALSE)==1 & sum(select_sites==FALSE)==1 & isTRUE(all(gene_list!=FALSE))) {
    stopifnot('Please enter gene column name in CpG annotation data' <- isFALSE(gene_col==FALSE)) # gene
    chr_dat <- chr_dat %>% rename(Gene = as.character(gene_col))
    chr_dat$Gene <- paste0(";", chr_dat$Gene, ";")
    genes <- paste0(";", gene_list, ";")
    gene_annot <- filter(chr_dat, grepl(paste(genes, collapse="|"), chr_dat$Gene))
    CpG_test <- db %>% filter(CpG %in% gene_annot$CpG) %>% dplyr::select(CpG) %>% as.data.frame()
  } else if (sum(select_chr==FALSE)==1 & sum(select_sites==FALSE)==1 & sum(gene_list==FALSE)==1) {
    stop("Please specify a selection argument for CpGs")
  } else if (sum(select_chr!=FALSE)==length(select_chr) &
             sum(select_sites!=FALSE)==length(select_sites) | isTRUE(all(gene_list!=FALSE))) {
    stop("Please specify only one argument at a time")
  }

  db_CpG <- db%>%filter(CpG %in% CpG_test$CpG) %>% as.data.frame() %>%
    dplyr::select(-CHR) %>% column_to_rownames(var='CpG') %>% as.matrix()
  if(NAs_to_zero) {
    db_CpG[is.na(db_CpG)] <- 0
  }
  db_CpG <- t(db_CpG)
  if (nrow(phen) != nrow(db_CpG)){
    message(paste( nrow(db_CpG) - nrow(phen), "samples have missing phenotypes, using", nrow(phen), "samples"))
    db_CpG <- db_CpG[rownames(phen),]
  }
  res <- as.data.frame(t(apply(db_CpG, 2, fun)))
  res <- rownames_to_column(res,var='CpG')
  if (!isFALSE(out_position)) {
    chr_dat <- db_obj$cpg_annot
    res <- merge(res,chr_dat[,c('CpG','CHR',out_position)])
  }
  return(res)
}
