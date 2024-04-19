
#' Association analysis with linear regression
#' @description lm_ewas_outcome performs association analysis with linear regression
#' using methylation as the outcome.

#' @param db_obj MethList created using \code{\link{create_methlist}}.
#' @param trait Exposure phenotype in `subject_annot` from `MethList`.
#' @param select_sites Character string consisting of CpG(s) to be tested. Default: `full` (all CpGs).
#'  Only one selection for CpG is allowed, `sites`, `chromosome`, or `gene`.
#' @param select_chr Character string of chromosome number to select CpGs from. Default: `FALSE`.
#' @param gene_list Character string of gene names to select CpGs from. Default: `FALSE`.
#' @param gene_col Column for gene names in the CpG annotation data from `Methlist`.
#' @param covariates_string Columns for covariates to adjust in model from `subject_annot` of `MethList`.
#' @param out_position Column name for CpG position in chromosome from `cpg_annot` of `MethList`.
#'  If given, output will contain these columns and chromosome number.
#' @param NAs_to_zero Convert missing values to 0 in methylation data. Default: `FALSE`.
#' @param block_size Integer to specify number of CpGs in each iteration block. Default: 50,000.
#' @param parallel Number of cores used for parallel processing (forking). For details please see
#' \code{\link[future]{plan}}.  Default: FALSE (sequential processes).
#'
#' @return A data frame with coefficient estimates for `trait` and testing statistics.
#' @export
#' @import tidyverse Rcpp furrr future RcppEigen
#' @importFrom Rcpp sourceCpp
#' @useDynLib MethParquet
#' @exportPattern "^[[:alpha:]]+"
#' @examples
#' library(tidyverse)
#' library(arrow)
#'
#' wdir <- getwd()
#' methpath <- paste0(wdir,'/inst/extdata/MethData.csv')
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
#' # Run linear regression model
#' ewas_lm <- lm_ewas_outcome(db_obj=mlist,trait='bmi',select_sites='full',select_chr=FALSE,
#' covariates_string = c('age','sex'),out_position='MAPINFO',block_size=100000)
#' head(ewas_lm)
#' unlink(path,recursive=TRUE)

lm_ewas_outcome <-function(db_obj, trait, select_sites='full',select_chr=FALSE,
                             gene_list=FALSE,gene_col=FALSE, covariates_string, parallel=FALSE,
                             out_position=FALSE, NAs_to_zero=FALSE, block_size=50000){
  db <- db_obj$db
  pheno <- db_obj$subject_annot
  chr_dat <- db_obj$cpg_annot

  # lm_ewas function
  lewas <- function(cpg_list){
    db_CpG<-db%>%filter(CpG %in% cpg_list) %>% as.data.frame() %>%
      dplyr::select(-CHR) %>% column_to_rownames(var='CpG') %>% as.matrix()

    if(NAs_to_zero) {
      db_CpG[is.na(db_CpG)] <- 0
    }
    db_CpG <- t(db_CpG) # Now row-samples, col-CpGs

    if (nrow(XX) != nrow(db_CpG)){
      message(paste( nrow(db_CpG) - nrow(XX), "samples have missing phenotypes, using", nrow(XX), "samples"))
      db_CpG <- db_CpG[rownames(XX),]
    }

    if (!class(pheno[,trait])%in%c('character','factor')) {
      # β = [(Xt*X)^−1]*Xt*y
      n.trait<-which(colnames(XX)==trait)-1
      lm.b <- calculateBetasAndSE(XX, db_CpG, n.trait)

      test_stats <- lm.b$betas/lm.b$se_betas
      t_stat_df <- nrow(db_CpG) - ncol(XX)
      t_pval <- 2*pt(abs(test_stats), lower.tail=FALSE, df = t_stat_df)
      res <- data.frame(CpG = colnames(db_CpG), estimate = lm.b$betas, se = lm.b$se_betas,
                        t_stat = test_stats,t_stat_df=t_stat_df, p_value = t_pval)

      rownames(res)<-NULL
      return(res)
    } else {          # Categorical trait
      ss_res_full <- solveAndCalculateResiduals(XX, db_CpG)

      # reduced model
      if (isFALSE(covariates_string)==F){
        model_string2 <- paste(c(covariates_string),collapse="+")
        XX2<-model.matrix(as.formula(paste('~',model_string2)),data=pheno)
      } else {
        XX2 <- model.matrix(as.formula(paste('~',1)),data=pheno)
      }
      XX2 <- as.matrix(XX2[rownames(XX),])
      ss_res_red <- solveAndCalculateResiduals(XX2, db_CpG)

      # F test
      df_full <- nrow(XX)-nrow(b_full)
      df_red <- nrow(XX2)-nrow(b_red)
      F_val <- ((ss_res_red - ss_res_full)/(df_red-df_full))/(ss_res_full/df_full)
      p_val <- pf(F_val, df_red-df_full, df_full, lower.tail = FALSE)
      res <- data.frame(CpG = colnames(db_CpG), F_estimate=F_val,  p_value = p_val)
      rownames(res)<-NULL
      return(res)
    }
  }

  # Select CpGs:
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
             sum(select_sites!=FALSE)==length(select_sites) & isTRUE(all(gene_list!=FALSE))) {
    stop("Please specify only one argument at a time")
  }

  # multiple or no covariates:
  if (isFALSE(covariates_string)==FALSE) {
    model_string <- paste(c(trait,covariates_string),collapse="+")
    XX<-model.matrix(as.formula(paste('~',model_string)),data=pheno)
  } else{
    XX <- model.matrix(as.formula(paste(c('~1',trait),collapse='+')),data=pheno)
  }

  if (nrow(CpG_test) <= block_size) {
    res<-lewas(CpG_test$CpG)
    res<- res %>% mutate(fdr_bh= p.adjust(p_value, method = "BH"))
  } else{
    nblock  <- rep(1:ceiling(nrow(CpG_test)/block_size),each=block_size)[1:nrow(CpG_test)]
    blocks <- split(CpG_test$CpG,nblock)
    if (isFALSE(parallel)==TRUE) {
     res_list <- lapply(blocks,function(x) {lewas(cpg_list=x)})
      # res_list <- list()
      #for(i in 1:length(blocks)){
     #   res_list[[i]]<-lewas(cpg_list=blocks[[i]])
      #}
    } else{
      oplan<-plan(multicore,workers=parallel) #
      res_list <- future_map(blocks, ~lewas(.x))
      on.exit(plan(oplan), add = TRUE)
    }

    res <- do.call(rbind, res_list)
    res<- res %>% mutate(fdr_bh= p.adjust(p_value, method = "BH"))
  }
  if (!isFALSE(out_position)) {
    chr_dat <- db_obj$cpg_annot
    res <- merge(res,chr_dat[,c('CpG',out_position)])
  }
  return(res)
}


