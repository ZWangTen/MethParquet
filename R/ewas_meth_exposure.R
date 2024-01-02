#' Association analysis with methylation as exposure
#' @description ewas_meth_exposure computes association analysis for each CpG site using the null model
#' fitted by function \code{\link{NullModel}}
#'
#' @param db_obj MethList created using \code{\link{create_methlist}}.
#' @param m_null Null model created by function \code{\link{NullModel}}
#' @param select_sites Character string consisting of CpG(s) to be tested. Default: `full` (all CpGs).
#'  Only one selection for CpG is allowed, `sites`, `chromosome`, or `gene`.
#' @param select_chr Character string of chromosome number to select CpGs from. Default: `FALSE`.
#' @param gene_list Character string of gene names to select CpGs from. Default: `FALSE`.
#' @param gene_col Column for gene names in the CpG annotation data from `Methlist`.
#' @param NAs_to_zero Convert missing values to 0 in methylation data. Default: `FALSE`.
#' @param out_position Column name for CpG position in chromosome from `cpg_annot` of `MethList`.
#'  If given, output will contain these columns and chromosome number.
#' @param block_size Integer to specify number of CpGs in each iteration block. Default: 50,000.
#'
#' @return A data frame with coefficient estimates for CpGs as exposure and testing statistics.
#' @export
#' @import tidyverse
#' @examples
#' library(tidyverse)
#' library(arrow)
#' data(phenoData)
#' data(chrAnnotation)
#' data(MethData)
#' wdir <- getwd()
#' path <- paste0(wdir,'/Parquet_Directory')
#'
#' # Create Parquet data in 'path' and MethList
#' MethData %>% group_by(CHR) %>% arrow::write_dataset(path,format = "parquet")
#'
#' mlist <- create_methlist(db_path = path,cpg_col_db='CpG',subject_annot = phenoData,
#' subject_col_keep='all',cpgAnnot_col_keep=c(1:2,12:13,16),cpg_annot = chrAnnotation,
#' subject_id='sample_id',cpg_col_annot='Name', gene_col_name = 'UCSC_RefGene_Name')
#'
#' # Fit linear null model (default)
#' m.null=NullModel(db_obj=mlist,trait='bmi',covariates_string = c('age','sex'))
#'
#' # Testing association between each CpG and trait using fitted null model
#' lm_v=ewas_meth_exposure(db_obj=mlist,m.null,select_sites = 'full',select_chr = FALSE,
#' gene_list=FALSE,out_position='MAPINFO')
#' head(lm_v)
#' unlink(path,recursive=TRUE)

ewas_meth_exposure <- function(db_obj,m_null,select_sites='full',select_chr=FALSE,
                               gene_list=FALSE,gene_col=FALSE,NAs_to_zero=FALSE,
                               out_position=FALSE,block_size=50000){
  db <- db_obj$db
  pheno <- db_obj$subject_annot
  chr_dat <- db_obj$cpg_annot
  mod_n <- m_null$NullModel
  X <- m_null$X
  Y <- m_null$Y
  #######function to get res
  mod_test=function(cpg_list){
    db_CpG <- db%>%filter(CpG %in% cpg_list) %>% as.data.frame() %>%
      dplyr::select(-CHR) %>% column_to_rownames(var='CpG') %>% as.matrix()
    db_CpG <- t(db_CpG)
    if (nrow(X) != nrow(db_CpG)){
      message(paste( nrow(db_CpG) - nrow(X)), "samples have missing phenotypes, using", nrow(X), "samples")
      db_CpG <- db_CpG[rownames(X),]
    }
    if(NAs_to_zero) {
      db_CpG[is.na(db_CpG)] <- 0
    }
    if(length(class(mod_n)) ==1) {
      if (grepl('GENESIS', class(mod_n), fixed = TRUE)==TRUE) {
        C <- unique(diag(as.matrix(mod_n$cholSigmaInv)))
        CX <- mod_n$CX
        CXCXI <- mod_n$CXCXI
      } else {
        C <- as.numeric(chol(chol2inv(summary(mod_n)$sigma)))
        CX <- C*X
        CXCXI <- tcrossprod(CX, chol2inv(chol(crossprod(CX))))
      }
      qrmod <- base::qr(CX)
      Ytilde <- base::qr.resid(qrmod, as.matrix(C*Y)) #Ytilde <- CY - tcrossprod(CXCXI, crossprod(CY, CX))
      resid <- C*Ytilde  #resid <- m.null$residuals*C^2
      CG <- suppressWarnings(base::apply(db_CpG,2,function(i){C*i}))
      Gtilde <- CG - tcrossprod(CXCXI, crossprod(CG, CX))
      GPG <- colSums(Gtilde^2)
      score_SE <- sqrt(GPG)
      score <- as.vector(crossprod(db_CpG, resid))
      Stat <- score/score_SE # t-value
      Score_pval = pchisq(Stat^2, df = 1, lower.tail = FALSE)
      Estimate=score/GPG
      #Est.SE = 1/score.SE,
      res <- data.frame(CpG = colnames(db_CpG),Estimate = Estimate, Score = Stat,p_value = Score_pval) #
      rownames(res)<-NULL
      res<- res %>% mutate(fdr_bh= p.adjust(p_value, method = "BH"))
      return(res)
    } else { #glm.scoretest
      w <- mod_n$weights
      r <- mod_n$residuals
      dispersion <- 1
      ws <- sqrt(w)
      x2_1w <- base::qr.resid(mod_n$qr,ws*db_CpG)
      zw <- ws*r
      z_vec <- colSums(as.matrix(x2_1w*zw))/sqrt(colSums(as.matrix(x2_1w * x2_1w)))/sqrt(dispersion)
      p_val <- pchisq(z_vec^2,df=1,lower.tail = FALSE)
      res <- data.frame(CpG = colnames(db_CpG),Score = z_vec,p_value = p_val) #
      rownames(res)<-NULL
      res <- res %>% mutate(fdr_bh= p.adjust(p_value, method = "BH"))
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
             sum(select_sites!=FALSE)==length(select_sites) | isTRUE(all(gene_list!=FALSE))) {
    stop("Please specify only one argument at a time")
  }

  if(nrow(CpG_test) <= block_size) {
    res <- mod_test(cpg_list=CpG_test$CpG)
  } else{
    nblock <- ceiling(nrow(CpG_test)/block_size)
    blocks <- unname(split(1:nrow(CpG_test), cut(1:nrow(CpG_test), nblock)))
    res_list <- list()
    for(i in 1:length(blocks)){
      res_list[[i]] <- mod_test(cpg_list=CpG_test$CpG[blocks[[i]]])
    }
    res <- do.call(rbind, res_list)
  }
  if (!isFALSE(out_position)) {
    chr_dat <- db_obj$cpg_annot
    res <- merge(res,chr_dat[,c('CpG','CHR',out_position)])
  }
  return(res)
}
