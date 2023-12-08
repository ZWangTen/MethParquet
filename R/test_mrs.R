
#' Test Methylation Risk Score against a phenotype chosen
#'
#' @param db_obj Methlist containing the phenotype data (subject annotation)
#' created using \code{\link{create_methlist}}
#' Sample names in phenotype data should be matching sample names in `mrs`.
#' Set `db_obj` to `FALSE` if `phe_data` is specified.
#' @param phe_data Phenotype data in data frame format. Default: `FALSE`
#' Sample names in `phe_data` should be matching sample names in `mrs`
#' @param mrs Methylation risk score, like the matrix obtained from function \code{\link{dev_meth_score}}
#' @param outcome Column for outcome inside phenotype data.
#' @param covariates covariates to be included in the model. Default: `FALSE`
#'
#' @return
#' A list containing fitted model and coefficients estimates of methylation risk score (MRS).
#' For multinomial logistic regression, the list will also contain likelihood ratio test for MRS.
#' @export
#' @import tidyverse
#' @importFrom nnet multinom
#' @importFrom car Anova
#' @examples
#' library(tidyverse)
#' library(arrow)
#' library(nnet)
#' library(car)
#' data(phenoData)
#' data(chrAnnotation)
#' data(ewas_bmi)
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
#' # Calculation of MRS associated with BMI
#' mrs_bmi <- dev_meth_score(db_obj=mlist,cpg_col='CpG',meth_effect_data = ewas_bmi,effect_col='Effect')
#'
#' height <- test_mrs(db_obj = mlist, outcome = 'height', covariates = c('age','sex'), mrs = mrs_bmi)
#' height[[2]]
#' unlink(path,recursive=TRUE)

test_mrs <- function(db_obj, phe_data=FALSE, mrs, outcome, covariates=FALSE){

  # Set phenotype data
  if (isFALSE(phe_data)==FALSE) {
    phen <- phe_data
  } else {
    phen <- db_obj$subject_annot
  }

  # Check if outcome and covariates are present in phenotype data, stop if not
  stopifnot('Your outcome variable is not in your phenotype data'=outcome %in% colnames(phen))

  # Check if sample names match in mrs and phen
  stopifnot('Samples in your phenotype data do not match those in your mrs '= suppressWarnings(isTRUE(all(row.names(phen) == row.names(mrs)))) )

  # Add mrs to phen
  phen <- phen %>% mutate(mrs=as.numeric(mrs))

  # Run regression:
  if (class(phen[,outcome])=='numeric' & !all(na.omit(phen[,outcome]) %in% 0:1)) {
    if (length(covariates)==1 & isFALSE(covariates)==TRUE) {
      x <- model.matrix(as.formula(paste(c('~1','mrs'),collapse='+')),data=phen)
      y <- phen %>% filter(!is.na(phen[,outcome])) %>% dplyr::select(outcome) %>% as.matrix()
      if (nrow(x) != nrow(y)){
        message(paste( nrow(x) - nrow(y), "samples have missing phenotypes, using", nrow(y), "samples"))
        x <- model.matrix(as.formula(paste('~mrs+',paste(covariates,collapse="+"))),data = phen[rownames(y),])
      }
    } else {
      x <- model.matrix(as.formula(paste('~mrs+',paste(covariates,collapse="+"))),data = phen)
      y <- phen %>% filter(!is.na(phen[,outcome])) %>% dplyr::select(outcome) %>% as.matrix()
      x <- x[which(row.names(x)%in%row.names(y)),]
      y <- y[which(row.names(y)%in%row.names(x)),]
    }
    mod <- lm(y~x)
    coef_df <- as.data.frame(summary(mod)$coef[2,])
    colnames(coef_df)<-'mrs'
    res <- list(model_fitted=mod,coefficients=coef_df)
  } else if (all(na.omit(phen[,outcome]) %in% 0:1)) {
    if (length(covariates)==1 & isFALSE(covariates)==TRUE) {
      mod <- glm(as.formula(paste(outcome, "~mrs" )), data = phen, family = binomial)
    } else {
      mod <- glm(as.formula(paste(outcome, "~mrs+", paste(covariates,collapse = "+"))), data = phen, family = binomial)
    }
    coef_df <- as.data.frame(summary(mod)$coef[2,])
    colnames(coef_df)<-'mrs'
    res <- list(model_fitted=mod,coefficients=coef_df)
  } else if (class(phen[,outcome]) %in% c('character','factor')) {
    if (length(covariates)==1 & isFALSE(covariates)==TRUE) {
      mod <- multinom(as.formula(paste(outcome, "~mrs")), data = phen)
    } else {
      mod <- multinom(as.formula(paste(outcome, "~mrs+", paste(covariates,collapse = "+"))), data = phen)
    }
    coef_df <- as.data.frame(summary(mod)$coefficients[,2])
    colnames(coef_df)<-'mrs'
    # Likelihood ratio test
    anov <- car::Anova(mod,type = "III")
    res <- list(model_fitted=mod,coefficients=coef_df,likelihoodratio_test=anov)
  }
  return(res)
}

