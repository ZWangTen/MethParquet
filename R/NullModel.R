#' Fit regression model under Null Hypothesis
#' @description NullModel fits a linear, generalized or mixed regression model with random effects
#' specified by covariance matrices (See \code{\link[GENESIS]{fitNullModel}}) such as a kinship matrix.
#' The output of NullModel can be passed to ewas_meth_exposure for association testing.
#'
#' @param db_obj MethList created using \code{\link{create_methlist}}.
#' @param trait Outcome phenotype in `subject_annot` from `MethList`.
#' @param covariates_string Columns for covariates to adjust in model from `subject_annot` of `MethList`.
#' @param cov.mat A matrix or list of matrices for random effect in mixed linear model. Default: NULL.
#' @param method `lm` for linear regression, `glm` for generalized linear regression,
#' `mixed` for mixed linear regression. Default: `lm`
#' @param family Which model to fit for generalized and mixed regression models.
#' See \code{\link[stat]{glm}} and \code{\link[GENESIS]{fitNullModel}}
#'
#' @return A list containing fitted null model, X and Y matrices, where Y is the phenotypic trait
#' and X created based on `covariates_string`.
#' @export
#' @import dplyr
#' @importFrom GENESIS fitNullModel
#' @examples
#' data(phenoData)
#' data(chrAnnotation)
#' data(ewas_bmi)
#' wdir <- getwd()
#' path <- paste0(wdir,'/data/Parquet_Directory')
#'
#' # Create MethList
#' mlist <- create_methlist(db_path = path,cpg_col_db='CpG',subject_annot = phenoData,
#' subject_col_keep='all',cpgAnnot_col_keep=c(1:2,12:13,16),cpg_annot = chrAnnotation,
#' subject_id='sample_id',cpg_col_annot='Name', gene_col_name = 'UCSC_RefGene_Name')
#'
#' # Fit linear null model (default)
#' m.null=NullModel(db_obj=mlist,trait='bmi',covariates_string = c('age','sex'))
#' names(m.null)
#' m.null$NullModel

NullModel <- function(db_obj, trait, covariates_string=NULL,cov.mat=NULL,
                         method='lm',family=FALSE){
  pheno <- db_obj$subject_annot

  stopifnot(is.element(trait, colnames(pheno)))
  y <- pheno[,trait]

  if(is.null(covariates_string)==FALSE) {
    model_string <- paste(c(covariates_string),collapse="+")
    x_covariates <- model.matrix(as.formula(paste('~',model_string)),data=pheno)
    if (nrow(x_covariates) != length(y)){
      message(paste( length(y) - nrow(x_covariates), "samples have missing phenotypes, using",
                     nrow(x_covariates), "samples"))
      y <- pheno[rownames(x_covariates),trait]
    }
    if (method=='lm' & is.null(cov.mat)) {
      mod_null <- lm(y~-1+x_covariates)
    } else if (method=='glm' & is.null(cov.mat)) {
      mod_null <- suppressWarnings(glm(y~-1+x_covariates,family = family))
    } else if (method=='mixed') {
      mod_null <- GENESIS::fitNullModel(pheno, outcome = trait, covars = covariates_string,
                                        cov.mat = cov.mat, family = family)
    }
    return(list(NullModel=mod_null,Y=y,X=x_covariates))
  } else {
    x_covariates <- model.matrix(as.formula(paste('~',1)),data=pheno)
    if (method=='lm'){
      mod_null <- lm(y~-1+x_covariates) # model with only intercept
    } else if (method=='glm' & is.null(cov.mat)) {
      mod_null <- suppressWarnings(glm(y~-1+x_covariates,family = family))
    } else if (method=='mixed') {
      mod_null <- GENESIS::fitNullModel(pheno, outcome = trait, covars = covariates_string,
                                        cov.mat = cov.mat, family = family)
    }
    return(list(NullModel=mod_null,Y=y,X=x_covariates))
  }

}
