% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NullModel.R
\name{NullModel}
\alias{NullModel}
\title{Fit regression model under Null Hypothesis}
\usage{
NullModel(
  db_obj,
  trait,
  covariates_string = NULL,
  cov.mat = NULL,
  method = "lm",
  family = FALSE
)
}
\arguments{
\item{db_obj}{MethList created using \code{\link{create_methlist}}.}

\item{trait}{Outcome phenotype in `subject_annot` from `MethList`.}

\item{covariates_string}{Columns for covariates to adjust in model from `subject_annot` of `MethList`.}

\item{cov.mat}{A matrix or list of matrices for random effect in mixed linear model. Default: NULL.}

\item{method}{`lm` for linear regression, `glm` for generalized linear regression,
`mixed` for mixed linear regression. Default: `lm`}

\item{family}{Which model to fit for generalized and mixed regression models.
See \code{\link[stat]{glm}} and \code{\link[GENESIS]{fitNullModel}}}
}
\value{
A list containing fitted null model, X and Y matrices, where Y is the phenotypic trait
and X created based on `covariates_string`.
}
\description{
NullModel fits a linear, generalized or mixed regression model with random effects
specified by covariance matrices (See \code{\link[GENESIS]{fitNullModel}}) such as a kinship matrix.
The output of NullModel can be passed to ewas_meth_exposure for association testing.
}
\examples{
library(tidyverse)
library(arrow)
library(GENESIS)
data(phenoData)
data(chrAnnotation)

wdir <- getwd()
methpath <- system.file('extdata','MethData.csv',package='MethParquet')

# Create Parquet data and MethList
path <- paste0(wdir,'/Parquet_Directory')
write_parquet_meth(data_path=methpath,format='csv',group_by='CHR',parquet_path = path)

mlist <- create_methlist(db_path = path,cpg_col_db='CpG',subject_annot = phenoData,
subject_col_keep='all',cpgAnnot_col_keep=c(1:2,12:13,16),cpg_annot = chrAnnotation,
subject_id='sample_id',cpg_col_annot='Name', gene_col_name = 'UCSC_RefGene_Name')

# Fit linear null model (default)
m.null=NullModel(db_obj=mlist,trait='bmi',covariates_string = c('age','sex'))
names(m.null)
m.null$NullModel
unlink(path,recursive=TRUE)
}
