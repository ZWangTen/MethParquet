#' Calculation of Methylation Risk Score (MRS)
#' @description dev_meth_score calculates Methylation Risk Score as a weighted sum of an individual’s
#' DNA methylation values and pre-calculated weights associated with a specific phenotype.
#'
#' @param db_obj MethList created using \code{\link{create_methlist}}
#' @param meth_effect_data Phenotype data containing CpG sites and effect
#' @param effect_col Column for effect in `meth_effect_data`
#' @param cpg_col Column for CpG site names in `meth_effect_data`
#' @param scale_score Scale score using \code{\link[base]{scale}} function. Default: `FALSE`
#' @param divide_by_site_n Divide scores by number of methylation sites Default: `FALSE`
#' @param NAs_to_zero Convert any missing values in methylation data (inside Methlist) to 0. Default: `FALSE`
#'
#' @return A matrix of methylation risk scores associated to the phenotype for each subject
#' @export
#' @import dplyr
#' @import rlang
#' @examples
#' library(tidyverse)
#' library(arrow)
#' data(phenoData)
#' data(chrAnnotation)
#' data(ewas_bmi)
#'
#' wdir <- getwd()
#' methpath <- system.file('extdata','MethData.csv',package='MethParquet')
#'
#' # Create Parquet data and MethList
#' path <- paste0(wdir,'/Parquet_Directory')
#' write_parquet_meth(data_path=methpath,format='csv',group_by='CHR',parquet_path = path)
#'
#' mlist <- create_methlist(db_path = path,cpg_col_db='CpG',subject_annot = phenoData,
#' subject_col_keep='all',cpgAnnot_col_keep=c(1:2,12:13,16),cpg_annot = chrAnnotation,
#' subject_id='sample_id',cpg_col_annot='Name', gene_col_name = 'UCSC_RefGene_Name')
#'
#' # Calculation of MRS associated with BMI
#' mrs_bmi <- dev_meth_score(db_obj=mlist,cpg_col='CpG',meth_effect_data = ewas_bmi,effect_col='Effect')
#' head(mrs_bmi)
#' class(mrs_bmi)
#' unlink(path,recursive=TRUE)


dev_meth_score <- function(db_obj, meth_effect_data, cpg_col, scale_score = FALSE, effect_col,
                           divide_by_site_n = FALSE, NAs_to_zero = FALSE) {

  stopifnot('cpg_col is not in meth_effect_data'=cpg_col %in% colnames(meth_effect_data))
  meth_effect_data <- meth_effect_data %>% rename(CpG = as.character(cpg_col))
  meth_effect_data <- meth_effect_data %>% rename(Effect = as.character(effect_col))

  if(length(class(db_obj)) > 1) {
    ds <- db_obj
  } else {
    ds <- db_obj[[1]]
  }

  #filter parquet file by CpGs in 'meth_data' and collect
  available_data <- ds %>%
    filter(CpG %in% meth_effect_data$CpG) %>% as.data.frame()

  #filter meth_data to exclude missing CpGs
  meth_data <- meth_effect_data %>% filter(CpG %in% available_data$CpG)

  #match order of meth_data with available_data
  meth_data <- meth_data[match(available_data$CpG, meth_data$CpG),]

  #change rownames in available_data to CpG
  rownames(available_data) <- available_data$CpG

  #remove CpG and chr columns from available_data, then convert to matrix, and flip so CpGs are columns
  if("CHR" %in% colnames(available_data)){
    available_data <- available_data %>%
      dplyr::select(-CpG, -CHR) %>%
      as.matrix() %>% t()
  } else {
    available_data <- available_data %>%
      dplyr::select(-CpG) %>%
      as.matrix() %>% t()
  }

  #replace NAs with 0 if selected
  if(NAs_to_zero) {
    available_data[is.na(available_data)] <- 0
  }


  #calculate methylation risk score
  meth_score <- available_data[,meth_data$CpG] %*% meth_data$Effect

  #option to scale score
  if(scale_score) {
    meth_score <- scale(meth_score)
  }

  #option to divide score by number of CpG sites used
  if(divide_by_site_n){
    meth_score[,1] <- meth_score[,1] / nrow(meth_data)
  }

  #return scores
  return(meth_score)
}

