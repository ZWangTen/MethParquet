#' A mock Effect data for BMI phenotype to calculate Methylation Risk Score
#'
#' For demonstration purposes, some CpG names from original data are masked to match those in MethData.
#' The original Effect dataset for BMI is obtained from the EWAS Catalog.
#'
#' @docType data
#'
#' @usage data(ewas_bmi)
#'
#' @format ## `ewas_bmi`
#' A data frame with 477 rows and 7 columns:
#' \describe{
#'   \item{CpG}{CpG names}
#'   \item{Effect}{Effect estimate associated with BMI}
#'   \item{chr}{Chromosome}
#'   ...
#' }
#' @references Geurts et al. (2018) Int J Obes (Lond) 42(4):887-896
#' (\href{https://pubmed.ncbi.nlm.nih.gov/29278407/}{PubMed})
#'
#' @source <https://www.ewascatalog.org/?trait=bmi>
"ewas_bmi"
