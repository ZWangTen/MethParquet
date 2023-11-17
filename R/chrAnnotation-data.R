#' Sample CpG Annotation Dataset
#'
#' A subset of CpG annotation data from the Illumina.
#'
#' @docType data
#'
#' @usage data(chrAnnotation)
#'
#' @format ## `chrAnnotation`
#' A data frame with 240 rows and 52 columns:
#' \describe{
#'   \item{IlmnID, Name}{CpG names}
#'   \item{CHR}{Chromosome}
#'   \item{MAPINFO}{Position for CpG in chromosome}
#'   \item{UCSC_RefGene_Name}{Gene names}
#'   ...
#' }
#' @source <https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html>
"chrAnnotation"
