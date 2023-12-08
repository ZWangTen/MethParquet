
#' Plot methylation against phenotype
#' @description Create a plot of methylation against phenotype for each CpG site
#'
#' @param db_obj MethList created using \code{\link{create_methlist}}
#' @param sites Character string for CpG names
#' @param trait Exposure phenotype in `subject_annot` from `MethList`.
#' @param output_plot File path to save the plot. Default: FALSE.
#' @param width Width in inch for the saved output plot when given path (`output_plot`). Default: 9
#' @param height Height in inch for the saved output plot when given path (`output_plot`). Default: 8
#'
#' @return A plot with methylation beta values on y axis against the phenotype on x axis
#' @export
#' @import ggplot2
#' @import tidyverse
#' @importFrom arrow open_dataset
#'
#' @examples
#' library(tidyverse)
#' library(arrow)
#' library(ggplot2)
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
#' # Run linear model and get top results
#' ewas_lm <- lm_ewas_outcome(db_obj=mlist,trait='bmi',select_sites='full',select_chr=FALSE,
#' covariates_string = c('age','sex'),out_position='MAPINFO',block_size=100000)
#' sites<-ewas_lm$CpG[ewas_lm$fdr_bh<0.05][1:3]
#'
#' # Plot top results
#' plot_MethPheno(mlist,sites=sites,trait='bmi',output_plot = FALSE)
#' unlink(path,recursive=TRUE)

plot_MethPheno <- function(db_obj,sites,trait,output_plot=F,width=9,height=8){
  pheno <- db_obj$subject_annot
  meth.dat <- db_obj$db %>% filter(CpG %in% sites) %>% as.data.frame() %>%
    dplyr::select(-CHR) %>% dplyr::select(c(CpG,pheno[,"sample_id"]))

  meth.for.plot <- gather(meth.dat, key = "Person", value = "Beta_value", as.character(pheno[,"sample_id"]))
  meth.for.plot$phenotype <- pheno[,trait][match(meth.for.plot[, "Person"], pheno[,"sample_id"])]

  p <- ggplot(meth.for.plot, aes(Beta_value, phenotype))
  g <- p + geom_point( ) + geom_smooth(method=lm, se=FALSE, fullrange=FALSE) +
    facet_wrap(~CpG,scales = "free") + theme_bw() +  labs(y= paste(trait))

  if (!isFALSE(output_plot)){
    ggsave(paste0(output_plot, "/Beta_versus_phenotype.png"), g, width = width, height = height, dpi = 300, units = "in", device='png')
  }
  return(g)
}
