% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ewas_bmi-data.R
\docType{data}
\name{ewas_bmi}
\alias{ewas_bmi}
\title{A mock Effect data for BMI phenotype to calculate Methylation Risk Score}
\format{
## `ewas_bmi`
A data frame with 477 rows and 7 columns:
\describe{
  \item{CpG}{CpG names}
  \item{Effect}{Effect estimate associated with BMI}
  \item{chr}{Chromosome}
  ...
}
}
\source{
<https://www.ewascatalog.org/?trait=bmi>
}
\usage{
data(ewas_bmi)
}
\description{
For demonstration purposes, some CpG names from original data are masked to match those in MethData.
The original Effect dataset for BMI is obtained from the EWAS Catalog.
}
\references{
Geurts et al. (2018) Int J Obes (Lond) 42(4):887-896
(\href{https://pubmed.ncbi.nlm.nih.gov/29278407/}{PubMed})
}
\keyword{datasets}
