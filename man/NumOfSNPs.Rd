\name{NumOfSNPs}
\alias{NumOfSNPs}
\title{NumOfSNPs}
\usage{
  NumOfSNPs(file = "example.vcf", chromList = c("NC_016131.3","NC_016132.3","NC_016133.3","NC_016134.3","NC_016135.3"), filename = "Model_organism_SNPs", filter = FALSE, windowsize = 1e6, var = "nSNPs", scaleChroms = TRUE, line = TRUE, WriteTable = TRUE}
\arguments{
  \item{A VCF File}{Provide a Variant Call File to Analyze}
  \item{chromList}{Provide a character vector specifying chromosome names i.e. chromList = c("NC_016131.3","NC_016132.3","NC_016133.3","NC_016134.3","NC_016135.3")}
  \item{filename}{Provide prefix to csv file name generated when WriteTable = TRUE}
  \item{filter}{Boolean Value to filter variant file by Filter Info field or not}
  \item{windowsize}{Specify window size to analyze variants}
  \item{var}{Specify which variable to plot, either "nSNPs","DP",or "QUAL"}
  \item{scaleChroms}{Scale Chromsomes or not}
  \item{line}{Plot a line or point graph}
  \item{WriteTable}{Boolean, TRUE or FALSE to write a csv file or not}
}
\description{
  Plots A line graph of number of SNPs called for each chromosome
}