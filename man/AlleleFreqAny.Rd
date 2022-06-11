\name{AlleleFreqAny}
\alias{AlleleFreqAny}
\title{AlleleFreqAnyModelOrganism}
\usage{
  AlleleFreqAny(SNPset, Chrom, k, step, scalar}
\arguments{
  \item{SNPset}{A SNP set generated from QTLseqr package after vcf has been imported with importFromVCF, filtered upstream with filterSNPs and analyzed with runQTLseqAnalysis and runGprimeAnalysis}
  \item{k}{numeric vector or list. In the simplest case this is an integer specifying the window width (in numbers of observations) which is aligned to the original sample according to the align argument. Alternatively, width can be a list regarded as offsets compared to the current time, see below for details.}
  \item{step}{calculate FUN at every by-th time point rather than every point. by is only used if width is length 1 and either a plain scalar or a list.}
  \item{Chrom}{Character name of any Chromosome/Contig 1 i.e. "Chr01"}
  \item{scalar}{A scalar quanity used in multiplication against number of SNPs in sliding window. If Points are too lare scale down, if points are too small scale up}
  \item{LowerBound}{Lower Bound on newly averaged SNP positions}
  \item{UpperBound}{Upper Bound on newly averaged SNP positions}
}
\description{
  Plots A line graph with points with respect to Bulk segregant analysis of High Bulk and Low Bulk Allelic Frequencies with respect to a model organism
}

