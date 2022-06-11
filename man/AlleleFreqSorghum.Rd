\name{AlleleFreqSorghum}
\alias{AlleleFreqSorghum}
\title{HighandLowBulkFrequencies}
\usage{
AlleleFreqSorghum(SNPset = frame, k = 4,step = 1,Chrom1 = "Chr01",Chrom2 = "Chr02", Chrom3 = "Chr03",Chrom4 = "Chr04", Chrom5 = "Chr05",Chrom6 = "Chr06",Chrom7 = "Chr07",Chrom8 = "Chr08",Chrom9 = "Chr09",Chrom10 = "Chr10" )
}
\arguments{
\item{SNPset}{A SNP set generated from QTLseqr package after vcf has been imported with importFromVCF, filtered upstream with filterSNPs and analyzed with runQTLseqAnalysis and runGprimeAnalysis}
\item{k}{numeric vector or list. In the simplest case this is an integer specifying the window width (in numbers of observations) which is aligned to the original sample according to the align argument. Alternatively, width can be a list regarded as offsets compared to the current time, see below for details.}
\item{step}{calculate FUN at every by-th time point rather than every point. by is only used if width is length 1 and either a plain scalar or a list.}
\item{Chrom1}{Character name of Chromosome 1 i.e. "Chr01"}
\item{Chrom2}{Character name of Chromosome 2 i.e. "Chr02"}
\item{Chrom3}{Character name of Chromosome 3 i.e. "Chr03"}
\item{Chrom4}{Character name of Chromosome 4 i.e. "Chr04"}
\item{Chrom5}{Character name of Chromosome 5 i.e. "Chr05"}
\item{Chrom6}{Character name of Chromosome 6 i.e. "Chr06"}
\item{Chrom7}{Character name of Chromosome 7 i.e. "Chr07"}
\item{Chrom8}{Character name of Chromosome 8 i.e. "Chr08"}
\item{Chrom9}{Character name of Chromosome 9 i.e. "Chr09"}
\item{Chrom10}{Character name of Chromosome 10 i.e. "Chr10"}
}
\description{
Plots A line graph with points for High and Low Bulk Allelic Frequencies for model organism Sorghum
}

