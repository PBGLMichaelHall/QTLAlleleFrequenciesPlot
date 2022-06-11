#' @title AlleleFreqSorghum
#' @description Plots A line graph with points for High and Low Bulk Allelic Frequencies for model organism Sorghum
#' @param SNPset A SNP set generated from QTLseqr package after vcf has been imported with importFromVCF, filtered upstream with filterSNPs and analyzed with runQTLseqAnalysis and runGprimeAnalysis
#' @param k numeric vector or list. In the simplest case this is an integer specifying the window width (in numbers of observations) which is aligned to the original sample according to the align argument. Alternatively, width can be a list regarded as offsets compared to the current time, see below for details.
#' @param step calculate FUN at every by-th time point rather than every point. by is only used if width is length 1 and either a plain scalar or a list.
#' @param Chrom1 Character name of Chromosome 1 i.e. "Chr01"
#' @param Chrom2 Character name of Chromosome 1 i.e. "Chr02"
#' @param Chrom3 Character name of Chromosome 1 i.e. "Chr03"
#' @param Chrom4 Character name of Chromosome 1 i.e. "Chr04"
#' @param Chrom5 Character name of Chromosome 1 i.e. "Chr05"
#' @param Chrom6 Character name of Chromosome 1 i.e. "Chr06"
#' @param Chrom7 Character name of Chromosome 1 i.e. "Chr07"
#' @param Chrom8 Character name of Chromosome 1 i.e. "Chr08"
#' @param Chrom9 Character name of Chromosome 1 i.e. "Chr09"
#' @param Chrom10 Character name of Chromosome 1 i.e. "Chr10"
#' @export AlleleFreqSorghum


AlleleFreqSorghum<- function(SNPset, Chrom1,Chrom2,Chrom3,Chrom4,Chrom5,Chrom6,Chrom7,Chrom8,Chrom9,Chrom10, k, step){
  frame <- SNPset %>% dplyr::mutate(LowRef = AD_REF.LOW, HighRef = AD_REF.HIGH, LowAlt = AD_ALT.LOW, HighAlt = AD_ALT.HIGH) %>% dplyr::select(CHROM,POS,DP.LOW,DP.HIGH,LowRef, HighRef, LowAlt, HighAlt,nSNPs)
  frame1 <- frame %>% dplyr::filter(CHROM == Chrom1)

  frame1$NewPOS <- rollapply(frame1$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame1$LowBulkFreq <- frame1$LowAlt/frame1$DP.LOW
  frame1$newLowBulkFreq <- rollapply(frame1$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame1$NewPOS <- rollapply(frame1$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame1$highBulkFreq <- frame1$HighAlt/frame1$DP.HIGH
  frame1$newhighBulkFreq <- rollapply(frame1$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame2 <- frame %>% dplyr::filter(CHROM == Chrom2)

  frame2$NewPOS <- rollapply(frame2$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame2$LowBulkFreq <- frame2$LowAlt/frame2$DP.LOW
  frame2$newLowBulkFreq <- rollapply(frame2$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame2$NewPOS <- rollapply(frame2$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame2$highBulkFreq <- frame2$HighAlt/frame2$DP.HIGH
  frame2$newhighBulkFreq <- rollapply(frame2$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame3 <- frame %>% dplyr::filter(CHROM == Chrom3)

  frame3$NewPOS <- rollapply(frame3$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame3$LowBulkFreq <- frame3$LowAlt/frame3$DP.LOW
  frame3$newLowBulkFreq <- rollapply(frame3$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame3$NewPOS <- rollapply(frame3$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame3$highBulkFreq <- frame3$HighAlt/frame3$DP.HIGH
  frame3$newhighBulkFreq <- rollapply(frame3$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame4 <- frame %>% dplyr::filter(CHROM == Chrom4)

  frame4$NewPOS <- rollapply(frame4$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame4$LowBulkFreq <- frame4$LowAlt/frame4$DP.LOW
  frame4$newLowBulkFreq <- rollapply(frame4$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame4$NewPOS <- rollapply(frame4$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame4$highBulkFreq <- frame4$HighAlt/frame4$DP.HIGH
  frame4$newhighBulkFreq <- rollapply(frame4$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame5 <- frame %>% dplyr::filter(CHROM == Chrom5)

  frame5$NewPOS <- rollapply(frame5$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame5$LowBulkFreq <- frame5$LowAlt/frame5$DP.LOW
  frame5$newLowBulkFreq <- rollapply(frame5$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame5$NewPOS <- rollapply(frame5$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame5$highBulkFreq <- frame5$HighAlt/frame5$DP.HIGH
  frame5$newhighBulkFreq <- rollapply(frame5$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame6 <- frame %>% dplyr::filter(CHROM == Chrom6)

  frame6$NewPOS <- rollapply(frame6$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame6$LowBulkFreq <- frame6$LowAlt/frame6$DP.LOW
  frame6$newLowBulkFreq <- rollapply(frame6$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame6$NewPOS <- rollapply(frame6$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame6$highBulkFreq <- frame6$HighAlt/frame6$DP.HIGH
  frame6$newhighBulkFreq <- rollapply(frame6$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame7 <- frame %>% dplyr::filter(CHROM == Chrom7)

  frame7$NewPOS <- rollapply(frame7$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame7$LowBulkFreq <- frame7$LowAlt/frame7$DP.LOW
  frame7$newLowBulkFreq <- rollapply(frame7$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame7$NewPOS <- rollapply(frame7$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame7$highBulkFreq <- frame7$HighAlt/frame7$DP.HIGH
  frame7$newhighBulkFreq <- rollapply(frame7$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame8 <- frame %>% dplyr::filter(CHROM == Chrom8)

  frame8$NewPOS <- rollapply(frame8$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame8$LowBulkFreq <- frame8$LowAlt/frame8$DP.LOW
  frame8$newLowBulkFreq <- rollapply(frame8$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame8$NewPOS <- rollapply(frame8$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame8$highBulkFreq <- frame8$HighAlt/frame8$DP.HIGH
  frame8$newhighBulkFreq <- rollapply(frame8$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame9 <- frame %>% dplyr::filter(CHROM == Chrom9)

  frame9$NewPOS <- rollapply(frame9$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame9$LowBulkFreq <- frame9$LowAlt/frame9$DP.LOW
  frame9$newLowBulkFreq <- rollapply(frame9$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame9$NewPOS <- rollapply(frame9$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame9$highBulkFreq <- frame9$HighAlt/frame9$DP.HIGH
  frame9$newhighBulkFreq <- rollapply(frame9$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame10 <- frame %>% dplyr::filter(CHROM == Chrom10)

  frame10$NewPOS <- rollapply(frame10$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame10$LowBulkFreq <- frame10$LowAlt/frame10$DP.LOW
  frame10$newLowBulkFreq <- rollapply(frame10$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame10$NewPOS <- rollapply(frame10$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame10$highBulkFreq <- frame10$HighAlt/frame10$DP.HIGH
  frame10$newhighBulkFreq <- rollapply(frame10$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  data <- rbind(frame1, frame2, frame3, frame4, frame5, frame6, frame7, frame8, frame9, frame10)
  str(data)
  e <- ggplot(data, aes(x = NewPOS)) + geom_line(aes(y = newLowBulkFreq), color = "orange") + geom_point(aes(y = newLowBulkFreq), color = "orange") + geom_line(aes(y = highBulkFreq), color = "black") + geom_point(aes(y = highBulkFreq), color = "black") + facet_grid(rows = vars(CHROM))
  print(e)
}

#' @title AlleleFreqAny
#' @description Plots A line graph with points for High and Low Bulk Allelic Frequencies for model organism Sorghum
#' @param SNPset A SNP set generated from QTLseqr package after vcf has been imported with importFromVCF, filtered upstream with filterSNPs and analyzed with runQTLseqAnalysis and runGprimeAnalysis
#' @param k numeric vector or list. In the simplest case this is an integer specifying the window width (in numbers of observations) which is aligned to the original sample according to the align argument. Alternatively, width can be a list regarded as offsets compared to the current time, see below for details.
#' @param step calculate FUN at every by-th time point rather than every point. by is only used if width is length 1 and either a plain scalar or a list.
#' @param Chrom Character name of Chromosome 1 i.e. "Chr01"
#' @param scalar A scalar quantity used in multiplication against nSNPs. If Points are too large scale down, if Points are too small scale up
#' @param LowerBound Lower Bound on newly averaged SNP positions
#' @param UpperBound Upper Bound on newly averaged SNP positions
#' @export AlleleFreqAny


AlleleFreqAny <- function(SNPset, Chrom, k, step, scalar, LowerBound, UpperBound){
  frame <- SNPset %>% dplyr::mutate(LowRef = AD_REF.LOW, HighRef = AD_REF.HIGH, LowAlt = AD_ALT.LOW, HighAlt = AD_ALT.HIGH) %>% dplyr::select(CHROM,POS,DP.LOW,DP.HIGH,LowRef, HighRef, LowAlt, HighAlt,nSNPs)
  frame1 <- frame %>% dplyr::filter(CHROM == Chrom)
  str(frame1)
  frame1$NewPOS <- rollapply(frame1$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame1$LowBulkFreq <- frame1$LowAlt/frame1$DP.LOW
  frame1$newLowBulkFreq <- rollapply(frame1$LowBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame1$NewPOS <- rollapply(frame1$POS, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)
  frame1$highBulkFreq <- frame1$HighAlt/frame1$DP.HIGH
  frame1$newhighBulkFreq <- rollapply(frame1$highBulkFreq, 2*k-1, function(x) max(rollmean(x,k)), by = step, partial = TRUE)

  frame1 <- frame1 %>% filter(LowerBound <= NewPOS & NewPOS <= UpperBound)


  e <- ggplot(frame1, aes(x = NewPOS)) + geom_line(aes(y = newLowBulkFreq), color = "orange") + geom_point(aes(y = newLowBulkFreq), color = "orange", size = scalar*frame1$nSNPs) + geom_line(aes(y = highBulkFreq), color = "black") + geom_point(aes(y = highBulkFreq), color = "black", size = scalar*frame1$nSNPs)
  print(e)
}

