# QTLAlleleFrequenciesPlot
*Sorghum*

```r

install.packages("vcfR")
library(vcfR)

install.packages("ggplot2")
library(ggplot2)

install.packages("tidyr")
library(tidyr)

install.packages("dplyr")
library(dplyr)

install.packages("zoo")
library(zoo)

install.packages("devtools")
library(devtools)

devtools::install_github("PBGLMichaelHall/QTLseqr",force=TRUE)
library(QTLseqr)

devtools::install_github("PBGLMichaelHall/QTLAlleleFrequenciesPlot")
library(QTLAlleleFrequenciesPlot)

#Set High bulk and Low bulk sample names and parser generated file name
#The file name is generated from the QTLParser_1_MH function in line 119


#Choose which chromosomes/contigs will be included in the analysis,

Chroms <- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")

importFromVCF(file = "freebayes_D2.filtered.vcf",highBulk = "D2_F2_tt", 
lowBulk = "D2_F2_TT",chromList = Chroms,filename = "HallSorghum")

HighBulk <- "D2_F2_tt"
LowBulk <- "D2_F2_TT"
file <- "HallSorghum.csv"

#Choose which chromosomes/contigs will be included in the analysis,

df <-
  importFromTable(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  )

head(df)

ggplot(data =df) + geom_histogram(aes(x = DP.LOW + DP.HIGH)) + xlim(0,400)

ggplot(data = df) + geom_histogram(aes(x = REF_FRQ))


df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.10,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    depthDifference = 100,
    minSampleDepth = 30,
    minGQ = 99
  )
  


df_filt <- runQTLseqAnalysis(df_filt,
                             windowSize = 5000000,
                             popStruc = "F2",
                             bulkSize = c(45, 38),
                             replications = 10000,
                             intervals = c(95, 99)
)




df_filt <- runGprimeAnalysis(df_filt,
                             windowSize = 5000000,
                             outlierFilter = "deltaSNP",
                             filterThreshold = 0.3)
                             
QTLAlleleFreq::AlleleFreqSorghum(SNPset = df_filt, k = 4,step = 1,
Chrom1 = "Chr01",Chrom2 = "Chr02", Chrom3 = "Chr03",Chrom4 = "Chr04",
Chrom5 = "Chr05",Chrom6 = "Chr06",Chrom7 = "Chr07",Chrom8 = "Chr08"
,Chrom9 = "Chr09",Chrom10 = "Chr10")
```



![Rplot01](https://user-images.githubusercontent.com/93121277/173185271-ae15afbd-6787-49d0-80ed-93e91ef6a685.png)


```r

QTLAlleleFreq::AlleleFreqAny(SNPset = df_filt,k = 8, step = 1,Chrom = "Chr04", scalar = 1, LowerBound = 0 ,UpperBound =90000000)

```


![Rplot04](https://user-images.githubusercontent.com/93121277/173185533-f9f11bab-69ae-45b2-99ba-bea1535495bc.png)

```r

VCFstat::FacetChromQual(vcf = "freebayes_D2.filtered.vcf", chromlist = chromlist,windowSize = 5000000, ncol=10)

```




*Brachypodium Model Organism a grass species*
```r

QTLAlleleFreq::NumOfSNPs(file="freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf"
,chromList = c("NC_016131.3","NC_016132.3","NC_016133.3","NC_016134.3","NC_016135.3"),
filename = "HallBrach",filter = FALSE, windowsize = 1e6,var = "nSNPs",
scaleChroms = TRUE, line = TRUE, WriteTable = TRUE)

```

![Rplot01](https://user-images.githubusercontent.com/93121277/173334622-ae289522-93f7-4d79-b87a-c66803af0837.png)




