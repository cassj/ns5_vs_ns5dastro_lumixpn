#!/usr/bin/Rscript

library(beadarray)

# because I can't be bothered to type all the quotes
qw <- function(...) {
  as.character(sys.call()[-1])
}

options(stringsAsFactors = FALSE);
options(scipen=10)

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]
outfile = args[2]


#filename<-"/mnt/data/1738277028_Sample_Probe_Profile_no_normalisation.txt"
#outfile <- "/mnt/data/limma_results.csv"

test <- read.csv(filename, skip=7)

BSData<-readBeadSummaryData(filename,
                            sep="\t",
                            skip=7,
                            ProbeID="ProbeID",
                            columns=list(
                              exprs="AVG_Signal",
                              se.exprs="BEAD_STDEV",
                              NoBeads="Avg_NBEADS",
                              Detection="Detection"
		   )
)	

E = normaliseIllumina(BSData, method="quantile", transform="log2")
data <- exprs(E)

library(limma)

cols <- rep(qw(ns, astro),4)
ns<-which(cols=="ns")
astro<-which(cols=="astro")

design<-matrix(0,nrow=(ncol(data)), ncol=2)
colnames(design)<-c("ns","astro")
design[ns,"ns"]<-1
design[astro,"astro"]<-1


fit<-lmFit(data, design)
cont.matrix<-makeContrasts(nsvsastro=astro-ns, levels=design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)

write.fit(ebFit, file=outfile , adjust="BH")
data<-read.table(outfile, sep="\t", header=T)

data<- topTable(ebFit, number=nrow(data))
write.csv(data,outfile)
















