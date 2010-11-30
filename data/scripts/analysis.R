wd<-setwd("/home/cassj/work/astrocyte_array_study/manu_data/raw")

### 
# Load the data


library(beadarray)

#We only have bead-summary level data:

dataFile<-"v3/raw_data.csv"
qcFile <- "v3/qc_info.csv"

BSData<-readBeadSummaryData(dataFile, qcFile=qcFile,
		   sep=",",
                   skip=1,
                   columns=list(
		       exprs="AVG_Signal",
		       se.exprs="BEAD_STDERR",
		       NoBeads="Avg_NBEADS",
		       Detection="Detection"
		   ),
                   qc.skip=1,
		   qc.sep=",",
		   qc.columns=list(
		       exprs="AVG_Signal",	
		       se.exprs="BEAD_STDERR",
                       NoBeads="Avg_NBEADS", 
                       Detection="Detection Pval"
		   ))	




# this is a BSData object, 
# which isa ExpressionSetIllumina, 
# which isa ExpressionSet

# It doesn't like the sample_sheet.csv file, 
# so need to defined the phenoData ourselves.

pd <- data.frame(cell.line=rep(c("NS","Astro"),4))
rownames(pd) = c("NS5.1","Astro.1","NS5.2", "Astro.2", "NS5.3", "Astro.3","NS5.4", "Astro.4")
metadata<-data.frame(labelDescription=c("Cell Line (Astrocyte or NS5"))
pd <- new("AnnotatedDataFrame", data = pd, varMetadata=metadata)

phenoData(BSData)<-pd

save(BSData, file="BSData.RData")

# Need to add phenoData to the qc object, but not sure how to get at it
# According to the Vignette, the rownames of the QC exprs matrix should
# be the target names, but you can't have duplicate rownames in a dataframe
# Perhaps they should go in the featureData slot.

qc.names<-as.character(read.csv("v3/qc_info.csv", skip=1)[,1])

                            
###########

load("BSData.RData")

 # > BSData
 # ExpressionSetIllumina (storageMode: list)
 # assayData: 24613 features, 8 samples 
 #   element names: exprs, se.exprs, NoBeads, Detection, Narrays, arrayStDev, DiffScore 
 # phenoData
 #   rowNames: NS5.1, Astro.1, ..., Astro.4  (8 total)
 #   varLabels and varMetadata description:
 #     cell.line: Cell Line (Astrocyte or NS5
 # featureData
 #   featureNames: 
 #   fvarLabels and fvarMetadata description: none
 # experimentData: use 'experimentData(object)'
 # Annotation:  
 # QC Information
 #  Available Slots:  exprs se.exprs Detection NoBeads controlType
 #   featureNames: 1, 2, ..., 857, 858
 #   sampleNames: X1738277028_A.AVG_Signal, X1738277028_B.AVG_Signal, ..., X1738277028_G.AVG_Signal, X1738277028_H.AVG_Signal
 
 # > dim(assayData(BSData)$exprs)
 # [1] 24613     8
 # > dim(assayData(BSData)$se.exprs)
 # [1] 24613     8

 # > pData(BSData)
 #         cell.line
 # NS5.1          NS
 # Astro.1     Astro
 # NS5.2          NS
 # Astro.2     Astro
 # NS5.3          NS
 # Astro.3     Astro
 # NS5.4          NS
 # Astro.4     Astro




###### Normalisation and Quality Control


## how many non-positive values are we dealing with?
#
#
#for (i in colnames(exprs(BSData))){
#  this.col<-exprs(BSData)[,i]
#  negs<-length(which(this.col<=0))
#  cat(i, negs, "of", length(this.col), "\n", sep=" ")
#}
#
# # X1.AVG_Signal 7489 of 24613 
# # X3.AVG_Signal 7307 of 24613 
# # X5.AVG_Signal 7553 of 24613 
# # X7.AVG_Signal 7282 of 24613 
# # X9.AVG_Signal 7555 of 24613 
# # X11.AVG_Signal 7428 of 24613 
# # X13.AVG_Signal 7505 of 24613 
# # X15.AVG_Signal 7354 of 24613 
#
#
## This is actually a pretty large amount. 
#
#ok<-logical(nrow(exprs(BSData)))
#for(i in 1:nrow(exprs(BSData))) {
#  ok[i]<-all(exprs(BSData)[i,]>0)
#}
#
#mat<-exprs(BSData)[ok,]
# # > dim(mat)
# # [1] 13785     8
#
#
#mat<-log2(mat)
#
#postscript(file="v3/raw_beadsummary_exprs_boxplot.eps", paper="special", height=10, width=10)
#  boxplot(mat, las=2, main="log2 Raw Bead Summary Data")
#dev.off()
#postscript(file="v3/number_of_beads_boxplot.eps", paper="special", width=10, height=10)
#  boxplot(NoBeads(BSData), las=2, main="Number of Beads")
#dev.off() 
#
#
#ns<-c(1,3,5,7)
#astro<-c(2,4,6,8)
#
#postscript(file="v3/ns5_log2_pairs.eps", paper="special", height=10, width=10)
#  pairs(mat[,ns], pch='.')
#dev.off()
#
#postscript(file="v3/astro_log2_pairs.eps", paper="special", height=10, width=10)
#   pairs(mat[,astro], pch='.')
#dev.off()
#



###
# Filter on Detection p.values

# note that I had to rename the cols from "Detection Pval" to Detection
# as the readBeadSummaryData didn't work with "Detection Pval". 
# Maybe doesn't like the space.

# The detection score in v1 is just 1-p.val for reasons that escape me
# v3 has a detection pvalue.

det<-assayData(BSData)$Detection
 # > dim(det)
 # [1] 24613     8

#Detection of >=0.8 in 4/8 
sig.count<-apply(det, 1, function(x){length(which(x<0.2))})
detected<-sig.count[sig.count>=4]

 # > length(detected)
 # [1] 13882

# Save filtered data
detected<-names(detected)





##############
# Kee-Yew's analysis

kee.yew<-BSData

#set values <10 to 10

for (i in colnames(exprs(kee.yew))){
  exprs(kee.yew)[which(exprs(kee.yew)[,i]<10),i]<-10
}

####
# normalise per chip to median 

chip.meds<-apply(exprs(kee.yew), 2, median)

 # > length(gene.meds)
 # [1] 24620

exprs(kee.yew)<-exprs(kee.yew)/rep(chip.meds, each=nrow(exprs(kee.yew)))


###
# normalise per gene to median

gene.meds<-apply(exprs(kee.yew), 1, median)

 # > length(chip.meds)
 # [1] 8

exprs(kee.yew)<-exprs(kee.yew)/gene.meds
save(kee.yew, file="v3/kee-yew_normalised_BSData.RData")



###
#de: just a 2-sided t-test with Welch's estimate 

library(stats)

#ns5=1,3,5,7
#astro=2,4,6,7

#ns5=y, astro=x => +ve tstats indicate higher expression in the astro

tests<-apply(exprs(kee.yew),1,function(x){ t.test(x[c(2,4,6,8)],x[c(1,3,5,7)] )  })
save(tests, file="v3/tests.RData")

t.stats<-unlist(lapply(tests,function(x){x$statistic}))
names(t.stats)<-rownames(exprs(kee.yew))
t.stats<-sort(t.stats, decreasing=TRUE)
save(t.stats, file="v3/ky.tstats.RData")

p.vals<-unlist(lapply(tests,function(x){x$p.value}))
names(p.vals)<-rownames(exprs(kee.yew))
p.vals<-sort(p.vals)
save(p.vals, file="v3/pvals.RData")
kee.yew.ord<-names(p.vals)

det.pvals<-p.vals[detected]
det.pvals<-sort(det.pvals)

det.tstats<-t.stats[detected]
det.tstats<-sort(det.tstats, decreasing=TRUE)

kee.yew.ord.detected<-names(det.pvals)


###
# Benjamini-Hochberg FDRs

library(multtest)

bh.fdr<-mt.rawp2adjp(p.vals, "BH")[[1]]
rownames(bh.fdr)<-names(p.vals)
save(bh.fdr, file="v3/keeyew.bh.fdr.RData")

 # > colnames(bh.fdr)
 # [1] "rawp" "BH"  

keeyew.fdr.05<-rownames(bh.fdr[which(bh.fdr[,"BH"]<0.05),])

 # > length(keeyew.fdr.05)
 # [1] 5231

write.csv(bh.fdr, file="v3/keeyew.bh.fdr.csv")

## And with detection p-val filtering:

bh.fdr.detected<-mt.rawp2adjp(det.pvals, "BH")[[1]]
rownames(bh.fdr.detected)<-names(det.pvals)
save(bh.fdr.detected, file="v3/keeyew.bh.fdr.detected.RData")

det.fdr<-bh.fdr.detected[,"BH"]

keeyew.fdr.05.detected<-rownames(bh.fdr.detected[which(bh.fdr.detected[,"BH"]<0.05),])

 # >length(ind)
 # [1] 6172




#########
## alternative: recommended by beadarray package

library(affy)
library(lumi)

#Can't log with -ve numbers. 
#set values <10 to 10

biocbead<-BSData
for (i in colnames(exprs(biocbead))){
  exprs(biocbead)[which(exprs(biocbead)[,i]<10),i]<-10
}

biocbead<-normaliseIllumina(biocbead,
                       method="quantile",
                       transform="log2"
)


#limma ebayes mod t-stat thing

astro<-rep(c(0,1),4)
ns5<-rep(c(1,0),4)
design<-cbind(astro, ns5)
rownames(design)<-rownames(pData(phenoData(BSData)))

 #         astro ns5
 # NS5.1       0   1
 # Astro.1     1   0
 # NS5.2       0   1
 # Astro.2     1   0
 # NS5.3       0   1
 # Astro.3     1   0
 # NS5.4       0   1
 # Astro.4     1   0

fit<-lmFit(exprs(biocbead), design)


# -ve fold change implies Astro less than ns5.

cont.matrix<-makeContrasts(AstroVNS5=astro-ns5, levels=design)

 # > cont.matrix
 #        Contrasts
 # Levels  AstroVNS5
 #   astro         1
 #   ns5          -1

fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)
limma.genelist<-topTable(ebFit, number=nrow(exprs(BSData)),sort.by="p" )

 # > limma.genelist[1:5,]
 #            ID     logFC  AveExpr         t      P.Value    adj.P.Val        B
 # 14110 3390156 -6.464009 6.553933 -195.1274 2.719074e-13 6.692458e-09 18.70561
 # 6864   610112  6.984955 6.814406  163.1581 8.516507e-13 9.207046e-09 18.31352
 # 4607  1740292 -9.304864 7.974360 -156.2536 1.122217e-12 9.207046e-09 18.20349
 # 13996 5570068 -5.650043 6.146949 -138.0265 2.476054e-12 1.030111e-08 17.85237
 # 16710  430079 -5.670105 6.156981 -137.7510 2.507824e-12 1.030111e-08 17.84627
 
limma.ord<-limma.genelist$ID

limma.fdr.05<-limma.genelist$ID[limma.genelist$adj <= 0.05]

limma.pvals<-limma.genelist$P.Value
names(limma.pvals)<-limma.genelist$ID

limma.tstats<-limma.genelist$t
names(limma.tstats)<-limma.genelist$ID
limma.tstats<-sort(limma.tstats, decreasing=TRUE)

 # > length(limma.fdr.05)
 # [1] 6844

write.fit(ebFit, file="v3/limma.fit.csv", adjust="BH", sep=",")

# same for detection pval filtered
fit.detected<-lmFit(exprs(biocbead)[detected,], design)
fit.detected<-contrasts.fit(fit.detected, cont.matrix)
ebFit.detected<-eBayes(fit.detected)
limma.genelist.detected<-topTable(ebFit.detected, number=nrow(exprs(biocbead)[detected,]),sort.by="p" )
limma.ord.detected<-limma.genelist.detected$ID
limma.fdr.05.detected<-limma.genelist.detected$ID[limma.genelist.detected$adj <= 0.05]
write.fit(ebFit.detected, file="v3/limma.fit.detected.csv", adjust="BH", sep=",")

limma.det.pvals<-limma.genelist.detected$P.Value
names(limma.det.pvals)<-limma.genelist.detected$ID
limma.det.tstats<-limma.genelist.detected$t
names(limma.det.tstats)<-limma.genelist.detected$ID
limma.det.tstats<-sort(limma.det.tstats, decreasing=TRUE)
limma.det.fdr<-limma.genelist.detected$adj.P.Val
names(limma.det.fdr)<-limma.genelist.detected$ID





###
# Check to make sure this isn't massively different from 
# the original genelist.

ky.orig<-read.csv("../NS5_vs_Astro.csv")
ky.orig.p05<-read.csv("../NS5_vs_Astro_p05.csv")

save(ky.orig, file="v3/ky_orig.RData")
save(ky.orig.p05, file="v3/ky_orig_p05.RData")

 # > colnames(ky.orig)
 #  [1] "Systematic.Name"                 "Synonyms"                       
 #  [3] "Genbank"                         "Symbol"                         
 #  [5] "Definition"                      "Transcript"                     
 #  [7] "Cell.Line.Astro..raw"            "Cell.Line.Astro..t.test.p.value"
 #  [9] "Cell.Line.NS5..raw"              "Cell.Line.NS5..t.test.p.value"  


 # > colnames(ky.orig.p05)
 #  [1] "Gene.Name"             "Fold.Change"           "Common"               
 # [4] "Synonyms"              "Genbank"               "Symbol"               
 #  [7] "Definition"            "Transcript"            "GO.Molecular.Function"
 # [10] "GO.Cellular.Component" "GO.Biological.Process"


 #Number of sig at p05:

 # > nrow(ky.orig.p05)
 # [1] 6463
 # > length(keeyew.fdr.05)
 # [1] 5231
 # > length(keeyew.fdr.05.detected)
 # [1] 6172
 # > length(limma.fdr.05)
 # [1] 6844
 # > length(limma.fdr.05.detected)
 # [1] 7279

v1<-read.csv("NB_NS5_Astrocytes_gene_profile.csv", skip=7)
annot<-read.delim("Illumina_MouseRef-8_V1_1_R1_11234312_A.csv")

 # > dim(v1)
 # [1] 24611    65
 # > dim(annot)
 # [1] 24613    27
 # > dim(exprs(BSData))
 # [1] 24613     8

#discrepancy is due to the existance of control probes in the v1 file.

#can't use the TargetID as input to the readBeadSummaryData function
#for some reason, so we need to do mapping 

v3<-read.csv("v3/full_raw_data.csv", skip=1)
rownames(v3)<-v3$ProbeID   #this is AAID so we can subset

targetids<-as.character(v3$TargetID)
names(targetids)<-as.character(v3$ProbeID)
 # > length(targetids)
 # [1] 24613
 # > length(unique(targetids))
 # [1] 24607

kee.yew.orig.targets<-as.character(ky.orig.p05[,"Gene.Name"])

kee.yew.targets<-targetids[kee.yew.ord]
kee.yew.p05.targets<-targetids[keeyew.fdr.05]
kee.yew.p05.detected.targets<-targetids[keeyew.fdr.05.detected]

limma.targets<-targetids[limma.ord]
limma.p05.targets<-targetids[limma.fdr.05]
limma.p05.detected.targets<-targetids[limma.fdr.05.detected]

#ok, how many of each are in ky's sig set (6463)


 # > length(intersect(kee.yew.orig.targets, kee.yew.p05.targets))
 # [1] 5144 (of 5231)
 # > length(intersect(kee.yew.orig.targets, kee.yew.p05.detected.targets))
 # [1] 5933 (of 6172)
 # > length(intersect(kee.yew.orig.targets, limma.p05.targets))
 # [1] 6095 (of 6844)
 # > length(intersect(kee.yew.orig.targets, limma.p05.detected.targets))
 # [1] 6157 (of 7279)


#just compare to random (on chip order) set: 

length(intersect(ky.orig.p05[,1], targetids[rownames(exprs(BSData))[1:5000]]))
 # [1] 1111 (of 5000)

#What about their intersect with each other
length(intersect(limma.p05.targets, kee.yew.p05.targets))
 # 5063 

length(intersect(limma.p05.detected.targets, kee.yew.p05.detected.targets))
 # 5885


# looks pretty good








####
# Save for GSea

#txt expression file:
gsea_mat<-(exprs(BSData))

desc<-as.character(v3$DEFINITION)
names(desc)<-rownames(v3)
desc<-desc[rownames(gsea_mat)]

tgs<-as.character(v3$TargetID)
names(tgs)<-rownames(v3)
missing<-names(which(desc==""))
desc[missing]<-tgs[missing]


gsea_mat<-cbind(rownames(gsea_mat), desc,gsea_mat)

colnames(gsea_mat)<-c("Name","Description","NS1","Astro1","NS2","Astro2","NS3","Astro3", "NS4","Astro4")

write.table(gsea_mat, file="v3/gsea_exprs_mat.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


#rnk files for limma and ky
write.table(limma.tstats, file="v3/gsea_limma_ordered.rnk", col.names=FALSE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(t.stats, file="v3/gsea_ky_ordered.rnk", col.names=FALSE, row.names=TRUE, quote=FALSE, sep="\t")







###
# OK, now export data in a variety of formats to csv files for other programs

#Using the GEO platform file, we need to map to the probe Ids we have,
#but I'm not sure which col corresponds

 # > as.character(annot$Search_Key)[1]
 # [1] "scl54150.7.1_69"
 # > rownames(bh.fdr)[1]
 # [1] "6520008"
 # > as.character(annot$Search_Key)[1]
 # [1] "scl54150.7.1_69"
 # > as.character(annot$Transcript)[1]
 # [1] "ILMN_216726"
 # > as.character(annot$ILMN_Gene)[1]
 # [1] "RENBP"
 # > as.character(annot$Source_Reference_ID)[1]
 # [1] "NM_023132.2"
 # > as.character(annot$Accession)[1]
 # [1] "NM_023132.2"
 # > as.character(annot$Probe_Id)[1]
 # [1] "ILMN_2674533"
 # > as.character(annot$Obsolete_Probe_Id)[1]
 # [1] "scl54150.7.1_69-S-10"
 # > as.character(annot$Array_Address_Id)[1]
 # [1] "5720332"


#could be array address or GI
 # > length(intersect(rownames(exprs(BSData)), annot$Array_Address_Id ))
 # [1] 24613


#WTF is this for?
#test<-rownames(exprs(BSData))[10]
# # > test
# # [1] "3170373"
#
#foo<-annot[annot$Array_Address_Id==test,-(c(1,2,13,14,16,17,18,19,20,21,22,23,24,25,26))]
#
#bar<-v3[v3$ProbeID==test, c(1,2,67,68,69,70,71,72,73)]


 # > foo
 #          Search_Key  Transcript ILMN_Gene Source_Reference_ID   RefSeq_ID
 # 3170373 NM_053220.1 ILMN_194800     V1RA5         NM_053220.1 NM_053220.1
 #         Unigene_ID Entrez_Gene_ID       GI   Accession Symbol Array_Address_Id
 # 3170373         NA         113847 16716526 NM_053220.1  V1ra5          3170373
 #         Obsolete_Probe_Id
 # 3170373   GI_16716526-S-1

 # > bar
 #         TargetID ProbeID    SEARCH_KEY        TARGET PROBEID         GID
 # 10 GI_16716526-S 3170373 GI_16716526-S GI_16716526-S 3170373 GI_16716526
 #     TRANSCRIPT   ACCESSION SYMBOL
 # 10 GI_16716526 NM_053220.1  V1ra5

# Accession numbers, obsolete probe Ids, gene symbols correstpond. looks good!








#######
# Map AAIDs to some annotation

rownames(annot)<-annot$Array_Address_Id


# Probe ID
pids<-as.character(annot$Probe_Id)
names(pids)<-rownames(annot)

# Target ID
tids<-as.character(annot$Transcript)
names(tids)<-rownames(annot)

# RefSeq ID
rids<-as.character(annot$RefSeq_ID)
names(rids)<-rownames(annot)

# EntrezGene ID
eids<-as.character(annot$Entrez_Gene_ID)
names(eids)<-rownames(annot)

# GI
gids<-as.character(annot$GI)
names(gids)<-rownames(annot)

# Gene Symbol
symbols<-as.character(annot$Symbol)
names(symbols)<-rownames(annot)















#####
# Create some subsets of AAIDs

#note that we already have 'detected'

# expression values, ordered by significance (except for raw, obv)

ky.xpn<-exprs(kee.yew)[rownames(bh.fdr),]
limma.xpn<-exprs(biocbead)[limma.genelist$ID,]
raw.xpn<-exprs(BSData)


# ave expression values in 2 cell lines

ky.astro<-apply(ky.xpn[,c(2,4,6,8)],1,mean)
names(ky.astro)<-rownames(ky.xpn)

ky.ns5<-apply(ky.xpn[,c(1,3,5,7)],1, mean)
names(ky.ns5)<-rownames(ky.xpn)

limma.astro<-apply(limma.xpn[,c(2,4,6,8)],1,mean)
names(limma.astro)<-rownames(limma.xpn)

limma.ns5<-apply(limma.xpn[,c(1,3,5,7)],1,mean)
names(limma.ns5)<-rownames(limma.xpn)

raw.astro<-apply(raw.xpn[,c(2,4,6,8)],1,mean)
raw.ns5<-apply(raw.xpn[,c(1,3,5,7)],1,mean)

ky.astro.ord<-names(ky.astro)
ky.ns5.ord<-names(ky.ns5)

limma.astro.ord<-names(limma.astro)
limma.ns5.ord<-names(limma.ns5)

raw.astro.ord<-names(raw.astro)
raw.ns5.ord<-names(raw.ns5)




#fold changes 
# actually, log2(fc) for limma and not even that for ky, 
# but in either case values >0 indicate astro is bigger
# Raw is actually a fc.

limma.fc<-limma.astro-limma.ns5     
ky.fc<-ky.astro-ky.ns5              
raw.fc<-raw.astro/raw.ns5


#Having got the fcs, order the expression averages

ky.astro<-sort(ky.astro, decreasing=TRUE)
ky.ns5<-sort(ky.ns5, decreasing=TRUE)
limma.astro<-sort(limma.astro, decreasing=TRUE)
limma.ns5<-sort(limma.ns5,decreasing=TRUE)
raw.astro<-sort(raw.astro, decreasing=TRUE)
raw.ns5<-sort(raw.ns5, decreasing=TRUE)


#and make top 1000, 2000 filters
ky.astro.1000<-names(ky.astro[1:1000])
ky.astro.2000<-names(ky.astro[1:2000])
ky.ns5.1000<-names(ky.ns5[1:1000])
ky.ns5.2000<-names(ky.ns5[1:1000])

limma.astro.1000<-names(limma.astro[1:1000])
limma.astro.2000<-names(limma.astro[1:2000])
limma.ns5.1000<-names(limma.ns5[1:1000])
limma.ns5.2000<-names(limma.ns5[1:2000])



#WTF?
 # > length(intersect(ky.astro.1000, limma.astro.1000))
 # [1] 89
 # > length(intersect(ky.astro.2000, limma.astro.2000))
 # [1] 312



# Directions
up.ky<-names(ky.fc[ky.fc>=0])
down.ky<-names(ky.fc[ky.fc<0])

up.limma<-names(limma.fc[limma.fc>=0])
down.limma<-names(limma.fc[limma.fc<0])

up.raw<-names(raw.fc[raw.fc>=1])
down.raw<-names(raw.fc[raw.fc<1])














######
# Save all of that information in a single file:

# Full dataset:

ky.ids<-rownames(bh.fdr)
limma.ids<-limma.genelist$ID
raw.ids<-rownames(exprs(BSData))
 
# Statistics

ky.stats<-cbind(t.stats[ky.ids], p.vals[ky.ids], bh.fdr[,"BH"] )
limma.stats<-cbind(limma.tstats[limma.ids], limma.pvals[limma.ids], limma.genelist$adj.P.Val)

colnames(ky.stats)<-c("t", "p", "FDR")
colnames(limma.stats)<-c("mod-t","p", "FDR")


# Annotation

col.names<-c("ArrayAddressID","ProbeID", "TargetID", "RefSeq","EntrezGene","GI","GeneSymbol", "Detection Pvalue")

ky.annot<-cbind(ky.ids,pids[ky.ids], tids[ky.ids], rids[ky.ids], eids[ky.ids], gids[ky.ids], symbols[ky.ids], det.pvals[ky.ids])
limma.annot<-cbind(limma.ids,pids[limma.ids], tids[limma.ids], rids[limma.ids],eids[limma.ids], gids[limma.ids], symbols[limma.ids], limma.det.pvals[limma.ids])

colnames(ky.annot)<-col.names
colnames(limma.annot)<-col.names


# Normalised Expression Measures

ky.xpn<-cbind(ky.ns5[ky.ids], ky.astro[ky.ids], ky.xpn[ky.ids,])  #don't bother w/fc, it's meaningless
colnames(ky.xpn)<-c("ave_ns5","ave_astro", "ns1", "astro1", "ns2","astro2","ns3","astro3","ns4","astro4")

limma.xpn<-cbind(limma.ns5[limma.ids], limma.astro[limma.ids], limma.fc[limma.ids], limma.xpn[limma.ids,])
colnames(limma.xpn)<-c("log2(ave_ns5)","log2(ave_astro)","log2(astro/ns5)","ns1", "astro1", "ns2","astro2","ns3","astro3","ns4","astro4")

raw.xpn<-cbind(raw.ids, raw.ns5[raw.ids], raw.astro[raw.ids], raw.fc[raw.ids], raw.xpn[raw.ids,])
colnames(raw.xpn)<-c("ArrayAddressID","ave_ns5","ave_astro","astro/ns5","ns1", "astro1", "ns2","astro2","ns3","astro3","ns4","astro4")


ky<-cbind(ky.stats,ky.annot,ky.xpn )
limma<-cbind(limma.stats,limma.annot, limma.xpn)

stop("here")

                            
#write everything out as comma sep values files
write.table(ky, file="v3/results/ky.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(limma, file="v3/results/limma.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(raw.xpn, file="v3/results/raw_xpn.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")


#up & down sep
write.table(ky[up.ky,], file="v3/results/ky_up.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(ky[down.ky,], file="v3/results/ky_down.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")

write.table(limma[up.limma,], file="v3/results/limma_up.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(limma[down.limma,], file="v3/results/limma_down.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")




####
# And the same but filtered on detection and including fdr


ky.ids<-names(det.fdr)
limma.ids<-names(limma.det.fdr)
raw.ids<-detected

# Statistics

ky.stats<-cbind(det.tstats[ky.ids], det.pvals[ky.ids], det.fdr[ky.ids])
limma.stats<-cbind(limma.det.tstats[limma.ids], limma.det.pvals[limma.ids], limma.det.fdr[limma.ids])

colnames(ky.stats)<-c("t", "p","FDR")
colnames(limma.stats)<-c("mod-t","p", "FDR")


# Annotation

col.names<-c("Array Address ID","ProbeID", "TargetID", "RefSeq","EntrezGene","GI","GeneSymbol", "Detection Pvalue")

ky.annot<-cbind(ky.ids, pids[ky.ids], tids[ky.ids], rids[ky.ids], eids[ky.ids], gids[ky.ids], symbols[ky.ids], det.pvals[ky.ids])
limma.annot<-cbind(limma.ids,pids[limma.ids], tids[limma.ids], rids[limma.ids],eids[limma.ids], gids[limma.ids], symbols[limma.ids], limma.det.pvals[limma.ids])

colnames(ky.annot)<-col.names
colnames(limma.annot)<-col.names


# Normalised Expression Measures

ky.xpn<-exprs(kee.yew)[ky.ids,]
limma.xpn<-exprs(biocbead)[limma.ids,]
raw.xpn<-exprs(BSData)[detected,]


ky.xpn<-cbind(ky.ns5[ky.ids], ky.astro[ky.ids], ky.xpn[ky.ids,])  #don't bother w/fc, it's meaningless
colnames(ky.xpn)<-c("ave_ns5","ave_astro", "ns1", "astro1", "ns2","astro2","ns3","astro3","ns4","astro4")

limma.xpn<-cbind(limma.ns5[limma.ids], limma.astro[limma.ids], limma.fc[limma.ids], limma.xpn[limma.ids,])
colnames(limma.xpn)<-c("log2(ave_ns5)","log2(ave_astro)","log2(astro/ns5)","ns1", "astro1", "ns2","astro2","ns3","astro3","ns4","astro4")

raw.xpn<-cbind(raw.ids, raw.ns5[raw.ids], raw.astro[raw.ids], raw.fc[raw.ids], raw.xpn[raw.ids,])
colnames(raw.xpn)<-c("ArrayAddressID","ave_ns5","ave_astro","astro/ns5","ns1", "astro1", "ns2","astro2","ns3","astro3","ns4","astro4")

ky<-cbind(ky.stats,ky.annot,ky.xpn )
limma<-cbind(limma.stats,limma.annot, limma.xpn)

#write everything out as comma sep values files
write.table(ky, file="v3/results/ky_detected.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(limma, file="v3/results/limma_detected.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(raw.xpn, file="v3/results/raw_detected_xpn.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")



#up & down sep
write.table(ky[intersect(up.ky,detected),], file="v3/results/ky_detected_up.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(ky[intersect(down.ky,detected),], file="v3/results/ky_detected_down.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")

write.table(limma[intersect(up.limma,detected),], file="v3/results/limma_detected_up.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(limma[intersect(down.limma,detected),], file="v3/results/limma_detected_down.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")





#######
# Now save various IDs and subsets of IDs for external programs.

# name the rows of the annotation table with AAID so that we can subset it
rownames(annot)<-annot$Array_Address_Id



###
# Array Address ID

write.table(rownames(bh.fdr), file="v3/results/keeyew_bh_fdr_ordered_array_address_id.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(keeyew.fdr.05, file="v3/results/keeyew_fdr0.05_array_address_ids.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(keeyew.fdr.05.detected, file="v3/results/keeyew_fdr_05_detected_array_address_ids.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#up
write.table(intersect(rownames(bh.fdr), up.ky), file="v3/results/keeyew_bh_fdr_ordered_array_address_id_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(intersect(keeyew.fdr.05, up.ky), file="v3/results/keeyew_fdr0.05_array_address_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(intersect(keeyew.fdr.05.detected, up.ky), file="v3/results/keeyew_fdr_05_detected_array_address_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#down
write.table(intersect(rownames(bh.fdr), down.ky), file="v3/results/keeyew_bh_fdr_ordered_array_address_id_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(intersect(keeyew.fdr.05, down.ky), file="v3/results/keeyew_fdr0.05_array_address_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(intersect(keeyew.fdr.05.detected, down.ky), file="v3/results/keeyew_fdr_05_detected_array_address_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)


#astro
write.table(names(ky.astro.ord), file="v3/results/keeyew_astro_ord.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(ky.astro.1000, file="v3/results/keeyew_astro_top_1000.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(ky.astro.2000, file="v3/results/keeyew_astro_top_2000.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(names(ky.ns5.ord), file="v3/results/keeyew_ns5_ord.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(ky.ns5.1000, file="v3/results/keeyew_ns5_top_1000.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(ky.ns5.2000, file="v3/results/keeyew_ns5_top_2000.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )




write.table(limma.genelist$ID, file="v3/results/limma_bh_fdr_ordered_array_address_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(limma.fdr.05, file="v3/results/limma_fdr0.05_array_address_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(limma.fdr.05.detected, file="v3/results/limma_fdr_05_detected_array_address_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)


#up
write.table(intersect(limma.genelist$ID, up.limma), file="v3/results/limma_bh_fdr_ordered_array_address_id_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(intersect(limma.fdr.05, up.limma), file="v3/results/limma_fdr0.05_array_address_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(intersect(limma.fdr.05.detected, up.limma), file="v3/results/limma_fdr_05_detected_array_address_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#down
write.table(intersect(limma.genelist$ID, down.limma), file="v3/results/limma_bh_fdr_ordered_array_address_id_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(intersect(limma.fdr.05, down.limma), file="v3/results/limma_fdr0.05_array_address_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(intersect(limma.fdr.05.detected, down.limma), file="v3/results/limma_fdr_05_detected_array_address_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)


#astro
write.table(names(limma.astro.ord), file="v3/results/limma_astro_ord.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(limma.astro.1000, file="v3/results/limma_astro_top_1000.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(limma.astro.2000, file="v3/results/limma_astro_top_2000.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(names(limma.ns5.ord), file="v3/results/limma_ns5_ord.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(limma.ns5.1000, file="v3/results/limma_ns5_top_1000.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(limma.ns5.2000, file="v3/results/limma_ns5_top_2000.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )







###
# Probe IDs
write.table(pids[rownames(exprs(BSData))], file="v3/results/population_probe_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[detected], file="v3/results/population_detected_probe_ids.csv",row.names=FALSE, col.names=FALSE,quote=FALSE)


#fdr
write.table(pids[rownames(bh.fdr)], file="v3/results/keeyew_bh_fdr_ordered_probe_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[keeyew.fdr.05], file="v3/results/keeyew_fdr0.05_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[keeyew.fdr.05.detected], file="v3/results/keeyew_fdr_05_detected_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



#up
write.table(pids[intersect(rownames(bh.fdr), up.ky)], file="v3/results/keeyew_bh_fdr_ordered_probe_id_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(pids[intersect(keeyew.fdr.05, up.ky)], file="v3/results/keeyew_fdr0.05_prob_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(pids[intersect(keeyew.fdr.05.detected, up.ky)], file="v3/results/keeyew_fdr_05_detected_probe_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#down
write.table(pids[intersect(rownames(bh.fdr), down.ky)], file="v3/results/keeyew_bh_fdr_ordered_probe_id_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(pids[intersect(keeyew.fdr.05, down.ky)], file="v3/results/keeyew_fdr0.05_probe_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(pids[intersect(keeyew.fdr.05.detected, down.ky)], file="v3/results/keeyew_fdr_05_detected_probe_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)


#astro
write.table(pids[names(ky.astro.ord)], file="v3/results/keeyew_astro_ord_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[ky.astro.1000], file="v3/results/keeyew_astro_top_1000_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[ky.astro.2000], file="v3/results/keeyew_astro_top_2000_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(pids[names(ky.ns5.ord)], file="v3/results/keeyew_ns5_ord_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[ky.ns5.1000], file="v3/results/keeyew_ns5_top_1000_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[ky.ns5.2000], file="v3/results/keeyew_ns5_top_2000_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )


#fdr
write.table(pids[limma.genelist$ID], file="v3/results/limma_bh_fdr_ordered_probe_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[limma.fdr.05], file="v3/results/limma_fdr0.05_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[limma.fdr.05.detected], file="v3/results/limma_fdr_05_deted_probe_ids_csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



#up
write.table(pids[intersect(limma.genelist$ID, up.limma)], file="v3/results/limma_bh_fdr_ordered_probe_id_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(pids[intersect(limma.fdr.05, up.limma)], file="v3/results/limma_fdr0.05_probe_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(pids[intersect(limma.fdr.05.detected, up.limma)], file="v3/results/limma_fdr_05_detected_probe_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#down
write.table(pids[intersect(limma.genelist$ID, down.limma)], file="v3/results/limma_bh_fdr_ordered_probe_id_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(pids[intersect(limma.fdr.05, down.limma)], file="v3/results/limma_fdr0.05_probe_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(pids[intersect(limma.fdr.05.detected, down.limma)], file="v3/results/limma_fdr_05_detected_probe_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)


#astro
write.table(pids[names(limma.astro.ord)], file="v3/results/limma_astro_ord_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[limma.astro.1000], file="v3/results/limma_astro_top_1000_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[limma.astro.2000], file="v3/results/limma_astro_top_2000_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(pids[names(limma.ns5.ord)], file="v3/results/limma_ns5_ord_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[limma.ns5.1000], file="v3/results/limma_ns5_top_1000_probe_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pids[limma.ns5.2000], file="v3/results/limma_ns5_top_2000_probe.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )





###  The above have a unique value for each probe on the array, the others don't so we need to
### save a unique version too.


###
# Target IDs 


write.table(tids[rownames(exprs(BSData))], file="v3/results/population_target_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[rownames(exprs(BSData))]), file="v3/results/population_target_id_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[detected]), file="v3/results/population_detected_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

#fdr
write.table(tids[rownames(bh.fdr)], file="v3/results/keeyew_bh_fdr_ordered_target_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(tids[keeyew.fdr.05], file="v3/results/keeyew_fdr0.05_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[rownames(bh.fdr)]), file="v3/results/keeyew_bh_fdr_ordered_target_id_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[keeyew.fdr.05]), file="v3/results/keeyew_fdr0.05_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[keeyew.fdr.05.detected]), file="v3/results/keeyew_fdr_05_detected_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)


#up
write.table(tids[intersect(rownames(bh.fdr), up.ky)], file="v3/results/keeyew_bh_fdr_ordered_target_id_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(tids[intersect(keeyew.fdr.05, up.ky)], file="v3/results/keeyew_fdr0.05_target_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(tids[intersect(keeyew.fdr.05.detected, up.ky)], file="v3/results/keeyew_fdr_05_detected_target_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(tids[intersect(rownames(bh.fdr), up.ky)]), file="v3/results/keeyew_bh_fdr_ordered_target_id_up_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(tids[intersect(keeyew.fdr.05, up.ky)]), file="v3/results/keeyew_fdr0.05_target_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(tids[intersect(keeyew.fdr.05.detected, up.ky)]), file="v3/results/keeyew_fdr_05_detected_target_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)


#down
write.table(tids[intersect(rownames(bh.fdr), down.ky)], file="v3/results/keeyew_bh_fdr_ordered_target_id_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(tids[intersect(keeyew.fdr.05, down.ky)], file="v3/results/keeyew_fdr0.05_target_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(tids[intersect(keeyew.fdr.05.detected, down.ky)], file="v3/results/keeyew_fdr_05_detected_target_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(tids[intersect(rownames(bh.fdr), down.ky)]), file="v3/results/keeyew_bh_fdr_ordered_target_id_down_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(tids[intersect(keeyew.fdr.05, down.ky)]), file="v3/results/keeyew_fdr0.05_target_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(tids[intersect(keeyew.fdr.05.detected, down.ky)]), file="v3/results/keeyew_fdr_05_detected_target_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)



#astro
write.table(tids[names(ky.astro.ord)], file="v3/results/keeyew_astro_ord_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(tids[ky.astro.1000], file="v3/results/keeyew_astro_top_1000_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(tids[ky.astro.2000], file="v3/results/keeyew_astro_top_2000_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(tids[names(ky.astro.ord)]), file="v3/results/keeyew_astro_ord_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[ky.astro.1000]), file="v3/results/keeyew_astro_top_1000_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[ky.astro.2000]), file="v3/results/keeyew_astro_top_2000_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(tids[names(ky.ns5.ord)], file="v3/results/keeyew_ns5_ord_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(tids[ky.ns5.1000], file="v3/results/keeyew_ns5_top_1000_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(tids[ky.ns5.2000], file="v3/results/keeyew_ns5_top_2000_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(tids[names(ky.ns5.ord)]), file="v3/results/keeyew_ns5_ord_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[ky.ns5.1000]), file="v3/results/keeyew_ns5_top_1000_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[ky.ns5.2000]), file="v3/results/keeyew_ns5_top_2000_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )




#fdr
write.table(tids[limma.genelist$ID], file="v3/results/limma_bh_fdr_ordered_target_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(tids[limma.fdr.05], file="v3/results/limma_fdr0.05_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[limma.genelist$ID]), file="v3/results/limma_bh_fdr_ordered_target_id_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[limma.fdr.05]), file="v3/results/limma_fdr0.05_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[limma.fdr.05.detected]), file="v3/results.limma_fdr_05_detected_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



#up
write.table(tids[intersect(limma.genelist$ID, up.limma)], file="v3/results/limma_bh_fdr_ordered_target_id_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(tids[intersect(limma.fdr.05, up.limma)], file="v3/results/limma_fdr0.05_target_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(tids[intersect(limma.fdr.05.detected, up.limma)], file="v3/results/limma_fdr_05_detected_target_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(tids[intersect(limma.genelist$ID, up.limma)]), file="v3/results/limma_bh_fdr_ordered_target_id_up_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(tids[intersect(limma.fdr.05, up.limma)]), file="v3/results/limma_fdr0.05_target_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(tids[intersect(limma.fdr.05.detected, up.limma)]), file="v3/results/limma_fdr_05_detected_target_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#down
write.table(tids[intersect(limma.genelist$ID, down.limma)], file="v3/results/limma_bh_fdr_ordered_target_id_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(tids[intersect(limma.fdr.05, down.limma)], file="v3/results/limma_fdr0.05_target_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(tids[intersect(limma.fdr.05.detected, down.limma)], file="v3/results/limma_fdr_05_detected_target_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(tids[intersect(limma.genelist$ID, down.limma)]), file="v3/results/limma_bh_fdr_ordered_target_id_down_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(tids[intersect(limma.fdr.05, down.limma)]), file="v3/results/limma_fdr0.05_target_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(tids[intersect(limma.fdr.05.detected, down.limma)]), file="v3/results/limma_fdr_05_detected_target_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#astro
write.table(tids[names(limma.astro.ord)], file="v3/results/limma_astro_ord_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(tids[limma.astro.1000], file="v3/results/limma_astro_top_1000_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(tids[limma.astro.2000], file="v3/results/limma_astro_top_2000_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(tids[names(limma.astro.ord)]), file="v3/results/limma_astro_ord_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[limma.astro.1000]), file="v3/results/limma_astro_top_1000_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[limma.astro.2000]), file="v3/results/limma_astro_top_2000_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(tids[names(limma.ns5.ord)], file="v3/results/limma_ns5_ord_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(tids[limma.ns5.1000], file="v3/results/limma_ns5_top_1000_target_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(tids[limma.ns5.2000], file="v3/results/limma_ns5_top_2000_target.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(tids[names(limma.ns5.ord)]), file="v3/results/limma_ns5_ord_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[limma.ns5.1000]), file="v3/results/limma_ns5_top_1000_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(tids[limma.ns5.2000]), file="v3/results/limma_ns5_top_2000_target_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )







###
# RefSeq


write.table(rids[rownames(exprs(BSData))], file="v3/results/population_refseq_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[rownames(exprs(BSData))]), file="v3/results/population_refseq_id_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[detected]), file="v3/results/population_detected_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

#fdr
write.table(rids[rownames(bh.fdr)], file="v3/results/keeyew_bh_fdr_ordered_refseq_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(rids[keeyew.fdr.05], file="v3/results/keeyew_fdr0.05_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[rownames(bh.fdr)]), file="v3/results/keeyew_bh_fdr_ordered_refseq_id_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[keeyew.fdr.05]), file="v3/results/keeyew_fdr0.05_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[keeyew.fdr.05.detected]), file="v3/results/keeyew_fdr_05_detected_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)


#up
write.table(rids[intersect(rownames(bh.fdr), up.ky)], file="v3/results/keeyew_bh_fdr_ordered_refseq_id_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(rids[intersect(keeyew.fdr.05, up.ky)], file="v3/results/keeyew_fdr0.05_refseq_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(rids[intersect(keeyew.fdr.05.detected, up.ky)], file="v3/results/keeyew_fdr_05_detected_refseq_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(rids[intersect(rownames(bh.fdr), up.ky)]), file="v3/results/keeyew_bh_fdr_ordered_refseq_id_up_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(rids[intersect(keeyew.fdr.05, up.ky)]), file="v3/results/keeyew_fdr0.05_refseq_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(rids[intersect(keeyew.fdr.05.detected, up.ky)]), file="v3/results/keeyew_fdr_05_detected_refseq_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)


#down
write.table(rids[intersect(rownames(bh.fdr), down.ky)], file="v3/results/keeyew_bh_fdr_ordered_refseq_id_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(rids[intersect(keeyew.fdr.05, down.ky)], file="v3/results/keeyew_fdr0.05_refseq_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(rids[intersect(keeyew.fdr.05.detected, down.ky)], file="v3/results/keeyew_fdr_05_detected_refseq_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(rids[intersect(rownames(bh.fdr), down.ky)]), file="v3/results/keeyew_bh_fdr_ordered_refseq_id_down_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(rids[intersect(keeyew.fdr.05, down.ky)]), file="v3/results/keeyew_fdr0.05_refseq_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(rids[intersect(keeyew.fdr.05.detected, down.ky)]), file="v3/results/keeyew_fdr_05_detected_refseq_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)



#astro
write.table(rids[names(ky.astro.ord)], file="v3/results/keeyew_astro_ord_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(rids[ky.astro.1000], file="v3/results/keeyew_astro_top_1000_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(rids[ky.astro.2000], file="v3/results/keeyew_astro_top_2000_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(rids[names(ky.astro.ord)]), file="v3/results/keeyew_astro_ord_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[ky.astro.1000]), file="v3/results/keeyew_astro_top_1000_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[ky.astro.2000]), file="v3/results/keeyew_astro_top_2000_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(rids[names(ky.ns5.ord)], file="v3/results/keeyew_ns5_ord_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(rids[ky.ns5.1000], file="v3/results/keeyew_ns5_top_1000_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(rids[ky.ns5.2000], file="v3/results/keeyew_ns5_top_2000_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(rids[names(ky.ns5.ord)]), file="v3/results/keeyew_ns5_ord_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[ky.ns5.1000]), file="v3/results/keeyew_ns5_top_1000_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[ky.ns5.2000]), file="v3/results/keeyew_ns5_top_2000_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )







#fdr
write.table(rids[limma.genelist$ID], file="v3/results/limma_bh_fdr_ordered_refseq_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(rids[limma.fdr.05], file="v3/results/limma_fdr0.05_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[limma.genelist$ID]), file="v3/results/limma_bh_fdr_ordered_refseq_id_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[limma.fdr.05]), file="v3/results/limma_fdr0.05_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[limma.fdr.05.detected]), file="v3/results/limma_fdr_05_detected_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



#up
write.table(rids[intersect(limma.genelist$ID, up.limma)], file="v3/results/limma_bh_fdr_ordered_refseq_id_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(rids[intersect(limma.fdr.05, up.limma)], file="v3/results/limma_fdr0.05_refseq_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(rids[intersect(limma.fdr.05.detected, up.limma)], file="v3/results/limma_fdr_05_detected_refseq_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(rids[intersect(limma.genelist$ID, up.limma)]), file="v3/results/limma_bh_fdr_ordered_refseq_id_up_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(rids[intersect(limma.fdr.05, up.limma)]), file="v3/results/limma_fdr0.05_refseq_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(rids[intersect(limma.fdr.05.detected, up.limma)]), file="v3/results/limma_fdr_05_detected_refseq_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#down
write.table(rids[intersect(limma.genelist$ID, down.limma)], file="v3/results/limma_bh_fdr_ordered_refseq_id_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(rids[intersect(limma.fdr.05, down.limma)], file="v3/results/limma_fdr0.05_refseq_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(rids[intersect(limma.fdr.05.detected, down.limma)], file="v3/results/limma_fdr_05_detected_refseq_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(rids[intersect(limma.genelist$ID, down.limma)]), file="v3/results/limma_bh_fdr_ordered_refseq_id_down_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(rids[intersect(limma.fdr.05, down.limma)]), file="v3/results/limma_fdr0.05_refseq_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(rids[intersect(limma.fdr.05.detected, down.limma)]), file="v3/results/limma_fdr_05_detected_refseq_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#astro
write.table(rids[names(limma.astro.ord)], file="v3/results/limma_astro_ord_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(rids[limma.astro.1000], file="v3/results/limma_astro_top_1000_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(rids[limma.astro.2000], file="v3/results/limma_astro_top_2000_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(rids[names(limma.astro.ord)]), file="v3/results/limma_astro_ord_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[limma.astro.1000]), file="v3/results/limma_astro_top_1000_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[limma.astro.2000]), file="v3/results/limma_astro_top_2000_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(rids[names(limma.ns5.ord)], file="v3/results/limma_ns5_ord_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(rids[limma.ns5.1000], file="v3/results/limma_ns5_top_1000_refseq_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(rids[limma.ns5.2000], file="v3/results/limma_ns5_top_2000_refseq.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(rids[names(limma.ns5.ord)]), file="v3/results/limma_ns5_ord_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[limma.ns5.1000]), file="v3/results/limma_ns5_top_1000_refseq_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(rids[limma.ns5.2000]), file="v3/results/limma_ns5_top_2000_refseq_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )





###
# EntrezGene

write.table(eids[rownames(exprs(BSData))], file="v3/results/population_entrezgene_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[rownames(exprs(BSData))]), file="v3/results/population_entrezgene_id_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[detected]), file="v3/results/population_detected_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(eids[rownames(bh.fdr)], file="v3/results/keeyew_bh_fdr_ordered_entrezgene_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(eids[keeyew.fdr.05], file="v3/results/keeyew_fdr0.05_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[rownames(bh.fdr)]), file="v3/results/keeyew_bh_fdr_ordered_entrezgene_id_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[keeyew.fdr.05]), file="v3/results/keeyew_fdr0.05_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[keeyew.fdr.05.detected]), file="v3/results/keeyew_fdr_05_detected_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



#up
write.table(eids[intersect(rownames(bh.fdr), up.ky)], file="v3/results/keeyew_bh_fdr_ordered_entrezgene_id_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(eids[intersect(keeyew.fdr.05, up.ky)], file="v3/results/keeyew_fdr0.05_entrezgene_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(eids[intersect(keeyew.fdr.05.detected, up.ky)], file="v3/results/keeyew_fdr_05_detected_entrezgene_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(eids[intersect(rownames(bh.fdr), up.ky)]), file="v3/results/keeyew_bh_fdr_ordered_entrezgene_id_up_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(eids[intersect(keeyew.fdr.05, up.ky)]), file="v3/results/keeyew_fdr0.05_entrezgene_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(eids[intersect(keeyew.fdr.05.detected, up.ky)]), file="v3/results/keeyew_fdr_05_detected_entrezgene_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)


#down
write.table(eids[intersect(rownames(bh.fdr), down.ky)], file="v3/results/keeyew_bh_fdr_ordered_entrezgene_id_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(eids[intersect(keeyew.fdr.05, down.ky)], file="v3/results/keeyew_fdr0.05_entrezgene_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(eids[intersect(keeyew.fdr.05.detected, down.ky)], file="v3/results/keeyew_fdr_05_detected_entrezgene_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(eids[intersect(rownames(bh.fdr), down.ky)]), file="v3/results/keeyew_bh_fdr_ordered_entrezgene_id_down_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(eids[intersect(keeyew.fdr.05, down.ky)]), file="v3/results/keeyew_fdr0.05_entrezgene_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(eids[intersect(keeyew.fdr.05.detected, down.ky)]), file="v3/results/keeyew_fdr_05_detected_entrezgene_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)



#astro
write.table(eids[names(ky.astro.ord)], file="v3/results/keeyew_astro_ord_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(eids[ky.astro.1000], file="v3/results/keeyew_astro_top_1000_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(eids[ky.astro.2000], file="v3/results/keeyew_astro_top_2000_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(eids[names(ky.astro.ord)]), file="v3/results/keeyew_astro_ord_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[ky.astro.1000]), file="v3/results/keeyew_astro_top_1000_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[ky.astro.2000]), file="v3/results/keeyew_astro_top_2000_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(eids[names(ky.ns5.ord)], file="v3/results/keeyew_ns5_ord_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(eids[ky.ns5.1000], file="v3/results/keeyew_ns5_top_1000_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(eids[ky.ns5.2000], file="v3/results/keeyew_ns5_top_2000_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(eids[names(ky.ns5.ord)]), file="v3/results/keeyew_ns5_ord_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[ky.ns5.1000]), file="v3/results/keeyew_ns5_top_1000_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[ky.ns5.2000]), file="v3/results/keeyew_ns5_top_2000_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )






write.table(eids[limma.genelist$ID], file="v3/results/limma_bh_fdr_ordered_entrezgene_id.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(eids[limma.fdr.05], file="v3/results/limma_fdr0.05_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[limma.genelist$ID]), file="v3/results/limma_bh_fdr_ordered_entrezgene_id_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[limma.fdr.05]), file="v3/results/limma_fdr0.05_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[limma.fdr.05.detected]), file="v3/results/limma_fdr_05_detected_entrezgene_ids_unique", row.names=FALSE, col.names=FALSE, quote=FALSE)


#up
write.table(eids[intersect(limma.genelist$ID, up.limma)], file="v3/results/limma_bh_fdr_ordered_entrezgene_id_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(eids[intersect(limma.fdr.05, up.limma)], file="v3/results/limma_fdr0.05_entrezgene_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(eids[intersect(limma.fdr.05.detected, up.limma)], file="v3/results/limma_fdr_05_detected_entrezgene_ids_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(eids[intersect(limma.genelist$ID, up.limma)]), file="v3/results/limma_bh_fdr_ordered_entrezgene_id_up_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(eids[intersect(limma.fdr.05, up.limma)]), file="v3/results/limma_fdr0.05_entrezgene_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(eids[intersect(limma.fdr.05.detected, up.limma)]), file="v3/results/limma_fdr_05_detected_entrezgene_ids_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#down
write.table(eids[intersect(limma.genelist$ID, down.limma)], file="v3/results/limma_bh_fdr_ordered_entrezgene_id_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(eids[intersect(limma.fdr.05, down.limma)], file="v3/results/limma_fdr0.05_entrezgene_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(eids[intersect(limma.fdr.05.detected, down.limma)], file="v3/results/limma_fdr_05_detected_entrezgene_ids_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(eids[intersect(limma.genelist$ID, down.limma)]), file="v3/results/limma_bh_fdr_ordered_entrezgene_id_down_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(eids[intersect(limma.fdr.05, down.limma)]), file="v3/results/limma_fdr0.05_entrezgene_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(eids[intersect(limma.fdr.05.detected, down.limma)]), file="v3/results/limma_fdr_05_detected_entrezgene_ids_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#astro
write.table(eids[names(limma.astro.ord)], file="v3/results/limma_astro_ord_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(eids[limma.astro.1000], file="v3/results/limma_astro_top_1000_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(eids[limma.astro.2000], file="v3/results/limma_astro_top_2000_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(eids[names(limma.astro.ord)]), file="v3/results/limma_astro_ord_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[limma.astro.1000]), file="v3/results/limma_astro_top_1000_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[limma.astro.2000]), file="v3/results/limma_astro_top_2000_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(eids[names(limma.ns5.ord)], file="v3/results/limma_ns5_ord_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(eids[limma.ns5.1000], file="v3/results/limma_ns5_top_1000_entrezgene_ids.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(eids[limma.ns5.2000], file="v3/results/limma_ns5_top_2000_entrezgene.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(eids[names(limma.ns5.ord)]), file="v3/results/limma_ns5_ord_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[limma.ns5.1000]), file="v3/results/limma_ns5_top_1000_entrezgene_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(eids[limma.ns5.2000]), file="v3/results/limma_ns5_top_2000_entrezgene_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )






###
# GI


write.table(gids[rownames(exprs(BSData))], file="v3/results/population_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[rownames(exprs(BSData))]), file="v3/results/population_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[detected]), file="v3/results/population_detected_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)


write.table(gids[rownames(bh.fdr)], file="v3/results/keeyew_bh_fdr_ordered_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gids[keeyew.fdr.05], file="v3/results/keeyew_fdr0.05_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[rownames(bh.fdr)]), file="v3/results/keeyew_bh_fdr_ordered_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[keeyew.fdr.05]), file="v3/results/keeyew_fdr0.05_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[keeyew.fdr.05.detected]), file="v3/results/keeyew_fdr_05_detected_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)


#up
write.table(gids[intersect(rownames(bh.fdr), up.ky)], file="v3/results/keeyew_bh_fdr_ordered_GI_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(gids[intersect(keeyew.fdr.05, up.ky)], file="v3/results/keeyew_fdr0.05_GI_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(gids[intersect(keeyew.fdr.05.detected, up.ky)], file="v3/results/keeyew_fdr_05_detected_GI_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(gids[intersect(rownames(bh.fdr), up.ky)]), file="v3/results/keeyew_bh_fdr_ordered_GI_up_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(gids[intersect(keeyew.fdr.05, up.ky)]), file="v3/results/keeyew_fdr0.05_GI_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(gids[intersect(keeyew.fdr.05.detected, up.ky)]), file="v3/results/keeyew_fdr_05_detected_GI_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)


#down
write.table(gids[intersect(rownames(bh.fdr), down.ky)], file="v3/results/keeyew_bh_fdr_ordered_GI_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(gids[intersect(keeyew.fdr.05, down.ky)], file="v3/results/keeyew_fdr0.05_GI_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(gids[intersect(keeyew.fdr.05.detected, down.ky)], file="v3/results/keeyew_fdr_05_detected_GI_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(gids[intersect(rownames(bh.fdr), down.ky)]), file="v3/results/keeyew_bh_fdr_ordered_GI_down_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(gids[intersect(keeyew.fdr.05, down.ky)]), file="v3/results/keeyew_fdr0.05_GI_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(gids[intersect(keeyew.fdr.05.detected, down.ky)]), file="v3/results/keeyew_fdr_05_detected_GI_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)



#astro
write.table(gids[names(ky.astro.ord)], file="v3/results/keeyew_astro_ord_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gids[ky.astro.1000], file="v3/results/keeyew_astro_top_1000_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gids[ky.astro.2000], file="v3/results/keeyew_astro_top_2000_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(gids[names(ky.astro.ord)]), file="v3/results/keeyew_astro_ord_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[ky.astro.1000]), file="v3/results/keeyew_astro_top_1000_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[ky.astro.2000]), file="v3/results/keeyew_astro_top_2000_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(gids[names(ky.ns5.ord)], file="v3/results/keeyew_ns5_ord_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gids[ky.ns5.1000], file="v3/results/keeyew_ns5_top_1000_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gids[ky.ns5.2000], file="v3/results/keeyew_ns5_top_2000_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(gids[names(ky.ns5.ord)]), file="v3/results/keeyew_ns5_ord_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[ky.ns5.1000]), file="v3/results/keeyew_ns5_top_1000_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[ky.ns5.2000]), file="v3/results/keeyew_ns5_top_2000_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )







write.table(gids[limma.genelist$ID], file="v3/results/limma_bh_fdr_ordered_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gids[limma.fdr.05], file="v3/results/limma_fdr0.05_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[limma.genelist$ID]), file="v3/results/limma_bh_fdr_ordered_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[limma.fdr.05]), file="v3/results/limma_fdr0.05_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[limma.fdr.05.detected]), file="v3/results/limma_fdr_05_detected_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



#up
write.table(gids[intersect(limma.genelist$ID, up.limma)], file="v3/results/limma_bh_fdr_ordered_GI_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(gids[intersect(limma.fdr.05, up.limma)], file="v3/results/limma_fdr0.05_GI_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(gids[intersect(limma.fdr.05.detected, up.limma)], file="v3/results/limma_fdr_05_detected_GI_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(gids[intersect(limma.genelist$ID, up.limma)]), file="v3/results/limma_bh_fdr_ordered_GI_up_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(gids[intersect(limma.fdr.05, up.limma)]), file="v3/results/limma_fdr0.05_GI_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(gids[intersect(limma.fdr.05.detected, up.limma)]), file="v3/results/limma_fdr_05_detected_GI_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#down
write.table(gids[intersect(limma.genelist$ID, down.limma)], file="v3/results/limma_bh_fdr_ordered_GI_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(gids[intersect(limma.fdr.05, down.limma)], file="v3/results/limma_fdr0.05_GI_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(gids[intersect(limma.fdr.05.detected, down.limma)], file="v3/results/limma_fdr_05_detected_GI_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(gids[intersect(limma.genelist$ID, down.limma)]), file="v3/results/limma_bh_fdr_ordered_GI_down_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(gids[intersect(limma.fdr.05, down.limma)]), file="v3/results/limma_fdr0.05_GI_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(gids[intersect(limma.fdr.05.detected, down.limma)]), file="v3/results/limma_fdr_05_detected_GI_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#astro
write.table(gids[names(limma.astro.ord)], file="v3/results/limma_astro_ord_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gids[limma.astro.1000], file="v3/results/limma_astro_top_1000_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gids[limma.astro.2000], file="v3/results/limma_astro_top_2000_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(gids[names(limma.astro.ord)]), file="v3/results/limma_astro_ord_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[limma.astro.1000]), file="v3/results/limma_astro_top_1000_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[limma.astro.2000]), file="v3/results/limma_astro_top_2000_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(gids[names(limma.ns5.ord)], file="v3/results/limma_ns5_ord_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gids[limma.ns5.1000], file="v3/results/limma_ns5_top_1000_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gids[limma.ns5.2000], file="v3/results/limma_ns5_top_2000_GI.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(gids[names(limma.ns5.ord)]), file="v3/results/limma_ns5_ord_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[limma.ns5.1000]), file="v3/results/limma_ns5_top_1000_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(gids[limma.ns5.2000]), file="v3/results/limma_ns5_top_2000_GI_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )



###
# Gene Symbol

write.table(symbols[rownames(exprs(BSData))], file="v3/results/population_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[rownames(exprs(BSData))]), file="v3/results/population_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[detected]), file="v3/results/population_detected_target_ids_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(symbols[rownames(bh.fdr)], file="v3/results/keeyew_bh_fdr_ordered_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(symbols[keeyew.fdr.05], file="v3/results/keeyew_fdr0.05_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[rownames(bh.fdr)]), file="v3/results/keeyew_bh_fdr_ordered_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[keeyew.fdr.05]), file="v3/results/keeyew_fdr0.05_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[keeyew.fdr.05.detected]), file="v3/results/keeyew_fdr_05_detected_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)




#up
write.table(symbols[intersect(rownames(bh.fdr), up.ky)], file="v3/results/keeyew_bh_fdr_ordered_symbol_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(symbols[intersect(keeyew.fdr.05, up.ky)], file="v3/results/keeyew_fdr0.05_symbol_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(symbols[intersect(keeyew.fdr.05.detected, up.ky)], file="v3/results/keeyew_fdr_05_detected_symbol_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(symbols[intersect(rownames(bh.fdr), up.ky)]), file="v3/results/keeyew_bh_fdr_ordered_symbol_up_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(symbols[intersect(keeyew.fdr.05, up.ky)]), file="v3/results/keeyew_fdr0.05_symbol_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(symbols[intersect(keeyew.fdr.05.detected, up.ky)]), file="v3/results/keeyew_fdr_05_detected_symbol_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)


#down
write.table(symbols[intersect(rownames(bh.fdr), down.ky)], file="v3/results/keeyew_bh_fdr_ordered_symbol_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(symbols[intersect(keeyew.fdr.05, down.ky)], file="v3/results/keeyew_fdr0.05_symbol_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(symbols[intersect(keeyew.fdr.05.detected, down.ky)], file="v3/results/keeyew_fdr_05_detected_symbol_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(symbols[intersect(rownames(bh.fdr), down.ky)]), file="v3/results/keeyew_bh_fdr_ordered_symbol_down_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(symbols[intersect(keeyew.fdr.05, down.ky)]), file="v3/results/keeyew_fdr0.05_symbol_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(symbols[intersect(keeyew.fdr.05.detected, down.ky)]), file="v3/results/keeyew_fdr_05_detected_symbol_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)



#astro
write.table(symbols[names(ky.astro.ord)], file="v3/results/keeyew_astro_ord_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(symbols[ky.astro.1000], file="v3/results/keeyew_astro_top_1000_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(symbols[ky.astro.2000], file="v3/results/keeyew_astro_top_2000_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(symbols[names(ky.astro.ord)]), file="v3/results/keeyew_astro_ord_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[ky.astro.1000]), file="v3/results/keeyew_astro_top_1000_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[ky.astro.2000]), file="v3/results/keeyew_astro_top_2000_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(symbols[names(ky.ns5.ord)], file="v3/results/keeyew_ns5_ord_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(symbols[ky.ns5.1000], file="v3/results/keeyew_ns5_top_1000_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(symbols[ky.ns5.2000], file="v3/results/keeyew_ns5_top_2000_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(symbols[names(ky.ns5.ord)]), file="v3/results/keeyew_ns5_ord_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[ky.ns5.1000]), file="v3/results/keeyew_ns5_top_1000_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[ky.ns5.2000]), file="v3/results/keeyew_ns5_top_2000_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )







write.table(symbols[limma.genelist$ID], file="v3/results/limma_bh_fdr_ordered_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(symbols[limma.fdr.05], file="v3/results/limma_fdr0.05_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[limma.genelist$ID]), file="v3/results/limma_bh_fdr_ordered_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[limma.fdr.05]), file="v3/results/limma_fdr0.05_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[limma.fdr.05.detected]), file="v3/results/limma_fdr_05_detected_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)








#up
write.table(symbols[intersect(limma.genelist$ID, up.limma)], file="v3/results/limma_bh_fdr_ordered_symbol_up.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(symbols[intersect(limma.fdr.05, up.limma)], file="v3/results/limma_fdr0.05_symbol_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(symbols[intersect(limma.fdr.05.detected, up.limma)], file="v3/results/limma_fdr_05_detected_symbol_up.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(symbols[intersect(limma.genelist$ID, up.limma)]), file="v3/results/limma_bh_fdr_ordered_symbol_up_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(symbols[intersect(limma.fdr.05, up.limma)]), file="v3/results/limma_fdr0.05_symbol_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(symbols[intersect(limma.fdr.05.detected, up.limma)]), file="v3/results/limma_fdr_05_detected_symbol_up_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#down
write.table(symbols[intersect(limma.genelist$ID, down.limma)], file="v3/results/limma_bh_fdr_ordered_symbol_down.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(symbols[intersect(limma.fdr.05, down.limma)], file="v3/results/limma_fdr0.05_symbol_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(symbols[intersect(limma.fdr.05.detected, down.limma)], file="v3/results/limma_fdr_05_detected_symbol_down.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

write.table(unique(symbols[intersect(limma.genelist$ID, down.limma)]), file="v3/results/limma_bh_fdr_ordered_symbol_down_unique.csv", row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(unique(symbols[intersect(limma.fdr.05, down.limma)]), file="v3/results/limma_fdr0.05_symbol_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(unique(symbols[intersect(limma.fdr.05.detected, down.limma)]), file="v3/results/limma_fdr_05_detected_symbol_down_unique.csv", row.names=FALSE, quote=FALSE, col.names=FALSE)

#astro
write.table(symbols[names(limma.astro.ord)], file="v3/results/limma_astro_ord_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(symbols[limma.astro.1000], file="v3/results/limma_astro_top_1000_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(symbols[limma.astro.2000], file="v3/results/limma_astro_top_2000_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(symbols[names(limma.astro.ord)]), file="v3/results/limma_astro_ord_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[limma.astro.1000]), file="v3/results/limma_astro_top_1000_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[limma.astro.2000]), file="v3/results/limma_astro_top_2000_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

#ns5
write.table(symbols[names(limma.ns5.ord)], file="v3/results/limma_ns5_ord_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(symbols[limma.ns5.1000], file="v3/results/limma_ns5_top_1000_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(symbols[limma.ns5.2000], file="v3/results/limma_ns5_top_2000_symbol.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(unique(symbols[names(limma.ns5.ord)]), file="v3/results/limma_ns5_ord_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[limma.ns5.1000]), file="v3/results/limma_ns5_top_1000_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(symbols[limma.ns5.2000]), file="v3/results/limma_ns5_top_2000_symbol_unique.csv", row.names=FALSE, col.names=FALSE, quote=FALSE )



