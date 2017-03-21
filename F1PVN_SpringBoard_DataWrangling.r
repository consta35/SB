setwd("E:\\PVN_Sequencing\\NEW_PVN_2016_Aligned83\\New_aligned_F1_Tophat83")
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
## Increase memory limit to allow processes to occur 
memory.limit(size=800000000)

## Use GenomicAlignments to read in aligned sequencing files
readsF1Br1=readGAlignments("F1_Beta_SGM4_Tophat.bam")
readsF1Br2=readGAlignments("F1_Beta_SGM8_Tophat.bam")
readsF1Br3=readGAlignments("F1_Beta_SGM12_Tophat.bam")
readsF1Br4=readGAlignments("F1_Beta_SGM16_Tophat.bam")
readsF1Br5=readGAlignments("F1_Beta_SGM20_Tophat.bam")
readsF1Cr1=readGAlignments("F1_Control_SGM3_Tophat.bam")
readsF1Cr2=readGAlignments("F1_Control_SGM7_Tophat.bam")
readsF1Cr3=readGAlignments("F1_Control_SGM11_Tophat.bam")
readsF1Cr4=readGAlignments("F1_Control_SGM15_Tophat.bam")
readsF1Cr5=readGAlignments("F1_Control_SGM19_Tophat.bam")
txdb=makeTxDbFromGFF("UNEWcuffmerge83_grep.gtf", format = "gtf")
# Use the function transcriptsBy(txdb,'gene') for the whole genetic region instead of just exons
ex_by_gene<- exonsBy(txdb,'gene')

# Count how many sequenced reads overlap with exons
countsF1Br1= countOverlaps(ex_by_gene,readsF1Br1)
countsF1Br2= countOverlaps(ex_by_gene,readsF1Br2)
countsF1Br3= countOverlaps(ex_by_gene,readsF1Br3)
countsF1Br4= countOverlaps(ex_by_gene,readsF1Br4)
countsF1Br5= countOverlaps(ex_by_gene,readsF1Br5)
countsF1Cr1= countOverlaps(ex_by_gene,readsF1Cr1)
countsF1Cr2= countOverlaps(ex_by_gene,readsF1Cr2)
countsF1Cr3= countOverlaps(ex_by_gene,readsF1Cr3)
countsF1Cr4= countOverlaps(ex_by_gene,readsF1Cr4)
countsF1Cr5= countOverlaps(ex_by_gene,readsF1Cr5)
## Make Count Table
countTable<- data.frame(Control1=countsF1Cr1,Control2=countsF1Cr2, Control3=countsF1Cr3, Control4=countsF1Cr4,Control5=countsF1Cr5, Beta1=countsF1Br1,Beta2=countsF1Br2,Beta3=countsF1Br3,Beta4=countsF1Br4,Beta5=countsF1Br5, stringsAsFactors=FALSE)
x <- rowSums(countTable<=0)!=ncol(countTable)
newCountTable <- countTable[x,]
write.table(newCountTable, file="F1PVN_Rawcount_zero_elim83_NEW.csv", sep= ",")
data<-newCountTable
## Get PCA of data
pca<-prcomp(t(as.matrix(data)))
plot(pca$x)
text(pca$x[,1], pca$x[,2], colnames(data), pos=2)
filter <- apply(data, 1, function(x) length(x[x>5])>=5)
filtered <- data[filter,]
write.table(filtered, "filtered.csv", sep=",")

## remove Outliers
library (DESeq2)
pheno<- read.csv("phenof1.csv", header=TRUE, row.names=1)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(filtered),
colData = pheno,
design = ~ Treatment)
dds <- DESeq(dds)
res <- results(dds)
## find out how many outliers
summary(res)
## pull out Cooks Scores
maxCooks <- apply(assays(dds)[["cooks"]], 1, max)
write.table(maxCooks, "maxcooksf1pvn.csv", sep=",")

### Use Join to get gene names (I do this in Linux terminal)
cd /scratch/m/matthew7/consta35/PVN83_NEW_TOPHAT/PVN_NEW_F1_83_TOPHAT

# Convert files into linux compatible format
dos2unix F1_PVN_NEW_filtered_outliersRM.txt F1_PVN_NEW_GeneIDs83.txt

#Join files on XLOC code to create a file with gene names
#F1_PVN Gene IDs come from running Cuffdiff on F1 data using UNEWcuffmerge83_grep.gtf as reference. The Cuffdiff output gives XLOC with their matching gene names. F1_PVN_NEW_GeneIDs83.txt is a file with 2 columns from the cuffdiff output - the IDs and XLOC codes. 

join <(sort F1_PVN_NEW_GeneIDs83.txt) <(sort F1_PVN_NEW_filtered_outliersRM.txt) > F1_PVN_NEW_filtered_outliersRM_Names83_TRY5.txt

### Read output from Join back into R. In this data, we will have multiple rows with the same gene name- for ex. we will have multiple rows with counts for Y_RNA, and Sno. Melt and Cast allow us to sum all the counts that have the same gene name.
library(reshape)
library(reshape2)
data<- read.csv("F1_PVN_NEW_filtered_outliersRM_Names83_TRY5.csv", header=TRUE)
md<-melt(data, id= "Gene")
cd<-cast(md,formula=Gene~variable,sum)

### pull out genes that have at least 5 counts in at least 5 samples
filter <- apply(cd, 1, function(x) length(x[x>5])>=5)
filtered <- cd[filter,]

write.table(filtered, "Filtered_data.csv", sep=",")

