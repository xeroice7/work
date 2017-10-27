library(tidyverse)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ISLR)
library(gplots)
library("Biobase")

nan.na <- function (x) {
  x[is.nan(x)] <- NA
  return(x)
}

dist_cor <- function(x) {
  as.dist(1 - cor(t(x), method = "pearson"))
}

clus_wd2 <- function(x) {
  hclust(x, method = "ward.D2")
}

clus_wd <- function(x) {
  hclust(x, method = "ward.D")
}

######
# cBio - TCGA Cell 2015
#####
#### READ IN PATIENT DATA AND TIDY
patientdata <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/brca_tcga_pub/brca_tcga_pub/data_clinical.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

colnames(patientdata) <- lapply(patientdata[5,], as.character)

patientdata <- patientdata[-(1:5),]

patientdata <- subset(patientdata, select = c(PATIENT_ID, Tumor, `AJCC Stage`, Metastasis)) 

splits <- str_split_fixed(patientdata$`AJCC Stage`, " ", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '2'] <- 'Stages'
patientdata <- subset(patientdata, select=-1)

splits <- str_split_fixed(patientdata$Stage, "A|B|C", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'Stage'
patientdata <- subset(patientdata, select=-c(2, `AJCC Stage`))


#### READ IN CLINICAL DATA AND TIDY
expressiondata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/brca_tcga_pub/brca_tcga_pub/data_expression_median.csv", sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

expressiondata[expressiondata == "null"] <- NA

expressiondata <- subset(expressiondata, select = -Entrez_Gene_Id)

expressiondata <- expressiondata %>% drop_na()

expressionDF <- melt(expressiondata, id.vars = "Hugo_Symbol")

colnames(expressionDF) <- c("GENE", "PATIENT_ID", "EXPRESSION_LEVEL")

#### MERGING CLINICAL AND PATIENT DATA 
patientDF <- merge(x = expressionDF, y = patientdata, by = "PATIENT_ID")

#### BRINGING IN GENE CRITERIA WE WANT TO FILTER OUR DATASETS
#Total iPS overlaps (including somatic targets - 93 genes)
huvfips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/HUViPSvsFiPS/HUViPSFiPSOverlap.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
huvfips_overlap <- as.factor(huvfips_overlap$x)

#All Human Cell surface-specific targets indicated by GO analysis on GO site 
#data1 <- read.table("~/Desktop/HumanGOCellSurface.txt", sep = "" , header = F , nrows = 2000, na.strings ="", stringsAsFactors= F, fill = T) #Download genes from text file that were created by GO analysis of proteins 
                                                                                                                                          #that included "Cell Surface" as part of its GO
#write.csv(data1, "HumanGOCellSurface.csv") #Write surface data to a file making it a CSV

data <- read.csv("~/Desktop/HumanGOCellSurface.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) # Import CSV of surface expression
surface_genes <- data$V1 
surface_genes <- unique(surface_genes) 

#Unique iPS targets(no somatic source - 34 genes)
gene <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/CSVs/TotalUniqueiPSGENEandPID.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

# Examine stage differences between our MCF7/iPS overlap without a somatic filter (29 genes)
#iPS to MCF7 Comparisons (No Somatic Filters - can include somatic hits)
mcf7huvips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsHUViPS/HUViPSMCF7OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
mcf7fips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsFiPS/FiPSMCF7OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
mcf7ips_overlap <- intersect(mcf7fips_overlap$V1, mcf7huvips_overlap$V1)

#A smaller subset of genes pulled from FiPS and HUViPS MS sets to test clustering
testgenes <- read.csv("~/Desktop/TestGeneSetUnique.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) # Import CSV of surface expression
testgenes <- as.factor(testgenes$x)
c_testgenes <- as.character(testgenes)
c_gene <- as.character(gene)
totaltestgenes <- c(c_testgenes, c_gene)
totaltestgenes <- unique(totaltestgenes)
totaltestgenes <- as.factor(totaltestgenes)

#### FILTERING DATA WITH OUR CRITERIA OF INTEREST
meta <- factor(c("M0", "M1"))
meta <- ordered(meta, levels = c("M0", "M1"))
grades <- factor(c("T1", "T2", "T3"))
grades <- ordered(grades, levels = c("T1", "T2", "T3"))
stages <- factor(c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IIIC", "IV"))
stages <- ordered(stages, levels = c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IIIC", "IV"))
stage <- factor(c("I", "II", "III", "IV"))
stage <- ordered(stages, levels = c("I", "II", "III", "IV"))

#stage_diffs <- subset(patientDF, GENE %in% huvfips_overlap) #Subset rows that are matches to those in the targets
stage_diffs <- subset(patientDF, GENE %in% totaltestgenes)
#stage_diffs <- subset(stage_diffs, GENE %in% surface_genes)
stage_diffs <- subset(stage_diffs, Metastasis %in% meta)
stage_diffs <- subset(stage_diffs, Tumor %in% grades)
stage_diffs <- subset(stage_diffs, Stages %in% stages)
stage_diffs <- subset(stage_diffs, Stage %in% stage)
#stage_diffs <- filter(stage_diffs, ((Tumor == "T3") | (Tumor == "T2")))
#stage_diffs <- filter(stage_diffs, Stages == "IV")

patient <- dcast(stage_diffs, PATIENT_ID+Tumor+Stage+Stages+Metastasis ~ GENE, value.var = "EXPRESSION_LEVEL")

extras <- cbind.data.frame(patient$PATIENT_ID, patient$Metastasis, patient$Tumor, patient$Stages, patient$Stage)
colnames(extras) <- c("PATIENT_ID", "METASTASIS", "TUMOR", "STAGES", "STAGE")
#
splits <- str_split_fixed(patient$PATIENT_ID, pattern = "-", "3")

patient <- cbind.data.frame(splits, patient)

patient <- subset(patient, select = -c(`1`, `2`, PATIENT_ID)) # Getting rid of Stage Column to prepare for a calculation matrix

names(patient)[names(patient) == '3'] <- 'ID'

patient <- arrange(patient, Stages, Tumor, Metastasis)

patient <- patient %>% 
  unite(PATIENT_ID, ID, Stages, Tumor, Metastasis, sep = "_", remove = TRUE)

patient <- subset(patient, select = -Stage)

#
patient <- subset(patient, select = -c(Stages, Tumor, Metastasis, Stage))
patient[,1] <- as.character(patient[,1])

patient[,2:ncol(patient)] <- sapply(patient[,2:ncol(patient)],as.numeric)

sapply(patient, class)  #to check classes

patient <- patient %>% remove_rownames %>% column_to_rownames(var="PATIENT_ID") # Making the patient IDs the rownames, not the first column

redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)  # Making a color palette for the heatmap 

#Setting up the data for a heatmap/cluster analysis
patient <- patient[, apply(patient, 2, sum)!=0] ### Required to remove all the columns with 0 in them to get a distance matrix
test <- as.matrix(patient)
test1 <- t(test)

#heatmap(test1, Colv=NA, col=greenred(10),scale="none")
#cor(t(test1))
t <- 1-cor(t(test1))
hc <- hclust(as.dist(1-cor(t(test1))))
plot(hc)
heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=greenred(10), cexRow = 0.2)
heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=redblackgreen, cexRow = 0.2)

heatmap.2(test, 
          #distfun = dist_cor, hclustfun = clus_wd2, 
          # clustering
          #distfun = dist, 
          #hclust = clus,
          # scaling (genes are in rows)
          scale = "row",
          # color
          Rowv = FALSE,
          col = redblackgreen, 
          cexRow = 0.15,
          # labels
          #ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none")
 
X = t(scale(t(test),center=TRUE,scale=FALSE))
sv = t(X)
w <- which(is.na(as.matrix(sv)))
#sv <- sv[which(!is.finite(as.matrix(sv)))] <- 0
sv[w] <- 0
sv = svd(sv)
U = sv$u
V = sv$v
D = sv$d

aa<- grep("grey",colors())
bb<- grep("green",colors())
cc<-  grep("red",colors())
gcol2<- colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

## use the genes that drive the first PC1. This is the first major patter in the data
k=1
ord1<- order(abs(V[,k]),decreasing=TRUE)
x1<- as.matrix(X[ord1[1:nrow(X)],])
heatmap(x1,col=gcol2)
heatmap.2(x1, col=gcol2)

pheatmap(x1,col = gcol2,
          # clustering
          #distfun = dist, 
          #hclust = clus,
          # scaling (genes are in rows)
          #scale = "row",
          # color
          #col = redblackgreen, 
         Rowv=FALSE, 
         # labels
          #labRow = "", 
          #ColSideColors = class_labels, 
          # tweaking
          #trace = "none",
          fontsize_row = 2,
          #density.info = "none")
          )

# PVClust analysis using bootstrapping analysis:
library(pvclust)
result3 <- pvclust(test, nboot = 1000)
plot(result3, cex = 0.35, print.pv = FALSE)
pvrect(result3, alpha = 0.95)

##Site 2 for PCA 
#pca <- prcomp(t(gene_matrix))
pca <- prcomp(t(test)) # To see how genes may be clustered (genes become rows)
pca <- prcomp(test) #To see how Patients/Stage/Grades Cluster
#plot(pca$x[,1:2])
ipsgene <- c(c_gene, mcf7ips_overlap)
ipsgene <- unique(ipsgene)

#Adding back parameter of interest
pca$x <- rownames_to_column(as.data.frame(pca$x))
names(pca$x)[names(pca$x) == 'rowname'] <- 'PATIENT_ID'
pca$x <- merge(pca$x, extras, by="PATIENT_ID")

ggplot(as.data.frame(pca$x[,1:10])) + 
  geom_point(aes(x=PC1, y=PC2, color = pca$x$METASTASIS), size = 2) #+ 
  #geom_text(aes(x=PC1, y=PC2, label = rownames(pca$x)), size = 4, hjust=-0.25, vjust=-0.5, color = ifelse(rownames(pca$x) %in% ipsgene, "red", "black"))

#Notes:

#  1.) In general, matrices of gene data are usually samples in columns
#and genes in rows, which is the transpose of what prcomp() expects, so you
#have to use t().

#2.) Usually when I plot the results, I also use pch, col, xlab, ylab,
#main, etc. to make the plotting symbols for each group different
#shapes and colors, add reasonable axis labels, a main title, etc. See ?par
#for other variables you can pass to plot.

#3.) It is nice to follow up with legend() to add a nice legend to the
#plot. People seem to like that stuff ;-D.

#4.) The general recommendation is to set scale. = TRUE in the call to
#prcomp. I'm not sure if it will make much difference with microarray
#data because it is usually normalized anyway, so I usually don't
#bother.
#However, it is worth a try.














#Seurat testing for PCA:
pbmc <- CreateSeuratObject(raw.data = test, min.cells = 3, min.genes = 200, 
                           project = "Metabric")

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)


#### PCA Analysis
X = t(scale(test,center=TRUE,scale=FALSE)) #One has to be aware that in the microarray data, columns are usually samples(observations n), 
#rows are genes(features p). We need to first transpose X for the microarray data for the svd analysis 
# I HAVE REMOVED THE INTERIOR t FUNCTION, WHICH WAS: X = t(scale(t(test),center=TRUE,scale=FALSE))
sv = t(X)
w <- which(is.na(as.matrix(sv)))
#sv <- sv[which(!is.finite(as.matrix(sv)))] <- 0
sv[w] <- 0
sv = svd(sv)
U = sv$u
V = sv$v
D = sv$d

cols = as.numeric(as.factor(colnames(test)))
plot(U[,1],U[,2],type="n",xlab="PC1",ylab="PC2")
text(U[,1],U[,2],colnames(X),col=cols)

par(mfrow=c(1,1))
Z = t(X)%*%V

# plot PC1 vs PC2
plot(Z[,1], Z[,2], type ="n", xlab="PC1", ylab="PC2")
text(Z[,1], Z[,2], colnames(X), col=cols)

plot(Z[,2], Z[,3], type ="n", xlab = "PC2", ylab="PC3")
text(Z[,2], Z[,3], colnames(X), col=cols)

pc_dat<- data.frame(type = rownames(Z), PC1 = Z[,1], PC2= Z[,2])
library(ggplot2)
ggplot(pc_dat,aes(x=PC1, y=PC2, col=type)) + geom_point() + geom_text(aes(label = type), hjust=0, vjust=0)


######
# cBio - Metabric 
#####
#### READ IN PATIENT DATA AND TIDY
tumordata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/brca_metabric/data_clinical_supp_sample.csv", sep = ",", stringsAsFactors = FALSE)
patientdata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/brca_metabric/data_clinical_supp_patient.csv", sep = ",", stringsAsFactors = FALSE)

patientDF <- merge(x = tumordata, y = patientdata, by = "PATIENT_ID", all = TRUE)

patientDF <- subset(patientDF, select = c(PATIENT_ID, GRADE, TUMOR_STAGE)) 

#### READ IN CLINICAL DATA AND TIDY
expressiondata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/brca_metabric/data_expression.csv", sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

expressiondata[expressiondata == "null"] <- NA

expressiondata <- subset(expressiondata, select = -Entrez_Gene_Id)

expressiondata <- expressiondata %>% drop_na()

expressionDF <- melt(expressiondata, id.vars = "Hugo_Symbol")

colnames(expressionDF) <- c("GENE", "PATIENT_ID", "EXPRESSION_LEVEL")

#### MERGING CLINICAL AND PATIENT DATA 
patientDF <- merge(x = expressionDF, y = patientDF, by = "PATIENT_ID")

#### BRINGING IN GENE CRITERIA WE WANT TO FILTER OUR DATASETS
#Total iPS overlaps (including somatic targets - 93 genes)
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/HUViPSvsFiPS/")
huvfips_overlap <- data.frame(read.csv(paste(file.path, "HUViPSFiPSOverlap.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
huvfips_overlap <- as.factor(huvfips_overlap$x)

#All Human Cell surface-specific targets indicated by GO analysis on GO site 
data1 <- read.table("~/Desktop/HumanGOCellSurface.txt", sep = "" , header = F , nrows = 2000,
                    na.strings ="", stringsAsFactors= F, fill = T) #Download genes from text file that were created by GO analysis of proteins 
#that included "Cell Surface" as part of its GO

write.csv(data1, "data.csv") #Write surface data to a file making it a CSV

data <- read.csv("~/Desktop/data.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) # Import CSV of surface expression
surface_genes <- data$V1 
surface_genes <- unique(surface_genes) 

#Unique iPS targets(no somatic source - 34 genes)
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/CSVs/")
gene <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

#iPS to MCF7 Comparisons (No Somatic Filters - can include somatic hits)
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsHUViPS/")
mcf7huvips_overlap <- data.frame(read.csv(paste(file.path, "HUViPSMCF7OverlapGENE.csv", sep=""), header = FALSE, sep = ",", stringsAsFactors = FALSE))

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsFiPS/")
mcf7fips_overlap <- data.frame(read.csv(paste(file.path, "FiPSMCF7OverlapGENE.csv", sep=""), header = FALSE, sep = ",", stringsAsFactors = FALSE))

mcf7ips_overlap <- intersect(mcf7fips_overlap$V1, mcf7huvips_overlap$V1)

#A smaller subset of genes pulled from FiPS and HUViPS MS sets to test clustering
testgenes <- read.csv("~/Desktop/TestGeneSetUnique.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) # Import CSV of surface expression
testgenes <- as.factor(testgenes$x)
c_testgenes <- as.character(testgenes)
c_gene <- as.character(gene)
totaltestgenes <- c(c_testgenes, c_gene)
totaltestgenes <- unique(totaltestgenes)
totaltestgenes <- as.factor(totaltestgenes)

#### FILTERING DATA WITH OUR CRITERIA OF INTEREST
grades <- factor(c("1", "2", "3"))
grades <- ordered(grades, levels = c("1", "2", "3"))
stages <- factor(c("1", "2", "3", "4"))
stages <- ordered(stages, levels = c("1", "2", "3", "4"))
stage_diffs <- subset(patientDF, GENE %in% totaltestgenes)
#stage_diffs <- subset(patientDF, GENE %in% huvfips_overlap) #Subset rows that are matches to those in the targets
#stage_diffs <- subset(stage_diffs, GENE %in% gene)
stage_diffs <- subset(stage_diffs, GRADE %in% grades)
stage_diffs <- subset(stage_diffs, TUMOR_STAGE %in% stages)

#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 1)
#stage_diffs <- filter(stage_diffs, ((GRADE == 1) | (GRADE == 2)))
#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 2)
#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 3)
stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 4)
#stage_diffs <- filter(stage_diffs, GRADE == 3)

patient <- dcast(stage_diffs, PATIENT_ID+GRADE+TUMOR_STAGE ~ GENE, value.var = "EXPRESSION_LEVEL")

extras <- cbind.data.frame(patient$PATIENT_ID, patient$GRADE, patient$TUMOR_STAGE)
colnames(extras) <- c("PATIENT_ID", "TUMOR", "STAGES")
#
splits <- str_split_fixed(patient$PATIENT_ID, pattern = "-", "2")

patient <- cbind.data.frame(splits, patient)

patient <- subset(patient, select = -c(`1`, PATIENT_ID)) # Getting rid of Stage Column to prepare for a calculation matrix

names(patient)[names(patient) == '2'] <- 'ID'

patient <- arrange(patient, TUMOR_STAGE, GRADE)

patient <- patient %>% 
  unite(PATIENT_ID, ID, TUMOR_STAGE, GRADE, sep = "_", remove = TRUE)

#
patient <- subset(patient, select = -c(GRADE, TUMOR_STAGE))
patient[,1] <- as.character(patient[,1])

patient[,2:ncol(patient)] <- sapply(patient[,2:ncol(patient)],as.numeric)

sapply(patient, class)  #to check classes

patient <- patient %>% remove_rownames %>% column_to_rownames(var="PATIENT_ID") # Making the patient IDs the rownames, not the first column

#patientDF <- t(patientDF)

redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)  # Making a color palette for the heatmap 

#Setting up the data for a heatmap/cluster analysis
patient <- patient[, apply(patient, 2, sum)!=0] ### Required to remove all the columns with 0 in them to get a distance matrix
test <- as.matrix(patient)
test1 <- t(test)
#heatmap(test1, Colv=NA, col=greenred(10),scale="none")
#cor(t(test1))
#t <- 1-cor(t(test1))
#t <- as.dist(1-cor(t(test1)))
hc <- hclust(as.dist(1-cor(t(test1))))
plot(hc)
heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=greenred(10), cexRow = 0.2)
heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=redblackgreen, cexRow = 0.2)

heatmap.2(test, 
          #distfun = dist_cor, hclustfun = clus_wd2, 
          # clustering
          #distfun = dist, 
          #hclust = clus,
          # scaling (genes are in rows)
          scale = "row",
          # color
          Rowv = FALSE,
          col = redblackgreen, 
          cexRow = 0.15,
          cexCol = 0.5,
          # labels
          #ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none")

# PVClust analysis using bootstrapping analysis:
library(pvclust)
result1 <- pvclust(test, nboot = 1000)
plot(result1, cex = 0.35, print.pv = FALSE)
pvrect(result1, alpha = 0.95)

##Site 2 for PCA 
#pca <- prcomp(t(gene_matrix))
pca <- prcomp(t(test)) # To see how genes may be clustered
pca <- prcomp(test) #To see how Patients/Stage/Grades Cluster
#plot(pca$x[,1:2])
ipsgene <- c(c_gene, mcf7ips_overlap)
ipsgene <- unique(ipsgene)

#Adding back parameter of interest
pca$x <- rownames_to_column(as.data.frame(pca$x))
names(pca$x)[names(pca$x) == 'rowname'] <- 'PATIENT_ID'
pca$x <- merge(pca$x, extras, by="PATIENT_ID")

ggplot(as.data.frame(pca$x[,1:10])) + 
  geom_point(aes(x=PC1, y=PC3, color = pca$x$TUMOR), size = 2) #+ 
#geom_text(aes(x=PC1, y=PC2, label = rownames(pca$x)), size = 4, hjust=-0.25, vjust=-0.5, color = ifelse(rownames(pca$x) %in% ipsgene, "red", "black"))

#Notes:

#  1.) In general, matrices of gene data are usually samples in columns
#and genes in rows, which is the transpose of what prcomp() expects, so you
#have to use t().

#2.) Usually when I plot the results, I also use pch, col, xlab, ylab,
#main, etc. to make the plotting symbols for each group different
#shapes and colors, add reasonable axis labels, a main title, etc. See ?par
#for other variables you can pass to plot.

#3.) It is nice to follow up with legend() to add a nice legend to the
#plot. People seem to like that stuff ;-D.

#4.) The general recommendation is to set scale. = TRUE in the call to
#prcomp. I'm not sure if it will make much difference with microarray
#data because it is usually normalized anyway, so I usually don't
#bother.
#However, it is worth a try.



#### PCA Analysis
X = t(scale(test,center=TRUE,scale=FALSE)) #One has to be aware that in the microarray data, columns are usually samples(observations n), 
#rows are genes(features p). We need to first transpose X for the microarray data for the svd analysis 
# I HAVE REMOVED THE INTERIOR t FUNCTION, WHICH WAS: X = t(scale(t(test),center=TRUE,scale=FALSE))
sv = t(X)
w <- which(is.na(as.matrix(sv)))
#sv <- sv[which(!is.finite(as.matrix(sv)))] <- 0
sv[w] <- 0
sv = svd(sv)
U = sv$u
V = sv$v
D = sv$d

cols = as.numeric(as.factor(colnames(test)))
plot(U[,1],U[,2],type="n",xlab="PC1",ylab="PC2")
text(U[,1],U[,2],colnames(X),col=cols)

par(mfrow=c(1,1))
Z = t(X)%*%V

# plot PC1 vs PC2
plot(Z[,1], Z[,2], type ="n", xlab="PC1", ylab="PC2")
text(Z[,1], Z[,2], colnames(X), col=cols)

plot(Z[,2], Z[,3], type ="n", xlab = "PC2", ylab="PC3")
text(Z[,2], Z[,3], colnames(X), col=cols)

pc_dat<- data.frame(type = rownames(Z), PC1 = Z[,1], PC2= Z[,2])
library(ggplot2)
ggplot(pc_dat,aes(x=PC1, y=PC2, col=type)) + geom_point() + geom_text(aes(label = type), hjust=0, vjust=0)

aa<- grep("grey",colors())
bb<- grep("green",colors())
cc<-  grep("red",colors())
gcol2<- colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

## use the genes that drive the first PC1. This is the first major patter in the data
k=1
ord1<- order(abs(V[,k]),decreasing=TRUE)
x1<- as.matrix(X[ord1[1:nrow(X)],])
heatmap(x1,col=gcol2, cexRow = 0.1, cexCol = 0.5 )
heatmap.2(x1, col=gcol2)

pheatmap(x1,col = gcol2,
         # clustering
         #distfun = dist, 
         #hclust = clus,
         # scaling (genes are in rows)
         #scale = "row",
         # color
         #col = redblackgreen, 
         Rowv=FALSE, 
         # labels
         #labRow = "", 
         #ColSideColors = class_labels, 
         # tweaking
         #trace = "none",
         fontsize_row = 1,
         #density.info = "none")
)


# Examine stage differences between our total unique 34 iPS targets' expression

stage_diffs <- t(patientDF)
stage_diffs <- subset(stage_diffs, rownames(stage_diffs) %in% gene) #Subset rows that are matches to those in the iPS targets
stage_diffs <- t(stage_diffs)
stage_diffs <- as.data.frame(stage_diffs)
stage_diffs <- rownames_to_column(stage_diffs, "ID")
stage_diffs <- merge(x = stage_diffs, y = patient_match, by = "ID")
stage_diffs <- stage_diffs %>%
                  select(ID, Stage, everything())
d <- melt(stage_diffs)#, id = 1:2, measure = 3:ncol(stage_diffs))

ggplot(d) +
  geom_tile(aes(x = Stage, y = ID, fill = value)) +
  scale_fill_gradient(low = "black", high = "red", name = "") +
  facet_wrap(~ variable)


# Examine stage differences between our MCF7/iPS overlap without a somatic filter (29 genes)
#iPS to MCF7 Comparisons (No Somatic Filters - can include somatic hits)
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsHUViPS/")
mcf7huvips_overlap <- data.frame(read.csv(paste(file.path, "HUViPSMCF7OverlapGENE.csv", sep=""), header = FALSE, sep = ",", stringsAsFactors = FALSE))

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsFiPS/")
mcf7fips_overlap <- data.frame(read.csv(paste(file.path, "FiPSMCF7OverlapGENE.csv", sep=""), header = FALSE, sep = ",", stringsAsFactors = FALSE))

mcf7ips_overlap <- intersect(mcf7fips_overlap$V1, mcf7huvips_overlap$V1)
mcf7ips_overlap <- mcf7fips_overlap$V1

stage_diffs <- t(patientDF)
stage_diffs <- subset(stage_diffs, rownames(stage_diffs) %in% mcf7ips_overlap) #Subset rows that are matches to those in the iPS targets
stage_diffs <- t(stage_diffs)
stage_diffs <- as.data.frame(stage_diffs)
stage_diffs <- rownames_to_column(stage_diffs, "ID")
stage_diffs <- merge(x = stage_diffs, y = patient_match, by = "ID")
stage_diffs <- stage_diffs %>%
  select(ID, Stage, everything())
d <- melt(stage_diffs)#, id = 1:2, measure = 3:ncol(stage_diffs))

ggplot(d) +
  geom_tile(aes(x = Stage, y = ID, fill = value)) +
  scale_fill_gradient(low = "black", high = "red", name = "") +
  facet_wrap(~ variable)

ggplot(d) + 
  geom_point(aes(x=Stage, y=value), position = position_jitter()) + 
  geom_boxplot(aes(x = Stage, y = value), fill = "yellow", color = "darkblue") + 
  facet_wrap(~variable)               


##########

set <- expressiondata[, ALL$mol.biol %in% c("BCR/ABL", "ALL1/AF4")]

heatmapdata = t(scale(t(expressiondata),center=TRUE,scale=FALSE))

X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
heatmapdata <- t(expressiondata)
heatmapdata <- t(scale(heatmapdata, center = TRUE, scale = TRUE))

colnames(expressionDF) <- c("GENE", "PATIENT_ID", "EXPRESSION_LEVEL")

#Merging Data frames
patientDF <- merge(x = expressionDF, y = patientdata, by = "PATIENT_ID") # all = TRUE)

patientDF <- arrange(patientDF, PATIENT_ID, `AJCC Stage`, Tumor, GENE)

######Hierarchial Clustering Analysis
zzz <- patientDF
zzz <- zzz[,-(1:2)]
colnames(zzz) <- NULL
#ncidat1 = t(NCI60$data)
ncidat <- t(zzz)
#colnames(ncidat1) = NCI60$labs
#unique(colnames(ncidat1))
#dim(ncidat)

#colnames(ncidat) <- lapply(ncidat[1,], as.character)
colnames(ncidat) <- expressiondata$Hugo_Symbol
#ncidat <- ncidat[-1,]  
#ncidat <- ncidat[,-1]  
X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
sv = t(X)
#is.na(X) <- do.call(cbind,lapply(X, is.infinite))
#length(which(!is.finite(as.matrix(X))))
w <- which(is.na(as.matrix(sv)))
#sv <- sv[which(!is.finite(as.matrix(sv)))] <- 0
sv[w] <- 0
sv = svd(sv)
U = sv$u
V = sv$v
D = sv$d


aa<- grep("grey",colors())
bb<- grep("green",colors())
cc<-  grep("red",colors())
gcol2<- colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

## use the genes that drive the first PC1. This is the first major patter in the data
k=1
ord1<- order(abs(V[,k]),decreasing=TRUE)
x1<- as.matrix(X[ord1[1:nrow(X)],])
heatmap(x1,col=gcol2)


expression$tumor_grade <- patientDF$GRADE[match(expression$variable, patientDF$PATIENT_ID)]
expression$stage <- patientDF$TUMOR_STAGE[match(expression$variable, patientDF$PATIENT_ID)]
#expression$value <- nan.na(ips.hits$value)

x <- subset(expression, select = -variable)

y <- dcast(x, Hugo_Symbol ~ stage, mean)
z <- dcast(x, Hugo_Symbol ~ stage, median)
#write.csv(y, "test.csv")

dy <- dist(as.matrix(y))   # find distance matrix 
dz <- dist(as.matrix(z))

hcy <- hclust(dy)                # apply hirarchical clustering 
hcz <- hclust(dz)

plot(hcy)




file.path <- ("~/Desktop/")
gene <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

split.names <- str_split_fixed(proteinDF$Composite.Element.REF,"[|]", 2) #special character, so needed [ ]
proteinDF <- data.frame(proteinDF, split.names)
proteinDF <- subset(proteinDF, select = -c(Composite.Element.REF, `X2`))
names(proteinDF)[names(proteinDF) == 'X1'] <- 'Gene'

ips.targets <- proteinDF[proteinDF$Gene %in% gene,]
ips.targets <- melt(ips.targets)
ips.targets <- subset(ips.targets, select = -variable)
ips.targets$value <- nan.na(ips.targets$value)
ips.targets <- ips.targets %>%
  group_by(Gene) %>%
  mutate(mean = mean(value, na.rm = TRUE))

ips.hits <- rnaDF[rnaDF$Hugo_Symbol %in% gene,]
ips.hits <- subset(ips.hits, select = -Entrez_Gene_Id)
ips.hits <- melt(ips.hits)
ips.hits <- subset(ips.hits, select = -variable)
ips.hits$value <- nan.na(ips.hits$value)

ggplot(ips.hits) + 
  geom_point(aes(x=Hugo_Symbol, y=value)) + 
  ggtitle("Median Z-Scores (mRNA) vs Gene by Tumor Grade") + 
  xlab("Median Z-score") + 
  ylab("Gene") + 
  theme(axis.text.x = element_text(angle = 90))

ips.hits <- ips.hits %>%
  group_by(Hugo_Symbol) %>%
  mutate(mean = mean(value, na.rm = TRUE))




######Hierarchial Clustering Analysis
install.packages("ISLR")
library(ISLR)
zzz <- expressiondata
zzz <- zzz[,-(1:2)]
colnames(zzz) <- NULL
#ncidat1 = t(NCI60$data)
ncidat <- t(zzz)
#colnames(ncidat1) = NCI60$labs
#unique(colnames(ncidat1))
#dim(ncidat)

#colnames(ncidat) <- lapply(ncidat[1,], as.character)
colnames(ncidat) <- expressiondata$Hugo_Symbol
#ncidat <- ncidat[-1,]  
#ncidat <- ncidat[,-1]  
#X = t(scale(t(ncidat), center=TRUE, scale = FALSE))
X = t(scale(t(test4),center=TRUE,scale=FALSE))
sv = t(X)
#is.na(X) <- do.call(cbind,lapply(X, is.infinite))
#length(which(!is.finite(as.matrix(X))))
w <- which(is.na(as.matrix(sv)))
#sv <- sv[which(!is.finite(as.matrix(sv)))] <- 0
sv[w] <- 0
sv = svd(sv)
U = sv$u
V = sv$v
D = sv$d


aa<- grep("grey",colors())
bb<- grep("green",colors())
cc<-  grep("red",colors())
gcol2<- colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

## use the genes that drive the first PC1. This is the first major patter in the data
k=1
ord1<- order(abs(V[,k]),decreasing=TRUE)
x1<- as.matrix(X[ord1[1:nrow(X)],])
heatmap(x1,col=gcol2)


expression$tumor_grade <- patientDF$GRADE[match(expression$variable, patientDF$PATIENT_ID)]
expression$stage <- patientDF$TUMOR_STAGE[match(expression$variable, patientDF$PATIENT_ID)]
#expression$value <- nan.na(ips.hits$value)

x <- subset(expression, select = -variable)

y <- dcast(x, Hugo_Symbol ~ stage, mean)
z <- dcast(x, Hugo_Symbol ~ stage, median)
#write.csv(y, "test.csv")

dy <- dist(as.matrix(y))   # find distance matrix 
dz <- dist(as.matrix(z))

hcy <- hclust(dy)                # apply hirarchical clustering 
hcz <- hclust(dz)

plot(hcy)

install.packages("factoextra")
install.packages("cluster")
install.packages("magrittr")
#####
