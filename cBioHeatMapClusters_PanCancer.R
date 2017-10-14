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

d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit) # display dendogram

######
# cBio - Cancer Cell Line Encyclopedia
#####
#### READ IN PATIENT DATA AND TIDY
file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/cellline_ccle_broad/cellline_ccle_broad/")
patientdata <- data.frame(read.csv(paste(file.path, "data_clinical.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

colnames(patientdata) <- lapply(patientdata[5,], as.character)

patientdata <- patientdata[-(1:5),]

patientdata <- subset(patientdata, select = c(SAMPLE_ID, PATIENT_ID, HISTOLOGICAL_SUBTYPE, TUMOR_TYPE)) 

#### READ IN CLINICAL DATA AND TIDY
file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/cellline_ccle_broad/cellline_ccle_broad/")
expressiondata <- read.csv(paste(file.path, "data_expression_median.csv", sep=""), check.names=FALSE, sep = ",", stringsAsFactors = FALSE)

expressiondata[expressiondata == "null"] <- NA

expressiondata <- subset(expressiondata, select = -Entrez_Gene_Id)

expressiondata <- expressiondata %>% drop_na()

expressionDF <- melt(expressiondata, id.vars = "Hugo_Symbol")

colnames(expressionDF) <- c("GENE", "SAMPLE_ID", "EXPRESSION_LEVEL")


#### MERGING CLINICAL AND PATIENT DATA 
patientDF <- merge(x = expressionDF, y = patientdata, by = "SAMPLE_ID")
patientDF$HISTOLOGICAL_SUBTYPE <- as.factor(patientDF$HISTOLOGICAL_SUBTYPE)
patientDF$TUMOR_TYPE <- as.factor(patientDF$TUMOR_TYPE)

#### BRINGING IN GENE CRITERIA WE WANT TO FILTER OUR DATASETS
#Total iPS overlaps (including somatic targets - 93 genes)
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/HUViPSvsFiPS/")
huvfips_overlap <- data.frame(read.csv(paste(file.path, "HUViPSFiPSOverlap.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
huvfips_overlap <- as.factor(huvfips_overlap$x)

#All Human Cell surface-specific targets indicated by GO analysis on GO site 
#Download genes from text file that were created by GO analysis of proteins 
#that included "Cell Surface" as part of its GO

data <- read.csv("~/Desktop/data.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) # Import CSV of surface expression
surface_genes <- data$V1 
surface_genes <- unique(surface_genes) 

#Unique iPS targets(no somatic source - 34 genes)
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/CSVs/")
gene <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

#### FILTERING DATA WITH OUR CRITERIA OF INTEREST
cancer_types <- c("breast", "lung")
stage_diffs <- subset(patientDF, GENE %in% huvfips_overlap) #Subset rows that are matches to those in the targets
stage_diffs <- subset(stage_diffs, GENE %in% surface_genes)
stage_diffs <- subset(stage_diffs, TUMOR_TYPE %in% cancer_types)

patient <- dcast(stage_diffs, SAMPLE_ID+PATIENT_ID+HISTOLOGICAL_SUBTYPE+TUMOR_TYPE ~ GENE, value.var = "EXPRESSION_LEVEL")

patient <- subset(patient, select = -c(HISTOLOGICAL_SUBTYPE, SAMPLE_ID)) # Getting rid of Stage Column to prepare for a calculation matrix

patient <- arrange(patient, TUMOR_TYPE)

patient <- patient %>% 
  unite(PATIENT_ID, PATIENT_ID, TUMOR_TYPE, sep = "_", remove = TRUE)

patient[,1] <- as.character(patient[,1])
patient[,2:ncol(patient)] <- sapply(patient[,2:ncol(patient)],as.numeric)

sapply(patient, class)  #to check classes

patient <- patient %>% remove_rownames %>% column_to_rownames(var="PATIENT_ID") # Making the patient IDs the rownames, not the first column

#patientDF <- t(patientDF)

redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)  # Making a color palette for the heatmap 


#Setting up the data for a heatmap/cluster analysis
test <- as.matrix(patient)

heatmap.2(test, 
          #distfun = dist_cor, hclustfun = clus_wd2, 
          # clustering
          #distfun = dist, 
          #hclust = clus,
          # scaling (genes are in rows)
          scale = "row",
          # color
          #Rowv = FALSE,
          col = redblackgreen, 
          cexRow = 0.075,
          cexCol = 0.5,
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
heatmap(x1,col=gcol2, cexRow = 0.1, cexCol = 0.5)
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