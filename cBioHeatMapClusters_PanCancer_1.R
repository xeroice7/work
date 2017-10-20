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

patientdata <- patientdata %>% drop_na(TUMOR_TYPE)

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

#MCF7 + iPS Overlap
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsHUViPS/")
mcf7huvips_overlap <- data.frame(read.csv(paste(file.path, "HUViPSMCF7OverlapGENE.csv", sep=""), header = FALSE, sep = ",", stringsAsFactors = FALSE))

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsFiPS/")
mcf7fips_overlap <- data.frame(read.csv(paste(file.path, "FiPSMCF7OverlapGENE.csv", sep=""), header = FALSE, sep = ",", stringsAsFactors = FALSE))

mcf7ips_overlap <- intersect(mcf7fips_overlap$V1, mcf7huvips_overlap$V1)

#Unique iPS targets(no somatic source - 34 genes)
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/CSVs/")
gene <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

#### FILTERING DATA WITH OUR CRITERIA OF INTEREST
#cancer_types <- c("breast", "lung")
#stage_diffs <- subset(patientDF, GENE %in% huvfips_overlap) #Subset rows that are matches to those in the targets
#stage_diffs <- subset(stage_diffs, GENE %in% surface_genes)
#stage_diffs <- subset(stage_diffs, TUMOR_TYPE %in% cancer_types)

#ONLY TO BE DONE WITH ALL GENES FROM A SET
stage_diffs <- patientDF

stage_diffs$Number <- 1:nrow(stage_diffs)

stage_diffs <- stage_diffs %>% 
  unite(PATIENT_ID, PATIENT_ID, TUMOR_TYPE, sep = "_", remove = TRUE)

stage_diffs <- subset(stage_diffs, select=-c(HISTOLOGICAL_SUBTYPE, SAMPLE_ID))

stage_diffs <- arrange(stage_diffs, PATIENT_ID, GENE)

which(stage_diffs$GENE == "1-Dec")

stage_diffs <- stage_diffs[1:1605681,] #subsetting partial data when a new patient begins

patient <- dcast(stage_diffs, PATIENT_ID+Number ~ GENE, value.var = "EXPRESSION_LEVEL")

##### RESUME NORMAL CODE
patient <- dcast(stage_diffs, SAMPLE_ID+PATIENT_ID+HISTOLOGICAL_SUBTYPE+TUMOR_TYPE ~ GENE, value.var = "EXPRESSION_LEVEL")

patient <- subset(patient, select = -c(HISTOLOGICAL_SUBTYPE, SAMPLE_ID)) # Getting rid of Stage Column to prepare for a calculation matrix

patient <- arrange(patient, TUMOR_TYPE)

patient$cancertypes <- as.character(patient$TUMOR_TYPE)

for (i in 1:length(patient$cancertypes)) {
  if (patient$cancertypes[i] == "glioma" | patient$cancertypes[i] =="medulloblastoma" | patient$cancertypes[i] =="neuroblastoma") {
    patient$cancertypes[i] <- "Brain"
  } else if (patient$cancertypes[i] == "AML" | patient$cancertypes[i] == "B-cell_ALL" | patient$cancertypes[i] == "CML" | patient$cancertypes[i] == "leukemia_other" | patient$cancertypes[i] == "multiple_myeloma" | patient$cancertypes[i] == "T-cell_ALL") {
    patient$cancertypes[i] <- "Leukemia"
  } else if (patient$cancertypes[i] == "breast") {
    patient$cancertypes[i] <- "Breast"
  } else if (patient$cancertypes[i] == "lung_NSC" | patient$cancertypes[i] == "lung_small_cell") {
    patient$cancertypes[i] <- "Lung"
  } else if (patient$cancertypes[i] == "lymphoma_Burkitt" | patient$cancertypes[i] == "lymphoma_DLBCL" | patient$cancertypes[i] == "lymphoma_Hodgkin" | patient$cancertypes[i] == "lymphoma_other") {
    patient$cancertypes[i] <- "Lymphoma"
  } else if (patient$cancertypes[i] == "chondrosarcoma" | patient$cancertypes[i] == "Ewings_Sarcoma" | patient$cancertypes[i] ==  "osteosarcoma") {
    patient$cancertypes[i] <- "Bone"
  } else if (patient$cancertypes[i] == "pancreas") {
    patient$cancertypes[i] <- "Pancreas"
  } else if (patient$cancertypes[i] == "colorectal") {
    patient$cancertypes[i] <- "Colorectal"
  } else {
    patient$cancertypes[i] <- "Other"
  }
  print(patient$cancertypes[i])
} 

patient$cancertypes <- as.factor(patient$cancertypes)

patient <- patient %>% 
  unite(PATIENT_ID, PATIENT_ID, cancertypes, sep = "_", remove = TRUE)

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

# PVClust analysis using bootstrapping analysis:
#library(pvclust)
result <- pvclust(test, nboot = 10000)
plot(result)
pvrect(result, alpha = 0.95)



###PCA Analysis
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