library(tidyverse)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ISLR)
library(gplots)
library("Biobase")

######
# cBio - TCGA Cell 2015
#####
#### READ IN PATIENT DATA AND TIDY
patientdata <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_tcga_pub/brca_tcga_pub/data_clinical.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

colnames(patientdata) <- lapply(patientdata[5,], as.character)

patientdata <- patientdata[-(1:5),]

patientdata <- subset(patientdata, select = c(PATIENT_ID, Tumor, `AJCC Stage`, Metastasis)) 

splits <- str_split_fixed(patientdata$`AJCC Stage`, " ", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '2'] <- 'Specific_Stages'
patientdata <- subset(patientdata, select=-1)

splits <- str_split_fixed(patientdata$Specific_Stage, "A|B|C", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'General_Stages'
patientdata <- subset(patientdata, select=-c(2, `AJCC Stage`))

#### READ IN CLINICAL DATA AND TIDY
expressiondata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_tcga_pub/brca_tcga_pub/data_expression_median.csv", sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

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

#Combination of Jun's list (29 genes) and my unique list (34 genes) - total of 50 genes
ipsgene <- c(c_gene, mcf7ips_overlap)
ipsgene <- unique(ipsgene)

#### FILTERING DATA WITH OUR CRITERIA OF INTEREST
meta <- factor(c("M0", "M1"))
meta <- ordered(meta, levels = c("M0", "M1"))
grades <- factor(c("T1", "T2", "T3"))
grades <- ordered(grades, levels = c("T1", "T2", "T3"))
stages <- factor(c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IIIC", "IV"))
stages <- ordered(stages, levels = c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IIIC", "IV"))
stage <- factor(c("I", "II", "III", "IV"))
stage <- ordered(stages, levels = c("I", "II", "III", "IV"))

stage_diffs <- patientDF
stage_diffs <- subset(stage_diffs, GENE %in% huvfips_overlap) #Subset rows that are matches to those in the targets
#stage_diffs <- subset(stage_diffs, GENE %in% ipsgene)
#stage_diffs <- subset(stage_diffs, GENE %in% surface_genes)
#stage_diffs <- subset(stage_diffs, Metastasis %in% meta)
#stage_diffs <- subset(stage_diffs, Tumor %in% grades)
#stage_diffs <- subset(stage_diffs, Specific_Stages %in% stages)
stage_diffs <- subset(stage_diffs, General_Stages %in% stage)
#stage_diffs <- filter(stage_diffs, ((Tumor == "T3") | (Tumor == "T2")))
#stage_diffs <- filter(stage_diffs, Tumor == "T1")
stage_diffs <- filter(stage_diffs, General_Stages == "I")
#stage_diffs <- filter(stage_diffs, Metastasis == "M1")

patient <- dcast(stage_diffs, PATIENT_ID+Tumor+General_Stages+Specific_Stages+Metastasis ~ GENE, value.var = "EXPRESSION_LEVEL")

#extras <- cbind.data.frame(patient$PATIENT_ID, patient$Metastasis, patient$Tumor, patient$Stages, patient$Stage)
#colnames(extras) <- c("PATIENT_ID", "METASTASIS", "TUMOR", "STAGES", "STAGE")
#
splits <- str_split_fixed(patient$PATIENT_ID, pattern = "-", "3")

patient <- cbind.data.frame(splits, patient)

patient <- subset(patient, select = -c(`1`, `2`, PATIENT_ID)) # Getting rid of Stage Column to prepare for a calculation matrix

names(patient)[names(patient) == '3'] <- 'ID'

patient <- arrange(patient, Specific_Stages, Tumor, Metastasis)

patient <- patient %>% 
  unite(PATIENT_ID, ID, Specific_Stages, Tumor, Metastasis, sep = "_", remove = TRUE)

patient <- subset(patient, select = -Specific_Stages)

#
patient <- subset(patient, select = -c(General_Stages, Tumor, Metastasis, Specific_Stages))
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
x <- cor(t(test1))
hc <- hclust(as.dist(1-cor(t(test1))))
#plot(hc)
#heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=greenred(10), cexRow = 0.2)
#heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=redblackgreen, cexRow = 0.2)
heatmap.2(x, 
          Rowv=T, 
          Colv=T, 
          col=bluered(256), 
          #breaks = breaks,
          cexRow = 0.75, #0.2
          cexCol = 0.75, 
          scale = "none", 
          trace = "none", 
          ColSideColors = ifelse(rownames(x) %in% ipsgene, "red", "black"),
          RowSideColors = ifelse(rownames(x) %in% ipsgene, "red", "black")
)


######
# cBio - Metabric 
#####
#### READ IN PATIENT DATA AND TIDY
tumordata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/data_clinical_supp_sample.csv", sep = ",", stringsAsFactors = FALSE)
patientdata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/data_clinical_supp_patient.csv", sep = ",", stringsAsFactors = FALSE)

patientDF <- merge(x = tumordata, y = patientdata, by = "PATIENT_ID", all = TRUE)

patientDF <- subset(patientDF, select = c(PATIENT_ID, GRADE, TUMOR_STAGE)) 

#### READ IN CLINICAL DATA AND TIDY
expressiondata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/data_expression.csv", sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

expressiondata[expressiondata == "null"] <- NA

expressiondata <- subset(expressiondata, select = -Entrez_Gene_Id)

expressiondata <- expressiondata %>% drop_na()

expressionDF <- melt(expressiondata, id.vars = "Hugo_Symbol")

colnames(expressionDF) <- c("GENE", "PATIENT_ID", "EXPRESSION_LEVEL")

#### MERGING CLINICAL AND PATIENT DATA 
patientDF <- merge(x = expressionDF, y = patientDF, by = "PATIENT_ID")

#### BRINGING IN GENE CRITERIA WE WANT TO FILTER OUR DATASETS
#Total iPS overlaps (including somatic targets - 93 genes)
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

#Combination of Jun's list (29 genes) and my unique list (34 genes) - total of 50 genes
ipsgene <- c(c_gene, mcf7ips_overlap)
ipsgene <- unique(ipsgene)

#### FILTERING DATA WITH OUR CRITERIA OF INTEREST
grades <- factor(c("1", "2", "3"))
grades <- ordered(grades, levels = c("1", "2", "3"))
stages <- factor(c("1", "2", "3", "4"))
stages <- ordered(stages, levels = c("1", "2", "3", "4"))

stage_diffs <- patientDF
#stage_diffs <- subset(patientDF, GENE %in% totaltestgenes)
stage_diffs <- subset(stage_diffs, GENE %in% mcf7ips_overlap) #Subset rows that are matches to those in the targets
#stage_diffs <- subset(stage_diffs, GENE %in% gene)
#stage_diffs <- subset(stage_diffs, GRADE %in% grades)
stage_diffs <- subset(stage_diffs, TUMOR_STAGE %in% stages)

stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 4)
#stage_diffs <- filter(stage_diffs, ((GRADE == 1) | (GRADE == 2)))
#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 2)
#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 3)
#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 4)
#stage_diffs <- filter(stage_diffs, GRADE == 1)

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
grade1 <- test1
grade4 <- test1
#heatmap(test1, Colv=NA, col=greenred(10),scale="none")
#cor(t(test1))
t <- 1-cor(t(test1))
x <- cor(t(test1))
y <- cor(t(grade1), t(grade4))
hc <- hclust(as.dist(1-cor(t(test1))))
#plot(hc)
#heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=greenred(10), cexRow = 0.2)
#heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=redblackgreen, cexRow = 0.2)
heatmap.2(y, 
          Rowv=F, 
          Colv=F, 
          col=bluered(256), 
          #breaks = breaks,
          cexRow = 0.75, #0.2
          cexCol = 0.75, 
          scale = "none", 
          trace = "none", 
          ColSideColors = ifelse(rownames(x) %in% ipsgene, "red", "black"),
          RowSideColors = ifelse(rownames(x) %in% ipsgene, "red", "black")
)
