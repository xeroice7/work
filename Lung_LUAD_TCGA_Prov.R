library(tidyverse)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ISLR)
library(gplots)
library("Biobase")

#######################################
# cBio - Lung Adenocarcinoma TCGA (Provisional) - 522 Patients
#######################################

### 
# RNA Seq Data (too few recordings for microarray expression)
###

#### READ IN PATIENT DATA AND TIDY
patientdata <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/data_bcr_clinical_data_patient.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

colnames(patientdata) <- lapply(patientdata[4,], as.character)

patientdata <- patientdata[-(1:4),]

patientdata <- subset(patientdata, select = c(PATIENT_ID,  AJCC_TUMOR_PATHOLOGIC_PT, AJCC_METASTASIS_PATHOLOGIC_PM, AJCC_PATHOLOGIC_TUMOR_STAGE, OS_STATUS)) 

splits <- str_split_fixed(patientdata$AJCC_PATHOLOGIC_TUMOR_STAGE, " ", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '2'] <- 'Specific_Stage'
patientdata <- subset(patientdata, select=-c(1, AJCC_PATHOLOGIC_TUMOR_STAGE))

splits <- str_split_fixed(patientdata$Specific_Stage, "A|B|C", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'General_Stages'
patientdata <- subset(patientdata, select=-2)

#### READ IN CLINICAL DATA AND TIDY (cant use microarray - too few patients recorded)
expressiondata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/data_RNA_Seq_v2_expression_median_edited.csv", check.names=FALSE, sep = ",", stringsAsFactors = FALSE)

expressiondata[expressiondata == "null"] <- NA

expressiondata <- subset(expressiondata, select = -Entrez_Gene_Id)

expressiondata <- expressiondata %>% drop_na()

###### To check to see if there are duplicate patients or genes
e <- t(expressiondata)
splits <- str_split_fixed(rownames(e), "-", 4)
e <- cbind.data.frame(splits, e)
e <- e %>% 
  unite(ID, '1', '2', '3', sep = "-", remove = TRUE)
any(duplicated(e$ID))
which(duplicated(e$ID))
any(duplicated(expressiondata$Hugo_Symbol)) ## Duplicate Genes?
which(duplicated(expressiondata$Hugo_Symbol)) 
#Replace any genes duplicates with "Gene-1"
doubles <- which(duplicated(expressiondata$Hugo_Symbol)) 
names <- expressiondata$Hugo_Symbol[doubles]
new_names <- paste(names, "-1", sep="")
expressiondata$Hugo_Symbol[doubles] <- new_names
#######

expressionDF <- melt(expressiondata, id.vars = "Hugo_Symbol")

splits <- str_split_fixed(expressionDF$variable, "-", 4)
expressionDF <- cbind.data.frame(splits, expressionDF)
expressionDF <- expressionDF %>% 
  unite(ID, '1', '2', '3', sep = "-", remove = TRUE)
expressionDF <- subset(expressionDF, select=-c(`4`, variable))

colnames(expressionDF) <- c("PATIENT_ID", "GENE", "EXPRESSION_LEVEL")

#### MERGING CLINICAL AND PATIENT DATA 
patientDF <- merge(x = expressionDF, y = patientdata, by = "PATIENT_ID")

#### BRINGING IN GENE CRITERIA WE WANT TO FILTER OUR DATASETS
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
grades <- factor(c("T1", "T1a", "T1b", "T2", "T2a", "T2b", "T3"))
grades <- ordered(grades, levels = c("T1", "T1a", "T1b", "T2", "T2a", "T2b", "T3"))
stages <- factor(c("I", "IA", "IB", "II", "IIA", "IIB", "IIIA", "IIIB", "IV"))
stages <- ordered(stages, levels = c("I", "IA", "IB", "II", "IIA", "IIB", "IIIA", "IIIB", "IV"))

stage_diffs <- patientDF
stage_diffs <- subset(stage_diffs, GENE %in% huvfips_overlap) #Subset rows that are matches to those in the targets
#stage_diffs <- subset(stage_diffs, GENE %in% surface_genes)
#stage_diffs <- subset(stage_diffs, AJCC_METASTASIS_PATHOLOGIC_PM %in% meta)
#stage_diffs <- subset(stage_diffs, AJCC_TUMOR_PATHOLOGIC_PT %in% grades)
#stage_diffs <- subset(stage_diffs, Specific_Stage %in% stages)
#stage_diffs <- filter(stage_diffs, ((Tumor == "T3") | (Tumor == "T2")))
stage_diffs <- filter(stage_diffs, AJCC_TUMOR_PATHOLOGIC_PT == "T1")
#stage_diffs <- filter(stage_diffs, General_Stages == "I")
#stage_diffs <- filter(stage_diffs, AJCC_METASTASIS_PATHOLOGIC_PM == "M1")


patient <- dcast(stage_diffs, PATIENT_ID+AJCC_TUMOR_PATHOLOGIC_PT+Specific_Stage+AJCC_METASTASIS_PATHOLOGIC_PM ~ GENE, value.var = "EXPRESSION_LEVEL")

splits <- str_split_fixed(patient$PATIENT_ID, pattern = "-", "3")

patient <- cbind.data.frame(splits, patient)

patient <- subset(patient, select = -c(`1`, `2`, PATIENT_ID)) # Getting rid of Stage Column to prepare for a calculation matrix

names(patient)[names(patient) == '3'] <- 'ID'
patient <- arrange(patient, Specific_Stage, AJCC_TUMOR_PATHOLOGIC_PT, AJCC_METASTASIS_PATHOLOGIC_PM)

patient <- patient %>% 
  unite(PATIENT_ID, ID, Specific_Stage, AJCC_TUMOR_PATHOLOGIC_PT, AJCC_METASTASIS_PATHOLOGIC_PM, sep = "_", remove = TRUE)

#patient_match <- data.frame(patientDF$ID, patientDF$Stage)  # Creating a data frame to match patients with stages later

#colnames(patient_match) <- c("ID", "Stage")
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
grade3 <- test1
#heatmap(test1, Colv=NA, col=greenred(10),scale="none")
#cor(t(test1))
#t <- 1-cor(t(test1))
x <- cor(t(test1))
y <- cor(t(grade1), t(grade3))
#hc <- hclust(as.dist(1-cor(t(test1))))
#plot(hc)
#heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=greenred(10), cexRow = 0.2)
#heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=redblackgreen, cexRow = 0.2)
heatmap.2(y, 
          Rowv=T, 
          Colv=T, 
          col=bluered(256), 
          #breaks = breaks,
          cexRow = 0.75, #0.2
          cexCol = 0.75, 
          scale = "none", 
          trace = "none" 
          #ColSideColors = ifelse(rownames(x) %in% ipsgene, "red", "black"),
          #RowSideColors = ifelse(rownames(x) %in% ipsgene, "red", "black")
)
