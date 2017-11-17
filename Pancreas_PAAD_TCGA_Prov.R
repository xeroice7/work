library(tidyverse)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ISLR)
library(gplots)
library(pvclust)

stage_split1 <- function(x, y) {
  splits <- str_split_fixed(x$y, " ", 2)
  #x <- cbind.data.frame(splits, x)
  #names(x)[names(x) == '2'] <- 'Stages'
  #x <- subset(x, select=-c(1, y))
}

######
# Expression Data
#####
#### READ IN PATIENT DATA AND TIDY
patientdata <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Pancreas/paad/tcga/data_bcr_clinical_data_patient.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

colnames(patientdata) <- lapply(patientdata[4,], as.character)

patientdata <- patientdata[-(1:4),]

patientdata <- subset(patientdata, select = c(PATIENT_ID, AJCC_PATHOLOGIC_TUMOR_STAGE, AJCC_METASTASIS_PATHOLOGIC_PM, AJCC_TUMOR_PATHOLOGIC_PT, OS_STATUS)) 

splits <- str_split_fixed(patientdata$AJCC_PATHOLOGIC_TUMOR_STAGE, " ", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '2'] <- 'Specific_Stages'
patientdata <- subset(patientdata, select=-c(1, AJCC_PATHOLOGIC_TUMOR_STAGE))

splits <- str_split_fixed(patientdata$Specific_Stage, "A|B|C", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'General_Stages'
patientdata <- subset(patientdata, select=-2)


#### READ IN CLINICAL DATA AND TIDY
expressiondata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Pancreas/paad/tcga/data_RNA_Seq_v2_expression_median_edited.csv", sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
#Got rid of 3 duplicate patients and several duplicate genes in "edited" copy

expressiondata[expressiondata == "null"] <- NA

expressiondata <- subset(expressiondata, select = -Entrez_Gene_Id)

expressiondata <- expressiondata %>% drop_na()

###### To check to see if there are duplicate patients
e <- t(expressiondata) 
splits <- str_split_fixed(rownames(e), "-", 4)
e <- cbind.data.frame(splits, e)
e <- e %>% 
  unite(ID, '1', '2', '3', sep = "-", remove = TRUE)
any(duplicated(e$ID))
which(duplicated(e$ID))
any(duplicated(expressiondata$Hugo_Symbol))
which(duplicated(expressiondata$Hugo_Symbol))
#########

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
grades <- factor(c("T1", "T2", "T3", "T4"))
grades <- ordered(grades, levels = c("T1", "T2", "T3", "T4"))
spec.stages <- factor(c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IV"))
spec.stages <- ordered(spec.stages, levels = c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IV"))
gen.stages <- factor(c("I", "II", "III", "IV"))
gen.stages <- ordered(gen.stages, levels = c("I", "II", "III", "IV"))
status <- factor(c("LIVING", "DECEASED"))
status <- ordered(status, levels = c("LIVING", "DECEASED"))
meta <- factor(c("M0", "M1"))
meta <- ordered(meta, levels = c("M0", "M1"))

#Subset rows that are matches to those in the targets
stage_diffs <- patientDF
stage_diffs <- subset(stage_diffs, GENE %in% huvfips_overlap)
#stage_diffs <- subset(stage_diffs, OS_STATUS %in% status)
#stage_diffs <- subset(stage_diffs, AJCC_TUMOR_PATHOLOGIC_PT %in% grades)
#stage_diffs <- subset(stage_diffs, General_Stages %in% gen.stages)
#stage_diffs <- subset(stage_diffs, Specific_Stages %in% spec.stages)
#stage_diffs <- subset(stage_diffs, AJCC_METASTASIS_PATHOLOGIC_PM %in% meta)
#stage_diffs <- filter(stage_diffs, ((Tumor == "T3") | (Tumor == "T2")))
#stage_diffs <- filter(stage_diffs, General_Stages == "I")
#stage_diffs <- filter(stage_diffs, AJCC_TUMOR_PATHOLOGIC_PT == "T3")
stage_diffs <- filter(stage_diffs, AJCC_METASTASIS_PATHOLOGIC_PM == "M1")

patient <- dcast(stage_diffs, PATIENT_ID+General_Stages+Specific_Stages+AJCC_TUMOR_PATHOLOGIC_PT+AJCC_METASTASIS_PATHOLOGIC_PM+OS_STATUS ~ GENE, value.var = "EXPRESSION_LEVEL")
#if there is an error, there are duplicate genes in the array - check this with which(duplicated(expressiondata$Hugo_Symbol)) and change the names with X"-1"
extras <- cbind.data.frame(patient$PATIENT_ID, patient$Specific_Stages, patient$AJCC_METASTASIS_PATHOLOGIC_PM, patient$AJCC_TUMOR_PATHOLOGIC_PT, patient$General_Stages, patient$OS_STATUS)

colnames(extras) <- c("PATIENT_ID", "Specific_Stages", "Metastasis", "Grade", "General_Stages", "OS_STATUS")

splits <- str_split_fixed(patient$PATIENT_ID, pattern = "-", "3")

patient <- cbind.data.frame(splits, patient)

patient <- subset(patient, select = -c(`1`, `2`, PATIENT_ID)) # Getting rid of Stage Column to prepare for a calculation matrix

names(patient)[names(patient) == '3'] <- 'ID'

patient <- arrange(patient, General_Stages, Specific_Stages, AJCC_TUMOR_PATHOLOGIC_PT, AJCC_METASTASIS_PATHOLOGIC_PM, OS_STATUS)

patient <- patient %>% 
  unite(PATIENT_ID, ID, General_Stages, Specific_Stages, AJCC_TUMOR_PATHOLOGIC_PT, AJCC_METASTASIS_PATHOLOGIC_PM, sep = "_", remove = TRUE)

patient <- subset(patient, select = -OS_STATUS)

patient[,1] <- as.character(patient[,1])

patient[,2:ncol(patient)] <- sapply(patient[,2:ncol(patient)],as.numeric)

sapply(patient, class)  #to check classes

patient <- patient %>% remove_rownames %>% column_to_rownames(var="PATIENT_ID") # Making the patient IDs the rownames, not the first column

#patient <- patient[, sample(ncol(patient), 5000)] # Randomly sample n number of genes for analysis

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
hc <- hclust(as.dist(1-cor(t(test1))))
plot(hc)
heatmap(x, Rowv=as.dendrogram(hc), Colv=NA, col=greenred(10), cexRow = 0.2, scale = "column")
heatmap.2(y, 
          Rowv=NA, 
          Colv=NA, 
          col=bluered(256), 
          #breaks = breaks,
          cexRow = 0.75, 
          cexCol = 0.75, 
          scale = "none", 
          trace = "none"#, 
          #ColSideColors = ifelse(rownames(x) %in% ipsgene, "red", "black"),
          #RowSideColors = ifelse(rownames(x) %in% ipsgene, "red", "black")
)

# PVClust analysis using bootstrapping analysis:
result3 <- pvclust(test, nboot = 1000)
plot(result3, cex = 0.35, print.pv = FALSE)
pvrect(result3, alpha = 0.95)

##Site 2 for PCA 
#pca <- prcomp(t(gene_matrix))
pca <- prcomp(test1) # To see how genes may be clustered (genes become rows)
pca <- prcomp(test) #To see how Patients/Stage/Grades Cluster
#plot(pca$x[,1:2])
ipsgene <- c(c_gene, mcf7ips_overlap)
ipsgene <- unique(ipsgene)

#Adding back parameter of interest
pca$x <- rownames_to_column(as.data.frame(pca$x))
names(pca$x)[names(pca$x) == 'rowname'] <- 'PATIENT_ID'
pca$x <- merge(pca$x, extras, by="PATIENT_ID")

ggplot(as.data.frame(pca$x[,1:10])) + 
  geom_point(aes(x=PC1, y=PC2)) + #, color = pca$x$GRADE), size = 2) #+ 
  geom_text(aes(x=PC1, y=PC2, label = rownames(pca$x)), size = 4, hjust=-0.25, vjust=-0.5, color = ifelse(rownames(pca$x) %in% ipsgene, "red", "black"))

