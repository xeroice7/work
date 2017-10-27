library(tidyverse)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ISLR)
library(gplots)
library("Biobase")

######
# cBio - METABRIC Nature 2012 & Nat Com 2016 Dataset
#####

### READ IN PATIENT DATA AND TIDY
tumordata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/brca_metabric/data_clinical_supp_sample.csv", sep = ",", stringsAsFactors = FALSE)
patientdata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/brca_metabric/data_clinical_supp_patient.csv", sep = ",", stringsAsFactors = FALSE)

patientDF <- merge(x = tumordata, y = patientdata, by = "PATIENT_ID", all = TRUE)

patientDF <- subset(patientDF, select = c(PATIENT_ID, GRADE, TUMOR_STAGE)) 


#### READ IN CLINICAL DATA AND TIDY
cnadata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/brca_metabric/data_CNA.csv", sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

cna <- subset(cnadata, select = -Entrez_Gene_Id)
cna <- melt(cna)
colnames(cna) <- c("GENE", "PATIENT_ID", "VALUE")
cna$VALUE <- as.factor(cna$VALUE)


#### MERGING CLINICAL AND PATIENT DATA 
cancerCNADF <- merge(x = cna, y = patientDF, by = "PATIENT_ID")

cnaDF <- cancerCNADF

# Examine stage differences between our MCF7/iPS overlap without a somatic filter (29 genes)
#iPS to MCF7 Comparisons (No Somatic Filters - can include somatic hits)
mcf7huvips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsHUViPS/HUViPSMCF7OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
mcf7fips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsFiPS/FiPSMCF7OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
mcf7ips_overlap <- intersect(mcf7fips_overlap$V1, mcf7huvips_overlap$V1)

#Unique iPS targets(no somatic source - 34 genes)
gene <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/CSVs/TotalUniqueiPSGENEandPID.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)


huvfips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/HUViPSvsFiPS/HUViPSFiPSOverlap.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
huvfips_overlap <- as.factor(huvfips_overlap$x)

#All Human Cell surface-specific targets indicated by GO analysis on GO site 
data <- read.csv("~/Desktop/data.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) # Import CSV of surface expression
surface_genes <- data$V1 
surface_genes <- unique(surface_genes) 

#A smaller subset of genes pulled from FiPS and HUViPS MS sets to test clustering
testgenes <- read.csv("~/Desktop/TestGeneSetUnique.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) # Import CSV of surface expression
testgenes <- as.factor(testgenes$x)
c_testgenes <- as.character(testgenes)
c_gene <- as.character(gene)
totaltestgenes <- c(c_testgenes, c_gene)
totaltestgenes <- unique(totaltestgenes)
totaltestgenes <- as.factor(totaltestgenes)

#overlap <- gene
#overlap <- c(mcf7ips_overlap)
#overlap <- c(mcf7ips_overlap, "TP53")
#overlap <- c(mcf7ips_overlap, "MYC")
#overlap <- c(mcf7ips_overlap, "FLOT2", "SPTBN1", "TRAP1", "DHX9")
overlap <- huvfips_overlap

stage_diffs <- subset(cnaDF, GENE %in% overlap) #Subset rows that are matches to those in the iPS targets

grades <- c("2", "3")
stage_diffs <- filter(stage_diffs, GRADE %in% grades)

#values <- c("1", "2")
#values <- c("-1", "-2")
values <- c("-1", "-2", "1", "2")
stage_diffs <- filter(stage_diffs, VALUE %in% values)

#colors <- c("red", "darkred")
#colors <- c("dark green", "light green")
colors <- c("dark green", "light green", "red", "darkred")

ggplot(stage_diffs) + 
  geom_point(aes(x=PATIENT_ID, y=GRADE, color = VALUE), position = position_jitter()) +
 # geom_boxplot(aes(x = TUMOR_STAGE, y = VALUE), fill = "firebrick", color = "darkblue") + 
  scale_color_manual(values=colors) +
  theme_dark() +
  facet_wrap(~GENE)               

ggplot(stage_diffs) +
  geom_tile(aes(x = STAGE, y = PATIENT_ID, fill = VALUE)) +
  scale_fill_gradient(low = "black", high = "red", name = "") +
  facet_wrap(~GENE)

##### PCA Analysis with Metabric Dataset
#### FILTERING DATA WITH OUR CRITERIA OF INTEREST
grades <- factor(c("1", "2", "3"))
grades <- ordered(grades, levels = c("1", "2", "3"))
stages <- factor(c("1", "2", "3", "4"))
stages <- ordered(stages, levels = c("1", "2", "3", "4"))
stage_diffs <- subset(cnaDF, GENE %in% huvfips_overlap)
#stage_diffs <- subset(patientDF, GENE %in% huvfips_overlap) #Subset rows that are matches to those in the targets
#stage_diffs <- subset(stage_diffs, GENE %in% surface_genes)
stage_diffs <- subset(stage_diffs, GRADE %in% grades)
stage_diffs <- subset(stage_diffs, TUMOR_STAGE %in% stages)

#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 4)
#stage_diffs <- filter(stage_diffs, ((GRADE == 1) | (GRADE == 2)))
#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 2)
#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 3)
#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 4)
#stage_diffs <- filter(stage_diffs, (GRADE == 3))

patient <- dcast(stage_diffs, PATIENT_ID+GRADE+TUMOR_STAGE ~ GENE, value.var = "VALUE")

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
patient <- patient %>% drop_na()
patient <- subset(patient, select = -c(GRADE, TUMOR_STAGE))
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

#Adding back parameter of interest for patient clustering
pca$x <- rownames_to_column(as.data.frame(pca$x))
names(pca$x)[names(pca$x) == 'rowname'] <- 'PATIENT_ID'
pca$x <- merge(pca$x, extras, by="PATIENT_ID")


ggplot(as.data.frame(pca$x[,1:10])) + 
  geom_point(aes(x=PC1, y=PC2), color = pca$x$STAGES, size = 2) + 
  geom_text(aes(x=PC1, y=PC2, label = rownames(pca$x)), size = 3, hjust=-0.25, vjust=-0.5, color = ifelse(rownames(pca$x) %in% ipsgene, "red", "black"))

t <- as.data.frame(pca$x)
t <- column_to_rownames(t, var="PATIENT_ID")

t <- filter(t, ((TUMOR == 1) | (TUMOR == 3)))

#GGPairs Test
ggpairs(t, columns = 1:6, mapping = aes(color = TUMOR))

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


######
# cBio - TCGA (Cell 2015) Dataset
#####

### READ IN PATIENT DATA AND TIDY
patientdata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/brca_tcga_pub/brca_tcga_pub/data_clinical.csv", sep = ",", stringsAsFactors = FALSE)

colnames(patientdata) <- lapply(patientdata[5,], as.character)

patientdata <- patientdata[-(1:5),]

patientdata <- subset(patientdata, select = c(PATIENT_ID, Tumor, `AJCC Stage`, Metastasis, OS_STATUS)) 

splits <- str_split_fixed(patientdata$`AJCC Stage`, " ", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '2'] <- 'Stages'
patientdata <- subset(patientdata, select=-1)

splits <- str_split_fixed(patientdata$Stage, "A|B|C", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'Stage'
patientdata <- subset(patientdata, select=-c(2, `AJCC Stage`))


#### READ IN CLINICAL DATA AND TIDY
cnadata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/brca_tcga_pub/brca_tcga_pub/data_CNA.csv", sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

cna <- subset(cnadata, select = -Entrez_Gene_Id)
cna <- melt(cna)
colnames(cna) <- c("GENE", "PATIENT_ID", "VALUE")
cna$VALUE <- as.factor(cna$VALUE)


#### MERGING CLINICAL AND PATIENT DATA 
cancerCNADF <- merge(x = cna, y = patientdata, by = "PATIENT_ID")

# Examine stage differences between our MCF7/iPS overlap without a somatic filter (29 genes)
#iPS to MCF7 Comparisons (No Somatic Filters - can include somatic hits)
mcf7huvips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsHUViPS/HUViPSMCF7OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
mcf7fips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsFiPS/FiPSMCF7OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))

mcf7ips_overlap <- intersect(mcf7fips_overlap$V1, mcf7huvips_overlap$V1)

#overlap <- c(mcf7ips_overlap)
#overlap <- c(mcf7ips_overlap, "TP53")
#overlap <- c(mcf7ips_overlap, "MYC")
#overlap <- c(mcf7ips_overlap, "FLOT2", "SPTBN1", "TRAP1", "DHX9")

stage_diffs <- subset(cancerCNADF, GENE %in% huvfips_overlap) #Subset rows that are matches to those in the iPS targets

grades <- c("T1", "T2", "T3")
stage_diffs <- filter(stage_diffs, Tumor %in% grades)

values <- c("1", "2")
#values <- c("-1", "-2")
values <- c("-1", "-2", "0", "1", "2")
stage_diffs <- filter(stage_diffs, VALUE %in% values)

#colors <- c("red", "darkred")
#colors <- c("dark green", "light green")
colors <- c("dark green", "light green", "white", "red", "darkred")

ggplot(stage_diffs) + 
  geom_point(aes(x=PATIENT_ID, y=Tumor, color = VALUE), position = position_jitter()) +
  # geom_boxplot(aes(x = TUMOR_STAGE, y = VALUE), fill = "firebrick", color = "darkblue") + 
  facet_wrap(~GENE)               

ggplot(stage_diffs) +
  geom_tile(aes(x = STAGE, y = PATIENT_ID, fill = VALUE)) +
  scale_fill_gradient(low = "black", high = "red", name = "") +
  facet_wrap(~GENE)


# Metastasis
#overlap <- c(mcf7ips_overlap)
#overlap <- c(mcf7ips_overlap, "TP53")
overlap <- c(mcf7ips_overlap, "MYC")

stage_diffs <- subset(cancerCNADF, GENE %in% overlap) #Subset rows that are matches to those in the iPS targets

meta <- c("M0", "M1") 
#meta <- c("M0")
#meta <- c("M1")

stage_diffs <- subset(stage_diffs, Metastasis %in% meta)

stage_diffs$Metastasis <- as.factor(stage_diffs$Metastasis)
stage_diffs$OS_STATUS <- as.factor(stage_diffs$OS_STATUS)
shapes <- c(16,4)
#colors <- c("dark green", "light green", "red", "darkred")
colors <- c("red", "darkred")

ggplot(stage_diffs) + 
  geom_point(aes(x=PATIENT_ID, y=Tumor, color = VALUE, shape = Metastasis), position = position_jitter()) +
  # geom_boxplot(aes(x = TUMOR_STAGE, y = VALUE), fill = "firebrick", color = "darkblue") + 
  scale_shape_manual(values=shapes) + scale_color_manual(values=colors) + 
  theme_dark() + 
  facet_wrap(~GENE)               

ggplot(stage_diffs) +
  geom_tile(aes(x = STAGE, y = PATIENT_ID, fill = VALUE)) +
  scale_fill_gradient(low = "black", high = "red", name = "") +
  facet_wrap(~GENE)


patient <- dcast(stage_diffs, PATIENT_ID+Stage+Stages+Tumor+Metastasis+OS_STATUS ~ GENE, value.var = "VALUE")
#patient <- filter(patient, Stage == "IV")

extras <- cbind.data.frame(patient$PATIENT_ID, patient$Stage, patient$Stages, patient$Tumor, patient$Metastasis, patient$OS_STATUS)
colnames(extras) <- c("PATIENT_ID", "STAGE", "STAGES", "TUMOR", "METASTASIS", "OS_STATUS")

patient <- patient %>% drop_na()
patient <- subset(patient, select = -c(Stage, Stages, Tumor, Metastasis, OS_STATUS))
patient[,1] <- as.character(patient[,1])

patient[,2:ncol(patient)] <- sapply(patient[,2:ncol(patient)],as.numeric)

sapply(patient, class)  #to check classes

patient <- patient %>% remove_rownames %>% column_to_rownames(var="PATIENT_ID") # Making the patient IDs the rownames, not the first column

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

### To plot patient characteristics:
ggplot(as.data.frame(pca$x[,1:10])) + 
  geom_point(aes(x=PC1, y=PC2, color = pca$x$METASTASIS), size = 2) 

### To plot genes:
ggplot(as.data.frame(pca$x[,1:10])) + 
  geom_point(aes(x=PC1, y=PC2)) + 
  geom_text(aes(x=PC1, y=PC2, label = rownames(pca$x)), size = 3, hjust=-0.25, vjust=-0.5) #color = ifelse(rownames(pca$x) %in% ipsgene, "red", "black"))

t <- as.data.frame(pca$x)
t <- column_to_rownames(t, var="PATIENT_ID")

t <- filter(t, ((STAGE == "I") | (STAGE == "IV")))

#GGPairs Test
ggpairs(t, columns = 1:6, mapping = aes(color = STAGE))


######
# cBio - Lung TCGA Broad 2016
#####

### READ IN PATIENT DATA AND TIDY
file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/nsclc_tcga_broad_2016/nsclc_tcga_broad_2016/")
patientdata <- read.csv(paste(file.path, "data_clinical.csv", sep=""), sep = ",", stringsAsFactors = FALSE)

colnames(patientdata) <- lapply(patientdata[5,], as.character)

patientdata <- patientdata[-(1:5),]

patientdata <- subset(patientdata, select = c(SAMPLE_ID, STAGE, OS_STATUS, VITAL_STATUS)) 
names(patientdata)[names(patientdata) == 'SAMPLE_ID'] <- 'PATIENT_ID' 

splits <- str_split_fixed(patientdata$STAGE, "A|B|C", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'Stage'
patientdata <- subset(patientdata, select=-c(2, STAGE))


#### READ IN CLINICAL DATA AND TIDY
file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/nsclc_tcga_broad_2016/nsclc_tcga_broad_2016/")
cnadata <- read.csv(paste(file.path, "data_CNA.csv", sep=""), sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

cna <- subset(cnadata, select = -Entrez_Gene_Id)
cna <- melt(cna)
colnames(cna) <- c("GENE", "PATIENT_ID", "VALUE")
cna$VALUE <- as.factor(cna$VALUE)


#### MERGING CLINICAL AND PATIENT DATA 
cancerCNADF <- merge(x = cna, y = patientdata, by = "PATIENT_ID")

# Examine stage differences between our MCF7/iPS overlap without a somatic filter (29 genes)
#iPS to MCF7 Comparisons (No Somatic Filters - can include somatic hits)
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsHUViPS/")
mcf7huvips_overlap <- data.frame(read.csv(paste(file.path, "HUViPSMCF7OverlapGENE.csv", sep=""), header = FALSE, sep = ",", stringsAsFactors = FALSE))

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsFiPS/")
mcf7fips_overlap <- data.frame(read.csv(paste(file.path, "FiPSMCF7OverlapGENE.csv", sep=""), header = FALSE, sep = ",", stringsAsFactors = FALSE))

mcf7ips_overlap <- intersect(mcf7fips_overlap$V1, mcf7huvips_overlap$V1)

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/CSVs/")
gene <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

overlap <- gene
#overlap <- c(mcf7ips_overlap)
#overlap <- c(mcf7ips_overlap, "TP53")
#overlap <- c(mcf7ips_overlap, "MYC")
#overlap <- c(mcf7ips_overlap, "TP53", "MYC")

stage_diffs <- subset(cancerCNADF, GENE %in% overlap) #Subset rows that are matches to those in the iPS targets

grades <- c("1", "2", "3")
stage_diffs <- filter(stage_diffs, GRADE %in% grades)

values <- c("-2", "-1", "1", "2")
#values <- c("-1", "-2")
stage_diffs <- filter(stage_diffs, VALUE %in% values)

colors <- c("dark green", "light green", "red", "darkred")

ggplot(stage_diffs) + 
  geom_point(aes(x=PATIENT_ID, y=GRADE, color = VALUE), position = position_jitter()) +
  # geom_boxplot(aes(x = TUMOR_STAGE, y = VALUE), fill = "firebrick", color = "darkblue") + 
  scale_color_manual(values=colors) +
  theme_dark() +
  facet_wrap(~GENE)               

ggplot(stage_diffs) +
  geom_tile(aes(x = STAGE, y = PATIENT_ID, fill = VALUE)) +
  scale_fill_gradient(low = "black", high = "red", name = "") +
  facet_wrap(~GENE)

