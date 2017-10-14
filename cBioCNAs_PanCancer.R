library(tidyverse)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ISLR)
library(gplots)
library("Biobase")

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
cnadata <- read.csv(paste(file.path, "data_CNA.csv", sep=""), sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

cna <- subset(cnadata, select = -c(Entrez_Gene_Id, Cytoband))
cna <- melt(cna)
colnames(cna) <- c("GENE", "SAMPLE_ID", "VALUE")
cna$VALUE <- as.factor(cna$VALUE)

#### MERGING CLINICAL AND PATIENT DATA 
cancerCNADF <- merge(x = cna, y = patientdata, by = "SAMPLE_ID")

cnaDF <- cancerCNADF

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

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/HUViPSvsFiPS/")
huvfips_overlap <- data.frame(read.csv(paste(file.path, "HUViPSFiPSOverlap.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
huvfips_overlap <- as.factor(huvfips_overlap$x)

#overlap <- gene
#overlap <- c(mcf7ips_overlap)
#overlap <- c(mcf7ips_overlap, "TP53")
#overlap <- c(mcf7ips_overlap, "MYC")
#overlap <- c(mcf7ips_overlap, "FLOT2", "SPTBN1", "TRAP1", "DHX9")
overlap <- huvfips_overlap

stage_diffs <- subset(cnaDF, GENE %in% overlap) #Subset rows that are matches to those in the iPS targets

stage_diffs$cancertypes <- stage_diffs$TUMOR_TYPE

for (i in 1:length(stage_diffs$cancertypes)) {
  if (stage_diffs$cancertypes[i] == "glioma" | stage_diffs$cancertypes[i] =="medulloblastoma" | stage_diffs$cancertypes[i] =="neuroblastoma") {
    stage_diffs$cancertypes[i] <- "Brain"
  } else if (stage_diffs$cancertypes[i] == "AML" | stage_diffs$cancertypes[i] == "B-cell_ALL" | stage_diffs$cancertypes[i] == "CML" | stage_diffs$cancertypes[i] == "leukemia_other" | stage_diffs$cancertypes[i] == "multiple_myeloma" | stage_diffs$cancertypes[i] == "T-cell_ALL") {
    stage_diffs$cancertypes[i] <- "Leukemia"
  } else if (stage_diffs$cancertypes[i] == "breast") {
    stage_diffs$cancertypes[i] <- "Breast"
  } else if (stage_diffs$cancertypes[i] == "lung_NSC" | stage_diffs$cancertypes[i] == "lung_small_cell") {
    stage_diffs$cancertypes[i] <- "Lung"
  } else if (stage_diffs$cancertypes[i] == "lymphoma_Burkitt" | stage_diffs$cancertypes[i] == "lymphoma_DLBCL" | stage_diffs$cancertypes[i] == "lymphoma_Hodgkin" | stage_diffs$cancertypes[i] == "lymphoma_other") {
    stage_diffs$cancertypes[i] <- "Lymphoma"
  } else if (stage_diffs$cancertypes[i] == "chondrosarcoma" | stage_diffs$cancertypes[i] == "Ewings_Sarcoma" | stage_diffs$cancertypes[i] ==  "osteosarcoma") {
    stage_diffs$cancertypes[i] <- "Bone"
  } else if (stage_diffs$cancertypes[i] == "pancreas") {
    stage_diffs$cancertypes[i] <- "Pancreas"
  } else if (stage_diffs$cancertypes[i] == "colorectal") {
    stage_diffs$cancertypes[i] <- "Colorectal"
  } else {
    stage_diffs$cancertypes[i] <- "Other"
  }
  print(stage_diffs$cancertypes[i])
} 

stage_diffs$cancertypes <- as.factor(stage_diffs$cancertypes)

#values <- c("1", "2")
#values <- c("-1", "-2")
values <- c("-1", "-2", "1", "2")
stage_diffs <- filter(stage_diffs, VALUE %in% values)

#colors <- c("red", "darkred")
#colors <- c("dark green", "light green")
colors <- c("dark green", "light green", "red", "darkred")

ggplot(stage_diffs) + 
  geom_point(aes(x=PATIENT_ID, y=VALUE, color = cancertypes), position = position_jitter()) +
 # geom_boxplot(aes(x = TUMOR_STAGE, y = VALUE), fill = "firebrick", color = "darkblue") + 
  #scale_color_manual(values=colors) +
  theme_dark() +
  facet_grid(cancertypes ~ GENE)               

ggplot(stage_diffs) +
  geom_tile(aes(x = STAGE, y = PATIENT_ID, fill = VALUE)) +
  scale_fill_gradient(low = "black", high = "red", name = "") +
  facet_wrap(~GENE)



######
# cBio - Broad (luad_broad)
#####

### READ IN PATIENT DATA AND TIDY
file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/luad_broad/luad_broad/")
patientdata <- read.csv(paste(file.path, "data_clinical.csv", sep=""), sep = ",", stringsAsFactors = FALSE)

colnames(patientdata) <- lapply(patientdata[5,], as.character)

patientdata <- patientdata[-(1:5),]

patientdata <- subset(patientdata, select = c(PATIENT_ID, STAGE)) 

#### READ IN CLINICAL DATA AND TIDY
file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/luad_broad/luad_broad/")
cnadata <- read.csv(paste(file.path, "data_CNA.csv", sep=""), sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

cna <- melt(cnadata)
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

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/HUViPSvsFiPS/")
huvfips_overlap <- data.frame(read.csv(paste(file.path, "HUViPSFiPSOverlap.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
huvfips_overlap <- as.factor(huvfips_overlap$x)

overlap <- gene
#overlap <- c(mcf7ips_overlap)
#overlap <- c(mcf7ips_overlap, "TP53")
#overlap <- c(mcf7ips_overlap, "MYC")
#overlap <- c(mcf7ips_overlap, "FLOT2", "SPTBN1", "TRAP1", "DHX9")
#overlap <- huvfips_overlap

stage_diffs <- subset(cancerCNADF, GENE %in% overlap) #Subset rows that are matches to those in the iPS targets

stages <- factor(c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IIIC", "IV"))
stages <- ordered(stages, levels = c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IIIC", "IV"))
stage_diffs <- filter(stage_diffs, STAGE %in% stages)

#values <- c("1", "2")
#values <- c("-1", "-2")
values <- c("-1", "-2", "1", "2")
stage_diffs <- filter(stage_diffs, VALUE %in% values)

#colors <- c("red", "darkred")
#colors <- c("dark green", "light green")
colors <- c("dark green", "light green", "red", "darkred")

ggplot(stage_diffs) + 
  geom_point(aes(x=PATIENT_ID, y=STAGE, color = VALUE)) + #, position = position_jitter()) +
  # geom_boxplot(aes(x = TUMOR_STAGE, y = VALUE), fill = "firebrick", color = "darkblue") + 
  scale_color_manual(values=colors) +
  theme_dark() +
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

