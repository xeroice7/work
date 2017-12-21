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

#Split Stages
splits <- str_split_fixed(patientdata$AJCC_PATHOLOGIC_TUMOR_STAGE, " ", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '2'] <- 'Specific_Stage'
patientdata <- subset(patientdata, select=-c(1, AJCC_PATHOLOGIC_TUMOR_STAGE))

splits <- str_split_fixed(patientdata$Specific_Stage, "A|B|C", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'General_Stages'
patientdata <- subset(patientdata, select=-2)

#Split Tumor
splits <- str_split_fixed(patientdata$AJCC_TUMOR_PATHOLOGIC_PT, "a|b", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'General_Grade'
patientdata <- subset(patientdata, select=-2)
names(patientdata)[names(patientdata) == 'AJCC_TUMOR_PATHOLOGIC_PT'] <- 'Specific_Grade'

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

# Examine stage differences between our A549/iPS overlap without a somatic filter (37 genes)
#iPS to A549 Comparisons (No Somatic Filters - can include somatic hits)
a549huvips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/A549vsHUViPS/HUViPSA549OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
a549fips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/A549vsFiPS/FiPSA549OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
a549ips_overlap <- intersect(a549fips_overlap$V1, a549huvips_overlap$V1)

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
stage_diffs <- subset(stage_diffs, GENE %in% a549ips_overlap) #Subset rows that are matches to those in the targets
#stage_diffs <- subset(stage_diffs, GENE %in% surface_genes)
#stage_diffs <- subset(stage_diffs, AJCC_METASTASIS_PATHOLOGIC_PM %in% meta)
#stage_diffs <- subset(stage_diffs, AJCC_TUMOR_PATHOLOGIC_PT %in% grades)
#stage_diffs <- subset(stage_diffs, Specific_Stage %in% stages)
#stage_diffs <- filter(stage_diffs, ((Tumor == "T3") | (Tumor == "T2")))
stage_diffs <- filter(stage_diffs, General_Grade == "T1")
#stage_diffs <- filter(stage_diffs, General_Stages == "I")
#stage_diffs <- filter(stage_diffs, AJCC_METASTASIS_PATHOLOGIC_PM == "M0")

patient <- dcast(stage_diffs, PATIENT_ID+General_Grade+Specific_Grade+General_Stages+Specific_Stage+AJCC_METASTASIS_PATHOLOGIC_PM ~ GENE, value.var = "EXPRESSION_LEVEL")

#splits <- str_split_fixed(patient$PATIENT_ID, pattern = "-", "3")

#patient <- cbind.data.frame(splits, patient)

#patient <- subset(patient, select = -c(`1`, `2`, PATIENT_ID)) # Getting rid of Stage Column to prepare for a calculation matrix

#names(patient)[names(patient) == '3'] <- 'ID'
#patient <- arrange(patient, Specific_Stage, AJCC_TUMOR_PATHOLOGIC_PT, AJCC_METASTASIS_PATHOLOGIC_PM)

#patient <- patient %>% 
 # unite(PATIENT_ID, ID, Specific_Stage, AJCC_TUMOR_PATHOLOGIC_PT, AJCC_METASTASIS_PATHOLOGIC_PM, sep = "_", remove = TRUE)

#patient_match <- data.frame(patientDF$ID, patientDF$Stage)  # Creating a data frame to match patients with stages later

#colnames(patient_match) <- c("ID", "Stage")
#####
patient <- subset(patient, select = -c(General_Grade, Specific_Grade, General_Stages, Specific_Stage, AJCC_METASTASIS_PATHOLOGIC_PM))
patient[,1] <- as.character(patient[,1])
patient[,2:ncol(patient)] <- sapply(patient[,2:ncol(patient)],as.numeric)

sapply(patient, class)  #to check classes

patient <- patient %>% remove_rownames %>% column_to_rownames(var="PATIENT_ID") # Making the patient IDs the rownames, not the first column

#patientDF <- t(patientDF)

#redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)  # Making a color palette for the heatmap 


#Setting up the data for a heatmap/cluster analysis
patient <- patient[, apply(patient, 2, sum)!=0] ### Required to remove all the columns with 0 in them to get a distance matrix (use for "cor - standard deviation is zero' errors)
test <- as.matrix(patient)
test1 <- t(test)
grade1 <- test1
grade2 <- test1
grade3 <- test1
grade4 <- test1
#heatmap(test1, Colv=NA, col=greenred(10),scale="none")
#cor(t(test1))
#t <- 1-cor(t(test1))
#remove <- c("KRT9") #use if standard deviation is 0 on a row
#grade1 <- grade1[!rownames(grade1) %in% remove,] #use if standard deviation is 0 on a row
#grade4 <- grade4[!rownames(grade4) %in% remove,] #use if standard deviation is 0 on a row
#x <- cor(t(test1))
y <- cor(t(grade1), t(grade3))
y <- cor(t(grade1))
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
genecor <- y
output <- vector("double", ncol(genecor))
for (i in 1:ncol(genecor)){
  output[i] <- genecor[i,i]
}
extracted <- cbind.data.frame(rownames(y), output)
colnames(extracted) <- c("Gene", "Correlation")

#####To use heatmap.2 function:
extracted <- extracted %>% remove_rownames %>% column_to_rownames(var="Gene")
extracted <- cbind(extracted, extracted)
extracted <- as.matrix(extracted)
ext_clust <- hclust(dist(1-extracted))
heatmap.2(extracted,
          Rowv=as.dendrogram(ext_clust), 
          Colv=NA, 
          col=bluered(256), 
          #breaks = breaks,
          cexRow = 0.7, 
          cexCol = 0.7, 
          scale = "none", 
          trace = "none")

######To use geom_tile function:
extracted <- melt(extracted)
ggplot(extracted, aes(x=variable, y=Gene)) +
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") + 
  scale_y_discrete(name="", limits = rev(levels(extracted$Gene))) + 
  geom_text(aes(x=variable, y=Gene, label=round(value, digits = 2)), size=2) + 
  xlab("Correlation between Grade 1 vs Grade 4") + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=6),
        legend.text = element_text(size=8), 
        axis.title.x = element_text(size=8)) 

## For plot to sort numerically
extracted <- arrange(extracted, desc(value))
ggplot(extracted) +
  geom_tile(aes(x=variable, y=Gene, fill=value)) + 
  scale_fill_distiller(palette = "RdBu") + 
  scale_y_discrete(name="", limits = extracted$Gene) + 
  geom_text(aes(x=variable, y=Gene, label = round(value, digits = 2)), size=2) + 
  xlab("Correlation between Grade 1 vs Grade 4") + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=6),
        legend.text = element_text(size=8), 
        axis.title.x = element_text(size=8)) 

############
#To look at mean stage differences
############
t_patient <- t(patient)
mean <- rowMeans(t_patient, na.rm = FALSE)
t_patient <- cbind.data.frame(t_patient, mean)

t_patient_lo <- t_patient #For Low
t_patient_2 <- t_patient
t_patient_3 <- t_patient
t_patient_hi <- t_patient #For Hi 

####### T-tests for difference in Stages
stage1 <- data.frame(t(t_patient_lo))
stage4 <- data.frame(t(t_patient_hi))
results <- mapply(t.test, stage1, stage4)
results <- plyr::ldply(results["p.value",], data.frame)
colnames(results) <- c("Genes", "pvalues")
results <- arrange(results, desc(pvalues))
results$Stage <- "Grade3"

results$colorscale <- cut(results$pvalues, breaks = c(0,0.05,0.1,0.25,0.5,1),right = FALSE)

ggplot(results) + 
  geom_tile(aes(x=Stage, y=Genes, fill = colorscale), color = "white") + 
  scale_fill_brewer(palette = "PRGn") +
  scale_y_discrete(name="", limits = results$Genes) + 
  geom_text(aes(x=Stage, y=Genes, label = round(pvalues, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8)) 

write.csv(results, "Lung_LUAD_TCGA_Stats_Meta.csv")
#######

colnames(t_patient_lo)[colnames(t_patient_lo) == 'mean'] <- 'mean_1'
colnames(t_patient_2)[colnames(t_patient_2) == 'mean'] <- 'mean_2'
colnames(t_patient_3)[colnames(t_patient_3) == 'mean'] <- 'mean_3'
colnames(t_patient_hi)[colnames(t_patient_hi) == 'mean'] <- 'mean_4'

t_patient_m1 <- merge(t_patient_lo, t_patient_2, by = "row.names")
t_patient_m2 <- merge(t_patient_3, t_patient_hi, by = "row.names")
t_patient_means <- merge(t_patient_m1, t_patient_m2, by = "Row.names")
t_patient_means <- t_patient_m1

t_patient_means$mean_diff_1 <- t_patient_means$mean_1-t_patient_means$mean_1
t_patient_means$mean_diff_2 <- t_patient_means$mean_2-t_patient_means$mean_1
t_patient_means$mean_diff_3 <- t_patient_means$mean_3-t_patient_means$mean_1
t_patient_means$mean_diff_4 <- t_patient_means$mean_4-t_patient_means$mean_1

cast_patient <- cbind.data.frame(t_patient_means$Row.names, t_patient_means$mean_diff_1, t_patient_means$mean_diff_4)#, t_patient_means$mean_diff_3, t_patient_means$mean_diff_4)
colnames(cast_patient) <- c("Genes", "Stage1", "Stage2", "Stage3", "Stage4")
#colnames(cast_patient) <- c("Genes", "Grade1", "Grade2", "Grade3", "Grade4")
colnames(cast_patient) <- c("Genes", "Grade1", "Grade3")
melt_patient <- melt(cast_patient)
colnames(melt_patient) <- c("Genes", "Stages", "Difference")
melt_patient$Genes <- as.factor(melt_patient$Genes)
#plot_patient <- cbind.data.frame(t_patient$Row.names, t_patient$mean_diff)

ggplot(melt_patient, aes(x=Stages, y=Genes)) + 
  geom_tile(aes(fill=Difference)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="", limits = rev(levels(melt_patient$Genes))) + 
  geom_text(aes(x=Stages, y=Genes, label = round(Difference, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

# To see descending values
#melt_patient <- filter(melt_patient, Stages == "ER+/HER2-")
melt_patient <- filter(melt_patient, Stages == "Grade3")
melt_patient <- arrange(melt_patient, desc(Difference))
ggplot(melt_patient) + 
  geom_tile(aes(x=Stages, y=Genes, fill = Difference)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="", limits = melt_patient$Genes) + 
  geom_text(aes(x=Stages, y=Genes, label = round(Difference, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

write.csv(cast_patient, "Lung_LUAD_TCGA_MeanExpressionChange_MCF7iPSOverlap_ByStage.csv")


##### To find most likely/intersting/promising candidates - BY GRADE:
#MCF7/iPS Overlap - GRADE
hits <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/Analysis/MeanExpressionChanges/MCF7iPSOverlap/Lung_LUAD_TCGA_MeanExpressionChange_MCF7iPSOverlap_ByGrade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

#Total iPS Overlap - GRADE
hits <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/Analysis/MeanExpressionChanges/TotaliPSOverlap/Lung_LUAD_TCGA_MeanExpressionChange_TotaliPSOverlap_ByGrade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

#A549/iPS Overlap - GRADE
hits <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/Analysis/MeanExpressionChanges/A549iPSOverlap/Lung_LUAD_TCGA_MeanExpressionChange_A549iPSOverlap_ByGrade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

hits_up <- filter(hits, Grade2 > 0 & Grade3 > Grade2 & Grade4 > Grade3)
hits_down <- filter(hits, Grade2 < 0 & Grade3 < Grade2, Grade4 < Grade3)

hits_up$DeltaG3toG2 <- hits_up$Grade3 - hits_up$Grade2
hits_up$DeltaG4toG3 <- hits_up$Grade4 - hits_up$Grade3
hits_up <- arrange(hits_up, -DeltaG4toG3, -DeltaG3toG2, -Grade2)

hits_down$DeltaG3toG2 <- hits_down$Grade3 - hits_down$Grade2
hits_down$DeltaG4toG3 <- hits_down$Grade4 - hits_down$Grade3
hits_down <- arrange(hits_down, DeltaG4toG3, DeltaG3toG2, Grade2)

hits_up <- subset(hits_up, select = -X)
hits_up$Genes <- factor(hits_up$Genes, levels = hits_up$Genes)
m_hits_up <- melt(hits_up)

ggplot(m_hits_up, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(hits_up$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

hits_down <- subset(hits_down, select = -X)
hits_down$Genes <- factor(hits_down$Genes, levels = hits_down$Genes)
m_hits_down <- melt(hits_down)

ggplot(m_hits_down, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(hits_down$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

##### To find most likely/intersting/promising candidates - BY STAGE:
#MCF7/iPS Overlap - STAGE
st_hits <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/Analysis/MeanExpressionChanges/MCF7iPSOverlap/Lung_LUAD_TCGA_MeanExpressionChange_MCF7iPSOverlap_ByStage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

#Total iPS Overlap - STAGE
st_hits <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/Analysis/MeanExpressionChanges/TotaliPSOverlap/Lung_LUAD_TCGA_MeanExpressionChange_TotaliPSOverlap_ByStage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

#A549/iPS Overlap - STAGE
st_hits <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/Analysis/MeanExpressionChanges/A549iPSOverlap/Lung_LUAD_TCGA_MeanExpressionChange_A549iPSOverlap_ByStage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

st_hits_up <- filter(st_hits, Stage2 > 0 & Stage3 > Stage2 & Stage4 > Stage3)
st_hits_down <- filter(st_hits, Stage2 < 0 & Stage3 < Stage2 & Stage4 < Stage3)

st_hits_up$DeltaS3toS2 <- st_hits_up$Stage3 - st_hits_up$Stage2
st_hits_up$DeltaS4toS3 <- st_hits_up$Stage4 - st_hits_up$Stage3
st_hits_up <- arrange(st_hits_up, -DeltaS4toS3, -DeltaS3toS2, -Stage2)

st_hits_down$DeltaS3toS2 <- st_hits_down$Stage3 - st_hits_down$Stage2
st_hits_down$DeltaS4toS3 <- st_hits_down$Stage4 - st_hits_down$Stage3
st_hits_down <- arrange(st_hits_down, DeltaS4toS3, DeltaS3toS2, Stage2)

st_hits_up <- subset(st_hits_up, select = -X)
st_hits_up$Genes <- factor(st_hits_up$Genes, levels = st_hits_up$Genes)
m_st_hits_up <- melt(st_hits_up)

ggplot(m_st_hits_up, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(st_hits_up$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

st_hits_down <- subset(st_hits_down, select = -X)
st_hits_down$Genes <- factor(st_hits_down$Genes, levels = st_hits_down$Genes)
m_st_hits_down <- melt(st_hits_down)

ggplot(m_st_hits_down, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(st_hits_down$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

intersect(hits_up$Genes, st_hits_up$Genes)
intersect(hits_down$Genes, st_hits_down$Genes)
