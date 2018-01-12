####
# cBio - Breast Adenocarcinoma Metabric - 2509 Patients
####

#Necessary Packages ===========================================
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ISLR)
library(gplots)

#Load Gene Filters =========================================
#Total iPS overlaps (including somatic targets - 93 genes)
huvfips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/HUViPSvsFiPS/HUViPSFiPSOverlap.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
huvfips_overlap <- as.factor(huvfips_overlap$x)

#All Human Cell surface-specific targets indicated by GO analysis on GO site 
data <- read.csv("~/Desktop/HumanGOCellSurface.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) # Import CSV of surface expression
surface_genes <- data$V1 
surface_genes <- unique(surface_genes) 

#Unique iPS targets(no somatic source - 34 genes)
gene <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/CSVs/TotalUniqueiPSGENEandPID.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

#iPS to MCF7 Comparisons (No Somatic Filters - can include somatic hits) - 29 genes
mcf7huvips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsHUViPS/HUViPSMCF7OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
mcf7fips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsFiPS/FiPSMCF7OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
mcf7ips_overlap <- intersect(mcf7fips_overlap$V1, mcf7huvips_overlap$V1)

#A smaller subset of genes pulled from FiPS and HUViPS MS sets to test clustering - 390 genes
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

#Load Patient Data ============================================
tumordata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/data_clinical_supp_sample.csv", sep = ",", stringsAsFactors = FALSE)
patientdata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/data_clinical_supp_patient.csv", sep = ",", stringsAsFactors = FALSE)

patientDF <- merge(x = tumordata, y = patientdata, by = "PATIENT_ID", all = TRUE)

#Tidy Patient Data ============================================
splits <- str_split_fixed(patientDF$THREEGENE, " ", 2)
patientDF <- cbind.data.frame(splits, patientDF)
names(patientDF)[names(patientDF) == '2'] <- 'Proliferation'
names(patientDF)[names(patientDF) == '1'] <- 'Markers'

splits <- str_split_fixed(patientDF$Proliferation, " ", 2)
patientDF <- cbind.data.frame(splits, patientDF)
names(patientDF)[names(patientDF) == '1'] <- 'RATE_OF_PROLIF'
patientDF <- subset(patientDF, select= -c(2, Proliferation))

#Load Clinical Data ================================================
expressiondata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/data_expression.csv", sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

#Tidy Clinical Data ================================================
expressiondata[expressiondata == "null"] <- NA

expressiondata <- subset(expressiondata, select = -Entrez_Gene_Id)

expressiondata <- expressiondata %>% drop_na()

expressionDF <- melt(expressiondata, id.vars = "Hugo_Symbol")

colnames(expressionDF) <- c("GENE", "PATIENT_ID", "EXPRESSION_LEVEL")

#Merge Patient and Clinical Data ==================================== 
patientDF <- merge(x = expressionDF, y = patientDF, by = "PATIENT_ID")

#Additional Filter Parameters =======================================
grades <- factor(c("1", "2", "3"))
grades <- ordered(grades, levels = c("1", "2", "3"))
stages <- factor(c("1", "2", "3", "4"))
stages <- ordered(stages, levels = c("1", "2", "3", "4"))

stage_diffs <- patientDF
#stage_diffs <- subset(patientDF, GENE %in% totaltestgenes)
stage_diffs <- subset(stage_diffs, GENE %in% huvfips_overlap) #Subset rows that are matches to those in the targets
#stage_diffs <- subset(stage_diffs, GENE %in% gene)
#stage_diffs <- subset(stage_diffs, GRADE %in% grades)
#stage_diffs <- subset(stage_diffs, TUMOR_STAGE %in% stages)

#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 1)
#stage_diffs <- filter(stage_diffs, ((GRADE == 1) | (GRADE == 2)))
stage_diffs <- filter(stage_diffs, GRADE == 3)
#stage_diffs <- filter(stage_diffs, (ER_STATUS == "ER-" & HER_STATUS == "HER2-"))
#stage_diffs <- filter(stage_diffs, HER_STATUS == "HER2+")
#stage_diffs <- filter(stage_diffs, RATE_OF_PROLIF == "Low")

patient <- dcast(stage_diffs, PATIENT_ID+GRADE+TUMOR_STAGE+ER_STATUS+HER_STATUS+RATE_OF_PROLIF ~ GENE, value.var = "EXPRESSION_LEVEL")

#Combine Patients with Clinical Data in Names =====================================
#extras <- cbind.data.frame(patient$PATIENT_ID, patient$GRADE, patient$TUMOR_STAGE)
#colnames(extras) <- c("PATIENT_ID", "TUMOR", "STAGES")
#
#splits <- str_split_fixed(patient$PATIENT_ID, pattern = "-", "2")

#patient <- cbind.data.frame(splits, patient)

#patient <- subset(patient, select = -c(`1`, PATIENT_ID)) # Getting rid of Stage Column to prepare for a calculation matrix

#names(patient)[names(patient) == '2'] <- 'ID'

#patient <- arrange(patient, TUMOR_STAGE, GRADE)

#patient <- patient %>% 
#  unite(PATIENT_ID, ID, TUMOR_STAGE, GRADE, sep = "_", remove = TRUE)

#Setting up patient matrix ========================================================
patient <- subset(patient, select = -c(GRADE, TUMOR_STAGE, RATE_OF_PROLIF, ER_STATUS, HER_STATUS))
patient[,1] <- as.character(patient[,1])

patient[,2:ncol(patient)] <- sapply(patient[,2:ncol(patient)],as.numeric)

sapply(patient, class)  #to check classes

patient <- patient %>% remove_rownames %>% column_to_rownames(var="PATIENT_ID") # Making the patient IDs the rownames, not the first column

#patientDF <- t(patientDF)

#redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)  # Making a color palette for the heatmap 

#Making patient matrix ===========================================================
patient <- patient[, apply(patient, 2, sum)!=0] ### Required to remove all the columns with 0 in them to get a distance matrix
test <- as.matrix(patient)
test1 <- t(test)
grade1 <- test1
grade2 <- test1
grade3 <- test1
#grade4 <- test1
#heatmap(test1, Colv=NA, col=greenred(10),scale="none")
#cor(t(test1))
t <- 1-cor(t(test1))
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
          trace = "none", 
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
  xlab("Correlation between Stage 1 vs Stage 4") + 
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
  xlab("Correlation between Stage 1 vs Stage 4") + 
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

write.csv(results, "Breast_Metabric_Stats_Grade.csv")
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

cast_patient <- cbind.data.frame(t_patient_means$Row.names, t_patient_means$mean_diff_1, t_patient_means$mean_diff_2)#, t_patient_means$mean_diff_3)#, t_patient_means$mean_diff_4)
#colnames(cast_patient) <- c("Genes", "Stage1", "Stage4")#, "Stage2", "Stage3")#, "Stage4")
#colnames(cast_patient) <- c("Genes", "HER2+", "ER+/HER2-", "ER-/HER2-")
colnames(cast_patient) <- c("Genes", "Grade1", "Grade3")#, "Grade3")
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


write.csv(cast_patient, "TripNeg_TotaliPSOverlap_ByGrade.csv")

#To cluster similar trends:
cast_patient <- cast_patient %>% remove_rownames %>% column_to_rownames(var="Genes")
cast_patient <- as.matrix(cast_patient)

heatmap.2(cast_patient, 
          Rowv=T, 
          Colv=NA, 
          col=bluered(256), 
          #breaks = breaks,
          cexRow = 0.5, 
          cexCol = 0.7, 
          scale = "none", 
          trace = "none")


plot_patient <- plot_patient %>% remove_rownames %>% column_to_rownames(var="t_patient$Row.names")
plot_patient <- as.matrix(plot_patient)
plot_patient <- cbind(plot_patient, plot_patient )

heatmap.2(plot_patient, 
          Rowv=NA, 
          Colv=NA, 
          col=bluered(256), 
          #breaks = breaks,
          cexRow = 0.7, 
          cexCol = 0.7, 
          scale = "none", 
          trace = "none")

##### To find most likely/intersting/promising candidates - BY GRADE:
#MCF7/iPS Overlap - GRADE
tripneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/MCF7iPSOverlap/TripNeg_MCF7iPSOverlap_ByGrade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
doubleneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/MCF7iPSOverlap/DoubleNeg_MCF7iPSOverlap_ByGrade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
singleneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/MCF7iPSOverlap/Her2Pos_MCF7iPSOverlap_ByGrade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
#Total iPS Overlap - GRADE
tripneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/TotaliPSOverlap/TripNeg_TotaliPSOverlap_ByGrade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
doubleneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/TotaliPSOverlap/DoubleNeg_TotaliPSOverlap_ByGrade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
singleneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/TotaliPSOverlap/Her2Pos_TotaliPSOverlap_ByGrade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

tripneg_up <- filter(tripneg, Grade2 > 0 & Grade3 > Grade2)
tripneg_down <- filter(tripneg, Grade2 < 0 & Grade3 < Grade2)

tripneg_up$DeltaG3toG2 <- tripneg_up$Grade3 - tripneg_up$Grade2
tripneg_up <- arrange(tripneg_up, -DeltaG3toG2, -Grade2)
#tripneg_up <- tripneg_up[1:7,]

tripneg_down$DeltaG3toG2 <- tripneg_down$Grade3 - tripneg_down$Grade2
tripneg_down <- arrange(tripneg_down, DeltaG3toG2, Grade2)
#tripneg_down <- tripneg_down[1:7,]

doubleneg_up <- filter(doubleneg, Grade2 > 0 & Grade3 > Grade2)
doubleneg_down <- filter(doubleneg, Grade2 < 0 & Grade3 < Grade2)

doubleneg_up$DeltaG3toG2 <- doubleneg_up$Grade3 - doubleneg_up$Grade2
doubleneg_up <- arrange(doubleneg_up, -DeltaG3toG2, -Grade2)
#doubleneg_up <- doubleneg_up[1:7,]

doubleneg_down$DeltaG3toG2 <- doubleneg_down$Grade3 - doubleneg_down$Grade2
doubleneg_down <- arrange(doubleneg_down, DeltaG3toG2, Grade2)
#doubleneg_down <- doubleneg_down[1:7,]

singleneg_up <- filter(singleneg, Grade2 > 0 & Grade3 > Grade2)
singleneg_down <- filter(singleneg, Grade2 < 0 & Grade3 < Grade2)

singleneg_up$DeltaG3toG2 <- singleneg_up$Grade3 - singleneg_up$Grade2
singleneg_up <- arrange(singleneg_up, -DeltaG3toG2, -Grade2)
#singleneg_up <- singleneg_up[1:7,]

singleneg_down$DeltaG3toG2 <- singleneg_down$Grade3 - singleneg_down$Grade2
singleneg_down <- arrange(singleneg_down, DeltaG3toG2, Grade2)
#singleneg_down <- singleneg_down[1:7,]

veryinteresting_up <- intersect(tripneg_up$Genes, doubleneg_up$Genes)
veryinteresting_up <- intersect(veryinteresting_up, singleneg_up$Genes)
intersect(tripneg_up$Genes, doubleneg_up$Genes)
intersect(tripneg_up$Genes, singleneg_up$Genes)
intersect(doubleneg_up$Genes, singleneg_up$Genes)

veryinteresting_down <- intersect(tripneg_down$Genes, doubleneg_down$Genes)
veryinteresting_down <- intersect(veryinteresting_down, singleneg_down$Genes)
intersect(tripneg_down$Genes, doubleneg_down$Genes)
intersect(tripneg_down$Genes, singleneg_down$Genes)
intersect(doubleneg_down$Genes, singleneg_down$Genes)

tripneg_up <- subset(tripneg_up, select = -X)
tripneg_up$Genes <- factor(tripneg_up$Genes, levels = tripneg_up$Genes)
m_tripneg_up <- melt(tripneg_up)

ggplot(m_tripneg_up, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(tripneg_up$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

tripneg_down <- subset(tripneg_down, select = -X)
tripneg_down$Genes <- factor(tripneg_down$Genes, levels = tripneg_down$Genes)
m_tripneg_down <- melt(tripneg_down)

ggplot(m_tripneg_down, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(tripneg_down$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

doubleneg_up <- subset(doubleneg_up, select = -X)
doubleneg_up$Genes <- factor(doubleneg_up$Genes, levels = doubleneg_up$Genes)
m_doubleneg_up <- melt(doubleneg_up)

ggplot(m_doubleneg_up, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(doubleneg_up$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

doubleneg_down <- subset(doubleneg_down, select = -X)
doubleneg_down$Genes <- factor(doubleneg_down$Genes, levels = doubleneg_down$Genes)
m_doubleneg_down <- melt(doubleneg_down)

ggplot(m_doubleneg_down, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(doubleneg_down$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

singleneg_up <- subset(singleneg_up, select = -X)
singleneg_up$Genes <- factor(singleneg_up$Genes, levels = singleneg_up$Genes)
m_singleneg_up <- melt(singleneg_up)

ggplot(m_singleneg_up, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(singleneg_up$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

singleneg_down <- subset(singleneg_down, select = -X)
singleneg_down$Genes <- factor(singleneg_down$Genes, levels = singleneg_down$Genes)
m_singleneg_down <- melt(singleneg_down)

ggplot(m_singleneg_down, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(singleneg_down$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))


##### To find most likely/intersting/promising candidates - BY STAGE:
#MCF7/iPS Overlap - STAGE
st_tripneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/MCF7iPSOverlap/TripNeg_MCF7iPSOverlap_ByStage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
st_doubleneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/MCF7iPSOverlap/DoubleNeg_MCF7iPSOverlap_ByStage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
st_singleneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/MCF7iPSOverlap/Her2Pos_MCF7iPSOverlap_ByStage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
#Total iPS Overlap - STAGE
st_tripneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/TotaliPSOverlap/TripNeg_TotaliPSOverlap_ByStage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
st_doubleneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/TotaliPSOverlap/DoubleNeg_TotaliPSOverlap_ByStage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
st_singleneg <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/MeanExpressionChanges/TotaliPSOverlap/Her2Pos_TotaliPSOverlap_ByStage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

st_tripneg_up <- filter(st_tripneg, Stage2 > 0 & Stage3 > Stage2)
st_tripneg_down <- filter(st_tripneg, Stage2 < 0 & Stage3 < Stage2)

st_tripneg_up$DeltaS3toS2 <- st_tripneg_up$Stage3 - st_tripneg_up$Stage2
st_tripneg_up <- arrange(st_tripneg_up, -DeltaS3toS2, -Stage2)

st_tripneg_down$DeltaS3toS2 <- st_tripneg_down$Stage3 - st_tripneg_down$Stage2
st_tripneg_down <- arrange(st_tripneg_down, DeltaS3toS2, Stage2)

st_doubleneg_up <- filter(st_doubleneg, Stage2 > 0 & Stage3 > Stage2 & Stage4 > Stage3)
st_doubleneg_down <- filter(st_doubleneg, Stage2 < 0 & Stage3 < Stage2 & Stage4 < Stage3)

st_doubleneg_up$DeltaS3toS2 <- st_doubleneg_up$Stage3 - st_doubleneg_up$Stage2
st_doubleneg_up$DeltaS4toS3 <- st_doubleneg_up$Stage4 - st_doubleneg_up$Stage3
st_doubleneg_up <- arrange(st_doubleneg_up, -DeltaS4toS3, -DeltaS3toS2, -Stage2)

st_doubleneg_down$DeltaS3toS2 <- st_doubleneg_down$Stage3 - st_doubleneg_down$Stage2
st_doubleneg_down$DeltaS4toS3 <- st_doubleneg_down$Stage4 - st_doubleneg_down$Stage3
st_doubleneg_down <- arrange(st_doubleneg_down, -DeltaS4toS3, -DeltaS3toS2, -Stage2)

st_singleneg_up <- filter(st_singleneg, Stage2 > 0 & Stage3 > Stage2)
st_singleneg_down <- filter(st_singleneg, Stage2 < 0 & Stage3 < Stage2)

st_singleneg_up$DeltaS3toS2 <- st_singleneg_up$Stage3 - st_singleneg_up$Stage2
st_singleneg_up <- arrange(st_singleneg_up, -DeltaS3toS2, -Stage2)

st_singleneg_down$DeltaS3toS2 <- st_singleneg_down$Stage3 - st_singleneg_down$Stage2
st_singleneg_down <- arrange(st_singleneg_down, DeltaS3toS2, Stage2)

st_veryinteresting_up <- intersect(st_tripneg_up$Genes, st_doubleneg_up$Genes)
st_veryinteresting_up <- intersect(st_veryinteresting_up, st_singleneg_up$Genes)
intersect(st_tripneg_up$Genes, st_doubleneg_up$Genes)
intersect(st_tripneg_up$Genes, st_singleneg_up$Genes)
intersect(st_doubleneg_up$Genes, st_singleneg_up$Genes)

st_veryinteresting_down <- intersect(st_tripneg_down$Genes, st_doubleneg_down$Genes)
st_veryinteresting_down <- intersect(st_veryinteresting_down, st_singleneg_down$Genes)
intersect(st_tripneg_down$Genes, st_doubleneg_down$Genes)
intersect(st_tripneg_down$Genes, st_singleneg_down$Genes)
intersect(st_doubleneg_down$Genes, st_singleneg_down$Genes)

st_tripneg_up <- subset(st_tripneg_up, select = -X)
st_tripneg_up$Genes <- factor(st_tripneg_up$Genes, levels = st_tripneg_up$Genes)
m_st_tripneg_up <- melt(st_tripneg_up)

ggplot(m_st_tripneg_up, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(st_tripneg_up$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

st_tripneg_down <- subset(st_tripneg_down, select = -X)
st_tripneg_down$Genes <- factor(st_tripneg_down$Genes, levels = st_tripneg_down$Genes)
m_st_tripneg_down <- melt(st_tripneg_down)

ggplot(m_st_tripneg_down, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(st_tripneg_down$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

st_doubleneg_up <- subset(st_doubleneg_up, select = -X)
st_doubleneg_up$Genes <- factor(st_doubleneg_up$Genes, levels = st_doubleneg_up$Genes)
m_st_doubleneg_up <- melt(st_doubleneg_up)

ggplot(m_st_doubleneg_up, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(st_doubleneg_up$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

st_doubleneg_down <- subset(st_doubleneg_down, select = -X)
st_doubleneg_down$Genes <- factor(st_doubleneg_down$Genes, levels = st_doubleneg_down$Genes)
m_st_doubleneg_down <- melt(st_doubleneg_down)

ggplot(m_st_doubleneg_down, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(st_doubleneg_down$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

st_singleneg_up <- subset(st_singleneg_up, select = -X)
st_singleneg_up$Genes <- factor(st_singleneg_up$Genes, levels = st_singleneg_up$Genes)
m_st_singleneg_up <- melt(st_singleneg_up)

ggplot(m_st_singleneg_up, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(st_singleneg_up$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))

st_singleneg_down <- subset(st_singleneg_down, select = -X)
st_singleneg_down$Genes <- factor(st_singleneg_down$Genes, levels = st_singleneg_down$Genes)
m_st_singleneg_down <- melt(st_singleneg_down)

ggplot(m_st_singleneg_down, aes(x=variable, y=Genes)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="Genes", limits = rev(levels(st_singleneg_down$Genes))) + 
  geom_text(aes(x=variable, y=Genes, label = round(value, digits = 2)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))




intersect(st_veryinteresting_down, veryinteresting_down)
intersect(st_veryinteresting_up, veryinteresting_up)
















#tripneg_up <- filter(tripneg, Grade2 >= 0.2)
#tripneg_down <- filter(tripneg, Grade2 <= -0.2)
#tripneg_up <- filter(tripneg_up, Grade3 >= (Grade2*1.25))
#tripneg_down <- filter(tripneg_down, Grade3 <= (Grade2*1.25))

#doubleneg_up <- filter(doubleneg, Grade2 >= 0.2)
#doubleneg_down <- filter(doubleneg, Grade2 <= -0.2)
#doubleneg_up <- filter(doubleneg_up, Grade3 >= (Grade2*1.25))
#doubleneg_down <- filter(doubleneg_down, Grade3 <= (Grade2*1.25))

#singleneg_up <- filter(singleneg, Grade2 >= 0.2)
#singleneg_down <- filter(singleneg, Grade2 <= -0.2)
#singleneg_up <- filter(singleneg_up, Grade3 >= (Grade2*1.25))
#singleneg_down <- filter(singleneg_down, Grade3 <= (Grade2*1.25))
