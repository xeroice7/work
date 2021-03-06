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
patientdata <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Uterus/ucec/tcga/data_bcr_clinical_data_patient.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

colnames(patientdata) <- lapply(patientdata[4,], as.character)

patientdata <- patientdata[-(1:4),]

patientdata <- subset(patientdata, select = c(PATIENT_ID, CLINICAL_STAGE, GRADE, OS_STATUS)) 

splits <- str_split_fixed(patientdata$CLINICAL_STAGE, " ", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '2'] <- 'Specific_Stages'
patientdata <- subset(patientdata, select=-c(1, CLINICAL_STAGE))

splits <- str_split_fixed(patientdata$Specific_Stage, "A|B|C", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'General_Stages'
patientdata <- subset(patientdata, select=-2)


#### READ IN CLINICAL DATA AND TIDY
expressiondata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Uterus/ucec/tcga/data_expression_median.csv", sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

expressiondata <- expressiondata[-1,]

expressiondata[expressiondata == "null"] <- NA

expressiondata <- subset(expressiondata, select = -Entrez_Gene_Id)

expressiondata <- expressiondata %>% drop_na()

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

#### FILTERING DATA WITH OUR CRITERIA OF INTEREST
grades <- factor(c("G1", "G2", "G3", "High Grade"))
grades <- ordered(grades, levels = c("G1", "G2", "G3", "High Grade"))
spec.stages <- factor(c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IIIC", "IIIC1", "IIIC2", "IV", "IVA", "IVB"))
spec.stages <- ordered(spec.stages, levels = c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IIIC", "IIIC1", "IIIC2", "IV", "IVA", "IVB"))
gen.stages <- factor(c("I", "II", "III", "IV"))
gen.stages <- ordered(gen.stages, levels = c("I", "II", "III", "IV"))
status <- factor(c("LIVING", "DECEASED"))
status <- ordered(status, levels = c("LIVING", "DECEASED"))

#Subset rows that are matches to those in the targets
stage_diffs <- patientDF
stage_diffs <- subset(stage_diffs, GENE %in% huvfips_overlap)
#stage_diffs <- subset(stage_diffs, OS_STATUS %in% status)
#stage_diffs <- subset(stage_diffs, GRADE %in% grades)
#stage_diffs <- subset(stage_diffs, General_Stages %in% gen.stages)
#stage_diffs <- subset(stage_diffs, Specific_Stages %in% spec.stages)
#stage_diffs <- filter(stage_diffs, ((Tumor == "T3") | (Tumor == "T2")))
stage_diffs <- filter(stage_diffs, General_Stages == "III")
#stage_diffs <- filter(stage_diffs, GRADE == "G3")

patient <- dcast(stage_diffs, PATIENT_ID+General_Stages+Specific_Stages+GRADE+OS_STATUS ~ GENE, value.var = "EXPRESSION_LEVEL")

extras <- cbind.data.frame(patient$PATIENT_ID, patient$Specific_Stages, patient$GRADE, patient$General_Stages, patient$OS_STATUS)

colnames(extras) <- c("PATIENT_ID", "Specific_Stages", "GRADE", "General_Stages", "OS_STATUS")

splits <- str_split_fixed(patient$PATIENT_ID, pattern = "-", "3")

patient <- cbind.data.frame(splits, patient)

patient <- subset(patient, select = -c(`1`, `2`, PATIENT_ID)) # Getting rid of Stage Column to prepare for a calculation matrix

names(patient)[names(patient) == '3'] <- 'ID'

patient <- arrange(patient, General_Stages, Specific_Stages, GRADE, OS_STATUS)

patient <- patient %>% 
  unite(PATIENT_ID, ID, General_Stages, Specific_Stages, GRADE, sep = "_", remove = TRUE)

patient <- subset(patient, select = -OS_STATUS)

patient[,1] <- as.character(patient[,1])

patient[,2:ncol(patient)] <- sapply(patient[,2:ncol(patient)],as.numeric)

sapply(patient, class)  #to check classes

patient <- patient %>% remove_rownames %>% column_to_rownames(var="PATIENT_ID") # Making the patient IDs the rownames, not the first column

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
#heatmap(test1, Colv=NA, col=greenred(10),scale="none")
#cor(t(test1))
t <- 1-cor(t(test1))
hc <- hclust(as.dist(1-cor(t(test1))))
plot(hc)
heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=greenred(10), cexRow = 0.2)
heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=redblackgreen, cexRow = 0.2)
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
#######
t_patient <- t(patient)
mean <- rowMeans(t_patient, na.rm = FALSE)
t_patient <- cbind.data.frame(t_patient, mean)

t_patient_lo <- t_patient #For Low
t_patient_2 <- t_patient
t_patient_3 <- t_patient
t_patient_hi <- t_patient #For Hi 

colnames(t_patient_1)[colnames(t_patient_1) == 'mean'] <- 'mean_1'
colnames(t_patient_2)[colnames(t_patient_2) == 'mean'] <- 'mean_2'
colnames(t_patient_3)[colnames(t_patient_3) == 'mean'] <- 'mean_3'
colnames(t_patient_4)[colnames(t_patient_4) == 'mean'] <- 'mean_4'

t_patient_m1 <- merge(t_patient_1, t_patient_2, by = "row.names")
t_patient_m2 <- merge(t_patient_3, t_patient_4, by = "row.names")
t_patient_means <- merge(t_patient_m1, t_patient_m2, by = "Row.names")

t_patient_means$mean_diff_1 <- t_patient_means$mean_1-t_patient_means$mean_1
t_patient_means$mean_diff_2 <- t_patient_means$mean_2-t_patient_means$mean_1
t_patient_means$mean_diff_3 <- t_patient_means$mean_3-t_patient_means$mean_1
t_patient_means$mean_diff_4 <- t_patient_means$mean_4-t_patient_means$mean_1

cast_patient <- cbind.data.frame(t_patient_means$Row.names, t_patient_means$mean_diff_1, t_patient_means$mean_diff_2, t_patient_means$mean_diff_3, t_patient_means$mean_diff_4)
colnames(cast_patient) <- c("Genes", "Stage1", "Stage2", "Stage3", "Stage4")
melt_patient <- melt(cast_patient)
colnames(melt_patient) <- c("Genes", "Stages", "Difference")
melt_patient$Genes <- as.factor(melt_patient$Genes)
#plot_patient <- cbind.data.frame(t_patient$Row.names, t_patient$mean_diff)

ggplot(melt_patient) + 
  geom_tile(aes(x=Stages, y=rev(Genes), fill = Difference)) + 
  #scale_fill_gradient(low = "white", high = "firebrick", name = "") + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="", limits = rev(levels(melt_patient$Genes))) + 
  geom_text(aes(x=Stages, y=rev(Genes), label = round(Difference, digits = 1)), size=2) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8)) 

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

ggplot(plotcscDF) + 
  geom_tile(aes(x = Population, y = Target, fill = meannormalized)) +  
  scale_fill_gradient(low = "white", high = "royalblue", name = "") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 


heatmap.2(plot_patient, 
          Rowv=NA, 
          Colv=NA, 
          col=bluered(256), 
          #breaks = breaks,
          cexRow = 0.7, 
          cexCol = 0.7, 
          scale = "none", 
          trace = "none")





#heatmap(test1, Colv=NA, col=greenred(10),scale="none")
#cor(t(test1))
t <- 1-cor(t(test1))
hc <- hclust(as.dist(1-cor(t(test1))))
plot(hc)
heatmap(test1, Rowv=as.dendrogram(hc) , Colv=T, col=greenred(10), cexRow = 0.2)
heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=redblackgreen, cexRow = 0.2)

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

