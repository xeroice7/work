####
# cBio - Lung Adenocarcinoma TCGA (Provisional) - 522 Patients
####

#Necessary Packages =======================================
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ISLR)
library(gplots)
library(caret)
library(rpart)
library(randomForest)
library(mlbench)
library(rattle)
library(rpart.plot)
library(class)
library(gmodels)

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

#iPS to A549 Comparisons (No Somatic Filters - can include somatic hits, 37 genes)
a549huvips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/A549vsHUViPS/HUViPSA549OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
a549fips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/A549vsFiPS/FiPSA549OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
a549ips_overlap <- intersect(a549fips_overlap$V1, a549huvips_overlap$V1)

#iPS to H1299 Comparisons(No Somatic Filters - can include somatic hits, 28 genes)
ncihuvips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/NCIvsFiPS/FiPSNCIOverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
ncifips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/NCIvsHUViPS/HUViPSNCIOverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
nciips_overlap <- intersect(ncifips_overlap$V1, ncihuvips_overlap$V1)

#Load Patient Data ============================================
# RNA Seq Data (too few recordings for microarray expression)
patientdata <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/data_bcr_clinical_data_patient.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

colnames(patientdata) <- lapply(patientdata[4,], as.character)

patientdata <- patientdata[-(1:4),]

#Tidy Patient Data ============================================
# Split Stages
splits <- str_split_fixed(patientdata$AJCC_PATHOLOGIC_TUMOR_STAGE, " ", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '2'] <- 'Specific_Stage'
patientdata <- subset(patientdata, select=-1)

splits <- str_split_fixed(patientdata$Specific_Stage, "A|B|C", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'General_Stages'
patientdata <- subset(patientdata, select=-2)

# Split Tumor
splits <- str_split_fixed(patientdata$AJCC_TUMOR_PATHOLOGIC_PT, "a|b", 2)
patientdata <- cbind.data.frame(splits, patientdata)
names(patientdata)[names(patientdata) == '1'] <- 'General_Grade'
patientdata <- subset(patientdata, select=-2)
names(patientdata)[names(patientdata) == 'AJCC_TUMOR_PATHOLOGIC_PT'] <- 'Specific_Grade'

patientdata <- subset(patientdata, select=-OTHER_PATIENT_ID)

#Load Clinical Data =============================================
expressiondata <- read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/data_RNA_Seq_v2_expression_median_edited.csv", check.names=FALSE, sep = ",", stringsAsFactors = FALSE)

expressiondata <- subset(expressiondata, select = -Entrez_Gene_Id)

expressiondata <- expressiondata %>% drop_na()

#Tidy Clinical Data =============================================
#To check to see if there are duplicate patients or genes
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

#expression <- t(expressiondata)

#Melt data
expressionDF <- melt(expressiondata, id.vars = "Hugo_Symbol")

#Get rid of last "-X" from patient ID to match patient data
splits <- str_split_fixed(expressionDF$variable, "-", 4)
expressionDF <- cbind.data.frame(splits, expressionDF)
expressionDF <- expressionDF %>% 
  unite(ID, '1', '2', '3', sep = "-", remove = TRUE)
expressionDF <- subset(expressionDF, select=-c(`4`, variable))

#Rename columns
colnames(expressionDF) <- c("PATIENT_ID", "GENE", "EXPRESSION_LEVEL")

#Recasting Clinical Data ===============================================
expressionDF <- dcast(expressionDF, PATIENT_ID ~ GENE, value.var = "EXPRESSION_LEVEL")

#Merging Clinical and Patient Data ===================================== 
patientDF <- merge(x = expressionDF, y = patientdata, by = "PATIENT_ID")

#ML Algos ==============================================================
data <- patientDF
#stage <- filter(patientDF, General_Stages == "I" | General_Stages == "IV")
stageIV <- filter(data, General_Stages == "IV")
stageI <- filter(data, General_Stages == "I")
#stageI <- sample_n(stageI, 175)
stage300 <- rbind.data.frame(stageI, stageIV)
stage <- stage300[,c(1,50:10050,20533)]
stage$General_Stages <- ordered(stage$General_Stages, levels = c("I", "IV"))
#stage1 <- cbind.data.frame(stage$PATIENT_ID, stage$General_Stages)
#colnames(stage1) <- c("PATIENT_ID", "General_Stages")
gradeT1 <- filter(data, General_Grade == "T1")
gradeT3 <- filter(data, General_Grade == "T3")
gradeT4 <- filter(data, General_Grade == "T4")

gradeT1T3 <- rbind.data.frame(gradeT1, gradeT3)
grade <- gradeT1T3[,c(1,50:1050,20532)]
grade$General_Grade <- ordered(grade$General_Grade, levels = c("T1", "T3"))
#stage <- merge(x = stage, y = stage1, by = "PATIENT_ID")
#names(stage) <- gsub(x = names(stage), pattern = "-", replacement = "")  

#names.use <- names(data)[(names(data) %in% gene)]
#names.use <- names(data)[(names(data) %in% a549ips_overlap)]
names.use <- names(data)[(names(data) %in% totaltestgenes)]
len <- length(names.use)

#data.ips <- data[, c("PATIENT_ID", "General_Grade", names.use)]
#data.a549 <- data[, c("PATIENT_ID", "General_Grade", names.use)]
data.genes <- data[, c("PATIENT_ID", "General_Grade", names.use)]
#data.rand <- data[, sample(ncol(data), len)]
#data.rand <- cbind.data.frame(data$PATIENT_ID, data$General_Grade, data.rand)

randomdata.a549 <- data[, sample(ncol(data), len)]
randomdata.a549 <- cbind.data.frame(data$PATIENT_ID, data$General_Grade, randomdata.a549)
names(randomdata.a549)[names(randomdata.a549) == 'data$PATIENT_ID'] <- 'PATIENT_ID'
names(randomdata.a549)[names(randomdata.a549) == 'data$General_Grade'] <- 'General_Grade'

#names(data.rand)[names(data.rand) == 'data$PATIENT_ID'] <- 'PATIENT_ID'
#names(data.rand)[names(data.rand) == 'data$General_Grade'] <- 'General_Grade'
#data.ips <- filter(data.ips, General_Grade == "T1" | General_Grade == "T3")
#data.ips$General_Grade <- ordered(data.ips$General_Grade, levels = c("T1", "T3"))
data.genes <- filter(data.genes, General_Grade == "T1" | General_Grade == "T3")
data.genes$General_Grade <- ordered(data.genes$General_Grade, levels = c("T1", "T3"))
#data.a549 <- filter(data.a549, General_Grade == "T1" | General_Grade == "T3")
#data.a549$General_Grade <- ordered(data.a549$General_Grade, levels = c("T1", "T3"))
#randomdata.a549 <- filter(randomdata.a549, General_Grade == "T1" | General_Grade == "T3")
#randomdata.a549$General_Grade <- ordered(randomdata.a549$General_Grade, levels = c("T1", "T3"))

#data.rand <- filter(data.rand, General_Grade == "T1" | General_Grade == "T3")
#data.rand$General_Grade <- ordered(data.rand$General_Grade, levels = c("T1", "T3"))

# Set a random seed
set.seed(754)

#Find most important factors
#control <- trainControl(method="repeatedcv", number=2, repeats=2)
#fit <- randomForest(General_Stages~., data=stage)
#fit <- rpart(General_Stages ~ ., data = stage, method = "class")

# train the model
#model <- train(General_Stages~., data=stage, method="lvq", preProcess="scale", trControl=control)

# estimate variable importance
#importance <- varImp(model, scale=FALSE)
# summarize importance
#print(importance)
# plot importance
#plot(importance)

# Create index to split based on labels  
#index.ips <- createDataPartition(data.ips$General_Grade, p=0.75, list=FALSE)
#index.rand <- createDataPartition(data.rand$General_Grade, p=0.75, list=FALSE)
index.a549 <- createDataPartition(data.a549$General_Grade, p=0.75, list = FALSE)
index.randomdata.a549 <- createDataPartition(randomdata.a549$General_Grade, p=0.75, list = FALSE)
index.genes <- createDataPartition(data.genes$General_Grade, p=0.75, list=FALSE)
# Subset training set with index
#ips.training <- data.ips[index.ips,]
#rand.training <- data.rand[index.rand,]
a549.training <- data.a549[index.a549,]
randomdata.a549.training <- randomdata.a549[index.randomdata.a549,]
genes.train <- data.genes[index.genes,]
# Subset test set with index
#ips.test <- data.ips[-index.ips,]
#rand.test <- data.rand[-index.rand,]
genes.test <- data.genes[-index.genes,]
a549.test <- data.a549[-index.a549,]
randomdata.a549.test <- randomdata.a549[-index.randomdata.a549,]

# Overview of algos supported by caret
#names(getModelInfo())

# Train a model
#model_knn <- train(iris.training[, 1:4], iris.training[, 5], method='knn')

#Possible to make other models simply by changing the method argument
#model_cart <- train(iris.training[, 1:4], iris.training[, 5], method='rpart2')

# Predict the labels of the test set
#predictions<-predict(object=model_knn,iris.test[,1:4])

# Evaluate the predictions
#table(predictions)

# Confusion matrix 
#confusionMatrix(predictions,iris.test[,5]) # Allows us to see the stats and results behind the predictions

# Now lets try preprocessing data with scaling and centering 
# Train the model with preprocessing
#model_knn <- train(iris.training[, 1:4], iris.training[, 5], method='knn', preProcess=c("center", "scale"))

# Predict values
#predictions<-predict.train(object=model_knn,iris.test[,1:4], type="raw")

# Confusion matrix
#confusionMatrix(predictions,iris.test[,5])

# Random Forest
#set.seed(7)
#fit.rf <- train(Species~., data=dataset, method="rf", metric=metric, trControl=control)

# Build the model (note: not all possible variables are used)
#stage.training$General_Stages <- ordered(stage.training$General_Stages, levels = c("I", "IV"))

#model.rf <- train(General_Stages ~ ., method = "rf", data = stage.training)

#modelgrade.rf <- train(stage.training[,2:1002], stage.training[,1003], method = "rf")
#ips.model.rf <- train(ips.training[,3:24], ips.training[,2], method = "rpart")
#rand.model.rf <- train(rand.training[,3:24], rand.training[,2], method = "rpart")
a549.model.rf <- train(a549.training[,3:24], a549.training[,2], method = "rf")
randomdata.a549.model.rf <- train(randomdata.a549.training[,3:24], randomdata.a549.training[,2], method = "rf")
genes.model.rf <- train(genes.train[3:289], genes.train[,2], method = "rf")
#predictgrade.rf <- predict(object=modelgrade.rf, stage.test)
#ips.predict.rf <- predict(object=ips.model.rf, ips.test)
#rand.predict.rf <- predict(object=rand.model.rf, rand.test)
a549.predict.rf <- predict(object=a549.model.rf, a549.test)
randomdata.a549.predict.rf <- predict(object=randomdata.a549.model.rf, randomdata.a549.test)
genes.predict.rf <- predict(genes.model.rf, genes.test)

#Possible to make other models simply by changing the method argument
#model_cart <- train(iris.training[, 1:4], iris.training[, 5], method='rpart2')

# Predict the labels of the test set
#predictions<-predict(object=model_knn,iris.test[,1:4])

# Evaluate the predictions
#table(predictgrade.rf)
table(ips.predict.rf)
table(rand.predict.rf)

# Confusion matrix 
#confusionMatrix(predictgrade.rf, stage.test[,1003]) # Allows us to see the stats and results behind the predictions
a <- confusionMatrix(ips.predict.rf, ips.test[,2])
b <- confusionMatrix(rand.predict.rf, rand.test[,2])
x <- confusionMatrix(a549.predict.rf, a549.test[,2])
y <- confusionMatrix(randomdata.a549.predict.rf, randomdata.a549.test[,2])

tocsv_x <- cbind.data.frame(t(x$positive), t(x$table), t(x$overall), t(x$byClass), t(x$dots))
tocsv_y <- cbind.data.frame(t(y$positive), t(y$table), t(y$overall), t(y$byClass), t(y$dots))
write.csv(tocsv_x, "Lung_LUAD_A549Test_A549.csv")
write.csv(tocsv_y, "Lung_LUAD_A549Test_Random.csv")

tocsv_a <- cbind.data.frame(t(a$positive), t(a$table), t(a$overall), t(a$byClass), t(a$dots))
tocsv_b <- cbind.data.frame(t(b$positive), t(b$table), t(b$overall), t(b$byClass), t(b$dots))
write.csv(tocsv_a, "Lung_LUAD_iPSTest_iPS.csv")
write.csv(tocsv_b, "Lung_LUAD_iPSTest_Random.csv")

confusionMatrix(genes.predict.rf, genes.test[,2])

genes.train2 <- genes.train[,2:ncol(genes.train)]
names(genes.train2)[names(genes.train2) == '2-Sep'] <- 'Sep2'
genes.test2 <- genes.test[,2:ncol(genes.test)]
names(genes.test2)[names(genes.test2) == '2-Sep'] <- 'Sep2'
model.genes.rf <- randomForest(General_Grade ~., data = genes.train2, importance = TRUE)
#ips.test2 <- ips.test[,2:ncol(ips.test)]
#ips.train <- ips.training[,2:ncol(ips.training)]
#model.ips <- randomForest(General_Grade ~., data = ips.train, importance = TRUE)
#model.ips.rpart <- rpart(General_Grade ~., data = ips.train, method = "class")
a549.train2 <- a549.training[,2:ncol(a549.training)]
a549.test2 <- a549.test[,2:ncol(a549.test)]
model.a549 <- randomForest(General_Grade ~., data = a549.train2, importance = TRUE)

randomdata.a549.train2 <- randomdata.a549.training[,2:ncol(randomdata.a549.training)]
names(randomdata.a549.train2)[names(randomdata.a549.train2) == 'SNORD114-29'] <- 'SNORD11429'
randomdata.a549.test2 <- randomdata.a549.test[,2:ncol(randomdata.a549.test)]
names(randomdata.a549.test2)[names(randomdata.a549.test2) == 'SNORD114-29'] <- 'SNORD11429'
model.randomdata.a549 <- randomForest(General_Grade ~., data = randomdata.a549.train2, importance = TRUE)

#rand.test2 <- rand.test[,2:ncol(rand.test)]
#names(rand.test2)[names(rand.test2) == 'SNORD114-5'] <- 'SNORD1145'
#rand.train <- rand.training[,2:ncol(rand.training)]
#names(rand.train)[names(rand.train) == 'SNORD114-5'] <- 'SNORD1145'
#model.rand <- randomForest(General_Grade ~., data = rand.train, importance = TRUE)
#model.rand.rpart <- rpart(General_Grade ~., data = rand.train, method = "class")

#predict.ips <- predict(model.ips, ips.test2)
#predict.rand <- predict(model.rand, rand.test2)
predict.a549 <- predict(model.a549, a549.test2)
predict.randomdata.a549 <- predict(model.randomdata.a549, randomdata.a549.test2)
predict.genes <- predict(model.genes.rf, genes.test2)
#fancyRpartPlot(model.ips.rpart)

#fit <- rpart(General_Stages ~ ., data = stage, method = "class")

# Get importance
importance.genes <- varImp(genes.model.rf, scale = FALSE)
importance.genes.rf <- varImp(model.genes.rf, scale = FALSE)

#importance.ips <- varImp(model.ips, scale=FALSE)
#importance.rand <- varImp(model.rand, scale = FALSE)
importance.a549 <- varImp(model.a549, scale = FALSE)
importance.randomdata.a549 <- varImp(model.randomdata.a549, scale = FALSE)
print(importance.a549)
plot(importance.a549)
print(importance.randomdata.a549)
plot(importance.randomdata.a549)
print(importance.genes.rf)
plot(importance.genes.rf)

importance.ips <- importance(model.ips)
importance.rand <- importance(model.rand)
importance.tgenes <- importance(model.genes.rf)
importance.a549.2 <- importance(model.a549)
importance.randomdata.a549.2 <- importance(model.randomdata.a549)

# summarize and plot importance
print(importance.ips)
plot(importance.ips)
print(importance.rand)
plot(importance.rand)
plot(importance.tgenes)

varImportance.genes <- data.frame(Variables = row.names(importance.tgenes), 
                            Importance = round(importance.tgenes[ ,'MeanDecreaseGini'],2))

# Create a rank variable based on importance
rankImportance.genes <- varImportance.genes %>%
  mutate(Rank = paste0('#',dense_rank(desc(Importance))))

# Use ggplot2 to visualize the relative importance of variables
ggplot(rankImportance.genes, aes(x = reorder(Variables, Importance), 
                           y = Importance, fill = Importance)) +
  geom_bar(stat='identity') + 
  geom_text(aes(x = Variables, y = 0.5, label = Rank),
            hjust=0, vjust=0.55, size = 4, colour = 'darkblue') +
  labs(x = 'Variables') +
  theme(axis.text.y = element_text(size = 6, color = ifelse(rankImportance.genes$Variables %in% ipsgene, "red", "black"))) + 
  coord_flip() #+ 
  #theme_few()















#Additional Filter Parameters ==========================================
meta <- factor(c("M0", "M1"))
meta <- ordered(meta, levels = c("M0", "M1"))
grades <- factor(c("T1", "T1a", "T1b", "T2", "T2a", "T2b", "T3"))
grades <- ordered(grades, levels = c("T1", "T1a", "T1b", "T2", "T2a", "T2b", "T3"))
stages <- factor(c("I", "IA", "IB", "II", "IIA", "IIB", "IIIA", "IIIB", "IV"))
stages <- ordered(stages, levels = c("I", "IA", "IB", "II", "IIA", "IIB", "IIIA", "IIIB", "IV"))

stage_diffs <- patientDF
stage_diffs <- subset(stage_diffs, GENE %in% nciips_overlap) #Subset rows that are matches to those in the targets
#stage_diffs <- subset(stage_diffs, GENE %in% surface_genes)
#stage_diffs <- subset(stage_diffs, AJCC_METASTASIS_PATHOLOGIC_PM %in% meta)
#stage_diffs <- subset(stage_diffs, AJCC_TUMOR_PATHOLOGIC_PT %in% grades)
#stage_diffs <- subset(stage_diffs, Specific_Stage %in% stages)
#stage_diffs <- filter(stage_diffs, ((Tumor == "T3") | (Tumor == "T2")))
#stage_diffs <- filter(stage_diffs, General_Grade == "T4")
stage_diffs <- filter(stage_diffs, General_Stages == "IV")
#stage_diffs <- filter(stage_diffs, AJCC_METASTASIS_PATHOLOGIC_PM == "M0")

#Recast patient data ===================================================
patient <- dcast(stage_diffs, PATIENT_ID+General_Grade+Specific_Grade+General_Stages+Specific_Stage+AJCC_METASTASIS_PATHOLOGIC_PM ~ GENE, value.var = "EXPRESSION_LEVEL")

#Combine Patients with Clinical Data in Names ==========================
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
y <- cor(t(grade1), t(grade4))
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
results$Stage <- "Stage4"

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
t_patient_m1 <- merge(t_patient_lo, t_patient_hi, by = "row.names")
t_patient_means <- t_patient_m1

t_patient_means$mean_diff_1 <- t_patient_means$mean_1-t_patient_means$mean_1
t_patient_means$mean_diff_2 <- t_patient_means$mean_2-t_patient_means$mean_1
t_patient_means$mean_diff_3 <- t_patient_means$mean_3-t_patient_means$mean_1
t_patient_means$mean_diff_4 <- t_patient_means$mean_4-t_patient_means$mean_1

cast_patient <- cbind.data.frame(t_patient_means$Row.names, t_patient_means$mean_diff_1, t_patient_means$mean_diff_4)#, t_patient_means$mean_diff_3, t_patient_means$mean_diff_4)
colnames(cast_patient) <- c("Genes", "Stage1", "Stage2", "Stage3", "Stage4")
#colnames(cast_patient) <- c("Genes", "Grade1", "Grade2", "Grade3", "Grade4")
colnames(cast_patient) <- c("Genes", "Stage1", "Stage4")
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
melt_patient <- filter(melt_patient, Stages == "Stage4")
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
