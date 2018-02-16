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
library(caret)
library(rpart)
library(randomForest)
library(mlbench)
library(rattle)
library(rpart.plot)
library(class)
library(gmodels)
library(pROC)
library(plotROC)
library(ROCR)
library(MLmetrics)

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

expressionDF <- t(expressiondata) #Only use for ML - need this arrangement

expressionDF <- melt(expressiondata, id.vars = "Hugo_Symbol")

colnames(expressionDF) <- c("GENE", "PATIENT_ID", "EXPRESSION_LEVEL")

#Recasting Clinical Data ===============================================
expressionDF <- dcast(expressionDF, PATIENT_ID ~ GENE, value.var = "EXPRESSION_LEVEL") # Only use for ML - need this arrangement

#Merge Patient and Clinical Data ==================================== 
patientDF <- merge(x = expressionDF, y = patientDF, by = "PATIENT_ID")

#ML Algos ==============================================================
data <- patientDF
#Stage Filters
data$TUMOR_STAGE <- as.factor(data$TUMOR_STAGE)
#data$TUMOR_STAGE <- ordered(data$TUMOR_STAGE, c(1, 2, 3, 4) )
stage1 <- filter(data, TUMOR_STAGE == 1)
stage2 <- filter(data, TUMOR_STAGE == 2)
stage3 <- filter(data, TUMOR_STAGE == 3)
stage4 <- filter(data, TUMOR_STAGE == 4)
#stage1$Stage <- "Low"
#stage2$Stage <- "Low"
#stage3$Stage <- "High"
#stage4$Stage <- "High"
stage1$Stage <- "Low"
stage2$Stage <- "Low"
stage3$Stage <- "High"
stage4$Stage <- "High"
stageDF <- rbind.data.frame(stage1, stage2, stage3, stage4)
stageDF$Stage <- as.factor(stageDF$Stage)
stageDF$Stage <- ordered(stageDF$Stage, c("Low", "High"))
#stageDF$Stage <- ordered(stageDF$Stage, c("Low", "Mid", "High"))

#stage <- stage300[,c(1,50:10050,20533)]
#stage1 <- cbind.data.frame(stage$PATIENT_ID, stage$General_Stages)
#colnames(stage1) <- c("PATIENT_ID", "General_Stages")

#stage <- merge(x = stage, y = stage1, by = "PATIENT_ID")
#names(stage) <- gsub(x = names(stage), pattern = "-", replacement = "")  

names.use <- names(stageDF)[(names(stageDF) %in% gene)]
len <- length(names.use)
#rand <- sample(ncol(stageDF), len)
#data.mcf7 <- stageDF[, c("PATIENT_ID", "TUMOR_STAGE", "Stage", names.use)]
data.mcf7 <- stageDF[, c(names.use, "Stage")]

#USE FOR RANDOM Data
data.mcf7 <- stageDF[, sample(ncol(stageDF), len)]
rand1 <- colnames(data.mcf7)
rand2 <- colnames(data.mcf7)
rand3 <- colnames(data.mcf7)
rand4 <- colnames(data.mcf7)
rand5 <- colnames(data.mcf7)
data.mcf7 <- cbind(data.mcf7, stageDF$Stage)
colnames(data.mcf7)[colnames(data.mcf7) == 'stageDF$Stage'] <- 'Stage'

#Grade Filters
data$GRADE <- as.factor(data$GRADE)
grade1 <- filter(data, GRADE == 1)
grade2 <- filter(data, GRADE == 2)
grade3 <- filter(data, GRADE == 3)
grade1$tumor_grade <- "Low"
grade2$tumor_grade <- "Mid"
grade3$tumor_grade <- "High"
gradeDF <- rbind.data.frame(grade1, grade2, grade3)
gradeDF$tumor_grade <- as.factor(gradeDF$tumor_grade)
gradeDF$tumor_grade <- ordered(gradeDF$tumor_grade, c("Low", "Mid", "High"))
gradeDF <- filter(gradeDF, ER_STATUS == "+")
gradeDF <- filter(gradeDF, PR_STATUS == "+")
gradeDF <- filter(gradeDF, HER2_STATUS == "-")


names.use <- names(gradeDF)[(names(gradeDF) %in% mcf7ips_overlap)]
len <- length(names.use)
#rand <- sample(ncol(stageDF), len)
#data.mcf7 <- stageDF[, c("PATIENT_ID", "TUMOR_STAGE", "Stage", names.use)]
data.mcf7 <- gradeDF[, c(names.use, "tumor_grade")]

#USE FOR RANDOM Data
data.mcf7 <- gradeDF[, sample(ncol(gradeDF), len)]
#randGR1 <- colnames(data.mcf7)
#randGR2 <- colnames(data.mcf7)
#randGR3 <- colnames(data.mcf7)
#randGR4 <- colnames(data.mcf7)
#randGR5 <- colnames(data.mcf7)
data.mcf7 <- cbind(data.mcf7, gradeDF$tumor_grade)
colnames(data.mcf7)[colnames(data.mcf7) == 'gradeDF$tumor_grade'] <- 'tumor_grade'



#data.rand <- data[, sample(ncol(data), len)]
#data.rand <- cbind.data.frame(data$PATIENT_ID, data$General_Grade, data.rand)

#randomdata.a549 <- data[, sample(ncol(data), len)]
#randomdata.a549 <- cbind.data.frame(data$PATIENT_ID, data$General_Grade, randomdata.a549)
#names(randomdata.a549)[names(randomdata.a549) == 'data$PATIENT_ID'] <- 'PATIENT_ID'
#names(randomdata.a549)[names(randomdata.a549) == 'data$General_Grade'] <- 'General_Grade'

#names(data.rand)[names(data.rand) == 'data$PATIENT_ID'] <- 'PATIENT_ID'
#names(data.rand)[names(data.rand) == 'data$General_Grade'] <- 'General_Grade'
#data.ips <- filter(data.ips, General_Grade == "T1" | General_Grade == "T3")
#data.ips$General_Grade <- ordered(data.ips$General_Grade, levels = c("T1", "T3"))
#data.genes <- filter(data.genes, General_Grade == "T1" | General_Grade == "T3")
#data.genes$General_Grade <- ordered(data.genes$General_Grade, levels = c("T1", "T3"))
#data.a549 <- filter(data.a549, General_Grade == "T1" | General_Grade == "T3")
#data.a549$General_Grade <- ordered(data.a549$General_Grade, levels = c("T1", "T3"))
#randomdata.a549 <- filter(randomdata.a549, General_Grade == "T1" | General_Grade == "T3")
#randomdata.a549$General_Grade <- ordered(randomdata.a549$General_Grade, levels = c("T1", "T3"))

#data.rand <- filter(data.rand, General_Grade == "T1" | General_Grade == "T3")
#data.rand$General_Grade <- ordered(data.rand$General_Grade, levels = c("T1", "T3"))

#### Setting Model Parameters ======================================================
# Set a random seed
set.seed(754)

# Create index to split based on labels  
#index.mcf7 <- createDataPartition(data.mcf7$Stage, p=0.8, list = FALSE)
index.mcf <- createDataPartition(data.mcf7$tumor_grade, p=0.8, list = FALSE)

# Subset training set with index
mcf7.train<- data.mcf7[index.mcf7,]

# Subset test set with index
mcf7.test <- data.mcf7[-index.mcf7,]

# Creating and setting consistent training parameters for all models
ctrl <- trainControl(method = "repeatedcv",
                     number=10, 
                     repeats=10,
                     classProbs = TRUE,
                     summaryFunction = multiClassSummary,
                     verboseIter = TRUE,
                     savePredictions = TRUE,
                     sampling = "smote")

### Random Forest Model ===========================================================================
rf.model <- train(tumor_grade ~ ., data = mcf7.train,
                     method = "rf",
                     metric = "logLoss",
                     preProcess = c("scale", "center"),
                     na.action = na.omit,
                     trControl = ctrl)
rf.model$finalModel$confusion
rf.model$finalModel$importance
rfPredict <- predict(rf.model, newdata = mcf7.test)
rfOutput <- confusionMatrix(rfPredict, mcf7.test$tumor_grade)

selectedIndices <- rf.model$pred$mtry == 2

# Plot:
#plot.roc(rfFit$pred$obs[selectedIndices],
 #        rfFit$pred$M[selectedIndices])

plot.roc(rf.model$pred$obs[selectedIndices],
         rf.model$pred$Low[selectedIndices])

multiclass.roc(rf.model$pred$obs[selectedIndices],
         rf.model$pred$Low[selectedIndices])

g <- ggplot(rf.model$pred[selectedIndices, ], aes(m = Low, d = factor(obs, levels = c("Low", "High")))) + 
  geom_roc(hjust = -0.4, vjust = 1.5) + coord_equal() + style_roc()

#Updated ROC curve plot
g <- ggplot(rf.model$pred[selectedIndices, ], aes(m = Low, d = factor(obs, levels = c("Low", "High")))) + 
  geom_roc(n.cuts=0) +  
  coord_equal() +
  style_roc() #+
  #geom_ribbon(alpha=0.2)

g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4)))

rf1 <- varImp(rf.model, scale = FALSE)
plot(rf1)
### RPart/Decision Tree Model ====================================================================================
rpart.model <- train(tumor_grade ~ ., data = mcf7.train,
                  method = "rpart",
                  metric = "ROC",
                  preProcess = c("scale", "center"),
                  trControl = ctrl)
#rpart.model$finalModel$confusion
#rpart.model$finalModel$importance
rpartPredict <- predict(rpart.model, newdata = mcf7.test)
rpartOutput <- confusionMatrix(rpartPredict, mcf7.test$tumor_grade)


#Plot 1
rpart.pred <- predict(rpart.model, newdata = mcf7.test, type = "prob")[,2]
rpart.pred2 <- prediction(rpart.pred, mcf7.test$tumor_grade) 
roc <- performance(rpart.pred2, measure = "tpr", x.measure = "fpr")
plot(performance(rpart.pred2, "tpr", "fpr"))
plot(performance(rpart.pred2, "tpr", "fpr"), lwd=2,col="blue",
     main="ROC:  Classification Trees on Adult Dataset")
abline(0, 1, lty = 2)

#Plot 2
auc <- performance(rpart.pred2, measure = "auc")
auc <- auc@y.values[[1]]

roc.data <- data.frame(fpr=unlist(roc@x.values),
                       tpr=unlist(roc@y.values),
                       model="GLM")

ggplot(roc.data, aes(x=fpr, ymin=0, ymax=tpr)) +
  geom_ribbon(alpha=0.2) +
  geom_line(aes(y=tpr)) +
  ggtitle(paste0("ROC Curve w/ AUC=", auc)) +
  #geom_roc(n.cuts=0) +  
  coord_equal() +
  style_roc() + 
  annotate("text", x=0.75, y=0.25, label=paste("AUC =", round(auc, 4)))

rpart1 <- varImp(rpart.model, scale = FALSE)
plot(rpart1)
### K-Nearest Neighbors Model =================================================================
knn.model <- train(tumor_grade ~ ., data = mcf7.train,
                  method = "knn",
                  metric = "ROC",
                  preProcess = c("scale", "center"),
                  trControl = ctrl)

knnPredict <- predict(knn.model, newdata = mcf7.test)
knnOutput <- confusionMatrix(knnPredict, mcf7.test$tumor_grade)

#Plot 1
knn.pred <- predict(knn.model, newdata = mcf7.test, type = "prob")[,2]
knn.pred2 <- prediction(knn.pred, mcf7.test$tumor_grade) 
roc <- performance(knn.pred2, measure = "tpr", x.measure = "fpr")
plot(performance(knn.pred2, "tpr", "fpr"))
plot(performance(knn.pred2, "tpr", "fpr"), lwd=2,col="blue",
     main="ROC:  Classification Trees on Adult Dataset")
abline(0, 1, lty = 2)

#Plot 2
auc <- performance(knn.pred2, measure = "auc")
auc <- auc@y.values[[1]]

roc.data <- data.frame(fpr=unlist(roc@x.values),
                       tpr=unlist(roc@y.values),
                       model="knn")

ggplot(roc.data, aes(x=fpr, ymin=0, ymax=tpr)) +
  geom_ribbon(alpha=0.2) +
  geom_line(aes(y=tpr)) +
  ggtitle(paste0("ROC Curve w/ AUC=", auc)) +
  #geom_roc(n.cuts=0) +  
  coord_equal() +
  style_roc() + 
  annotate("text", x=0.75, y=0.25, label=paste("AUC =", round(auc, 4)))

knn1 <- varImp(knn.model, scale = FALSE)
plot(knn1)

### Logistic Regression model ==============================================================
logistic.model1 <- glm(tumor_grade ~ .,family=binomial(link='logit'), data = mcf7.train)
logistic.model2 <- train(tumor_grade ~ ., data = mcf7.train,
                         method = "glm",
                         metric = "ROC",
                         preProcess = c("scale", "center"),
                         trControl = ctrl)
logistic.predict1 <- predict(logistic.model1, newdata = mcf7.test)
logistic.predict2 <- predict(logistic.model2, newdata = mcf7.test)
logisticOutput1 <- confusionMatrix(logistic.predict1, mcf7.test$tumor_grade)
logisticOutput2 <- confusionMatrix(logistic.predict2, mcf7.test$tumor_grade)

fitted.results <- predict(logistic.model1, newdata = mcf7.test, type='response')
fitted.results <- ifelse(fitted.results > 0.5,"High","Low")
fitted.results <- as.factor(fitted.results)
fitted.results <- ordered(fitted.results, c("Low", "High"))
confusionMatrix(fitted.results, mcf7.test$tumor_grade)

misClasificError <- mean(fitted.results != mcf7.test$tumor_grade)
print(paste('Accuracy',1-misClasificError))

p <- predict(logistic.model1, newdata = mcf7.test, type="response")
pr <- prediction(p, mcf7.test$tumor_grade)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

roc.data <- data.frame(fpr=unlist(prf@x.values),
                       tpr=unlist(prf@y.values),
                       model="GLM")

ggplot(roc.data, aes(x=fpr, ymin=0, ymax=tpr)) +
  geom_ribbon(alpha=0.2) +
  geom_line(aes(y=tpr)) +
  ggtitle(paste0("ROC Curve w/ AUC=", auc)) +
  #geom_roc(n.cuts=0) +  
  coord_equal() +
  style_roc() + 
  annotate("text", x=0.75, y=0.25, label=paste("AUC =", round(auc, 4)))

lr1 <- varImp(logistic.model1, scale = FALSE)
lr1 <- data.frame(genes=rownames(lr1), lr1, row.names=NULL)
ggplot(lr1) + 
  geom_bar(aes(y=Overall, x=reorder(genes, Overall)), stat = "identity") + 
  xlab("Gene") + 
  ylab("Overall Importance") + 
  coord_flip()

#Plot for Model 2

logreg.pred <- predict(logistic.model2, newdata = mcf7.test, type = "prob")[,2]
logreg.pred2 <- prediction(logreg.pred, mcf7.test$tumor_grade) 
auc <- performance(logreg.pred2, measure = "auc")
auc <- auc@y.values[[1]]

roc.data <- data.frame(fpr=unlist(roc@x.values),
                       tpr=unlist(roc@y.values),
                       model="glm")

ggplot(roc.data, aes(x=fpr, ymin=0, ymax=tpr)) +
  geom_ribbon(alpha=0.2) +
  geom_line(aes(y=tpr)) +
  ggtitle(paste0("ROC Curve w/ AUC=", auc)) +
  #geom_roc(n.cuts=0) +  
  coord_equal() +
  style_roc() + 
  annotate("text", x=0.75, y=0.25, label=paste("AUC =", round(auc, 4)))


lr2 <- varImp(logistic.model2, scale = FALSE)
plot(lr2)

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
#index.mcf7 <- createDataPartition(data.mcf7$Stage, p=0.75, list = FALSE)

# Subset training set with index
#ips.training <- data.ips[index.ips,]
#rand.training <- data.rand[index.rand,]
#mcf7.train<- data.mcf7[index.mcf7,]
#randomdata.a549.training <- randomdata.a549[index.randomdata.a549,]
#genes.train <- data.genes[index.genes,]
# Subset test set with index
#ips.test <- data.ips[-index.ips,]
#rand.test <- data.rand[-index.rand,]
#genes.test <- data.genes[-index.genes,]
#mcf7.test <- data.mcf7[-index.mcf7,]
#randomdata.a549.test <- randomdata.a549[-index.randomdata.a549,]

ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     ## new option here:
                     sampling = "smote")

#set.seed(5627)
#down_inside <- train(Stage ~ ., data = mcf7.train,
#                     method = "treebag",
#                     nbagg = 50,
#                     metric = "ROC",
#                     trControl = ctrl)


set.seed(42)
model_rf <- caret::train(classes ~ .,
                         data = train_data,
                         method = "rf",
                         preProcess = c("scale", "center"),
                         trControl = trainControl(method = "repeatedcv", 
                                                  number = 10, 
                                                  repeats = 10, 
                                                  savePredictions = TRUE, 
                                                  verboseIter = FALSE))

model_rf$finalModel$confusion

imp <- model_rf$finalModel$importance
imp[order(imp, decreasing = TRUE), ]


#Subsampling for Class Imbalances
#DO the subsampling OUTSIDE of the train function
library(ROSE)
library(DMwR)
downtrain.mcf7 <- downSample(x = mcf7.train[,-ncol(mcf7.train)], y = mcf7.train$Stage)
table(downtrain.mcf7$Class)
uptrain.mcf7 <- upSample(x = mcf7.train[,-ncol(mcf7.train)], y = mcf7.train$Stage)
table(uptrain.mcf7$Class)
smotetrain.mcf7 <- SMOTE(Stage ~., data = mcf7.train)
table(smotetrain.mcf7$Stage)
#set.seed(9560)
#rosetrain.mcf7 <- ROSE(Stage ~ ., data  = mcf7.train)
#table(rosetrain.mcf7$Stage)
#rose_train <- ROSE(Class ~ ., data  = imbal_train)$data                        
#table(rose_train$Class) 
ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(5627)
orig_fit <- train(Stage ~ ., data = mcf7.train, 
                  method = "treebag",
                  nbagg = 50,
                  metric = "ROC",
                  trControl = ctrl)

set.seed(5627)
down_outside <- train(Class ~ ., data = downtrain.mcf7, 
                      method = "treebag",
                      nbagg = 50,
                      metric = "ROC",
                      trControl = ctrl)

set.seed(5627)
up_outside <- train(Class ~ ., data = uptrain.mcf7, 
                    method = "treebag",
                    nbagg = 50,
                    metric = "ROC",
                    trControl = ctrl)

set.seed(5627)
smote_outside <- train(Stage ~ ., data = smotetrain.mcf7, 
                      method = "treebag",
                      nbagg = 50,
                      metric = "ROC",
                      trControl = ctrl)



outside_models <- list(original = orig_fit,
                       down = down_outside,
                       up = up_outside,
                       SMOTE = smote_outside)

outside_resampling <- resamples(outside_models)

test_roc <- function(model, data) {
  #library(pROC)
  roc_obj <- roc(data$Stage, 
                 predict(model, data, type = "prob")[, "Low"],
                 levels = c("Low", "High"))
  ci(roc_obj)
}

outside_test <- lapply(outside_models, test_roc, data = mcf7.test)
outside_test <- lapply(outside_test, as.vector)
outside_test <- do.call("rbind", outside_test)
colnames(outside_test) <- c("lower", "ROC", "upper")
outside_test <- as.data.frame(outside_test)

summary(outside_resampling, metric = "ROC")

outside_test # Test and training sets for the area under the ROC curve do not appear to correlate.... 

#Subsampling during resampling - doing the subsampling during/inside the train function 
ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     ## new option here:
                     sampling = "down")

set.seed(5627)
down_inside <- train(Stage ~ ., data = mcf7.train,
                     method = "treebag",
                     nbagg = 50,
                     metric = "ROC",
                     trControl = ctrl)

## now just change that option
ctrl$sampling <- "up"

set.seed(5627)
up_inside <- train(Stage ~ ., data = mcf7.train,
                   method = "treebag",
                   nbagg = 50,
                   metric = "ROC",
                   trControl = ctrl)

#ctrl$sampling <- "rose"

#set.seed(5627)
#rose_inside <- train(Stage ~ ., data = mcf7.train,
                    # method = "treebag",
                     #nbagg = 50,
                     #metric = "ROC",
                     #trControl = ctrl)

ctrl$sampling <- "smote"

set.seed(5627)
smote_inside <- train(Stage ~ ., data = mcf7.train,
                      method = "treebag",
                      nbagg = 50,
                      metric = "ROC",
                      trControl = ctrl)


inside_models <- list(original = orig_fit,
                      down = down_inside,
                      up = up_inside,
                      SMOTE = smote_inside)

inside_resampling <- resamples(inside_models)

inside_test <- lapply(inside_models, test_roc, data = mcf7.test)
inside_test <- lapply(inside_test, as.vector)
inside_test <- do.call("rbind", inside_test)
colnames(inside_test) <- c("lower", "ROC", "upper")
inside_test <- as.data.frame(inside_test)

summary(inside_resampling, metric = "ROC")

inside_test


#Feature Plots to determine differences between Features (Finding the defining features)
library(AppliedPredictiveModeling)
transparentTheme(trans = .4)
library(caret)
featurePlot(x = iris[, 1:4], 
            y = iris$Species, 
            plot = "pairs",
            ## Add a key at the top
            auto.key = list(columns = 3))

###
featurePlot(x = mcf7[,20:ncol(mcf7)], 
            y = mcf7$Stage, 
            plot = "pairs",
            ## Add a key at the top
            auto.key = list(columns = 3))
###



featurePlot(x = mcf7[,25:ncol(mcf7)], 
            y = mcf7$Stage,  
            plot = "ellipse",
            ## Add a key at the top
            auto.key = list(columns = 3))

transparentTheme(trans = .9)
featurePlot(x = mcf7[,25:ncol(mcf7)], 
            y = mcf7$Stage, 
            plot = "density", 
            ## Pass in options to xyplot() to 
            ## make it prettier
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")), 
            adjust = 1.5, 
            pch = "|", 
            layout = c(4, 1), 
            auto.key = list(columns = 3))

featurePlot(x = iris[, 1:4], 
            y = iris$Species, 
            plot = "box", 
            ## Pass in options to bwplot() 
            scales = list(y = list(relation="free"),
                          x = list(rot = 90)),  
            layout = c(4,1 ), 
            auto.key = list(columns = 2))



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
mcf7.model.rf <- train(mcf7.train[,4:30], mcf7.train[,3], method = "rf")
mcf7.model.rf2 <- randomForest(Stage ~., data = mcf7.train[,3:ncol(mcf7.train)], importance = TRUE)
#randomdata.a549.model.rf <- train(randomdata.a549.training[,3:24], randomdata.a549.training[,2], method = "rf")
#genes.model.rf <- train(genes.train[3:289], genes.train[,2], method = "rf")
#predictgrade.rf <- predict(object=modelgrade.rf, stage.test)
#ips.predict.rf <- predict(object=ips.model.rf, ips.test)
#rand.predict.rf <- predict(object=rand.model.rf, rand.test)
mcf7.predict.rf <- predict(object=mcf7.model.rf, mcf7.test)
#randomdata.a549.predict.rf <- predict(object=randomdata.a549.model.rf, randomdata.a549.test)
#genes.predict.rf <- predict(genes.model.rf, genes.test)

#Possible to make other models simply by changing the method argument
#model_cart <- train(iris.training[, 1:4], iris.training[, 5], method='rpart2')

# Predict the labels of the test set
#predictions<-predict(object=model_knn,iris.test[,1:4])

# Evaluate the predictions
#table(predictgrade.rf)
table(mcf7.predict.rf)
table(rand.predict.rf)

# Confusion matrix 
#confusionMatrix(predictgrade.rf, stage.test[,1003]) # Allows us to see the stats and results behind the predictions
#a <- confusionMatrix(ips.predict.rf, ips.test[,2])
#b <- confusionMatrix(rand.predict.rf, rand.test[,2])
#x <- confusionMatrix(a549.predict.rf, a549.test[,2])
#y <- confusionMatrix(randomdata.a549.predict.rf, randomdata.a549.test[,2])

#tocsv_x <- cbind.data.frame(t(x$positive), t(x$table), t(x$overall), t(x$byClass), t(x$dots))
#tocsv_y <- cbind.data.frame(t(y$positive), t(y$table), t(y$overall), t(y$byClass), t(y$dots))
#write.csv(tocsv_x, "Lung_LUAD_A549Test_A549.csv")
#write.csv(tocsv_y, "Lung_LUAD_A549Test_Random.csv")

#tocsv_a <- cbind.data.frame(t(a$positive), t(a$table), t(a$overall), t(a$byClass), t(a$dots))
#tocsv_b <- cbind.data.frame(t(b$positive), t(b$table), t(b$overall), t(b$byClass), t(b$dots))
#write.csv(tocsv_a, "Lung_LUAD_iPSTest_iPS.csv")
#write.csv(tocsv_b, "Lung_LUAD_iPSTest_Random.csv")

confusionMatrix(mcf7.predict.rf, mcf7.test[,3])

#genes.train2 <- genes.train[,2:ncol(genes.train)]
#names(genes.train2)[names(genes.train2) == '2-Sep'] <- 'Sep2'
#genes.test2 <- genes.test[,2:ncol(genes.test)]
#names(genes.test2)[names(genes.test2) == '2-Sep'] <- 'Sep2'
#model.genes.rf <- randomForest(General_Grade ~., data = genes.train2, importance = TRUE)
#ips.test2 <- ips.test[,2:ncol(ips.test)]
#ips.train <- ips.training[,2:ncol(ips.training)]
#model.ips <- randomForest(General_Grade ~., data = ips.train, importance = TRUE)
#model.ips.rpart <- rpart(General_Grade ~., data = ips.train, method = "class")
#a549.train2 <- a549.training[,2:ncol(a549.training)]
#a549.test2 <- a549.test[,2:ncol(a549.test)]
#model.a549 <- randomForest(General_Grade ~., data = a549.train2, importance = TRUE)

#randomdata.a549.train2 <- randomdata.a549.training[,2:ncol(randomdata.a549.training)]
#names(randomdata.a549.train2)[names(randomdata.a549.train2) == 'SNORD114-29'] <- 'SNORD11429'
#randomdata.a549.test2 <- randomdata.a549.test[,2:ncol(randomdata.a549.test)]
#names(randomdata.a549.test2)[names(randomdata.a549.test2) == 'SNORD114-29'] <- 'SNORD11429'
#model.randomdata.a549 <- randomForest(General_Grade ~., data = randomdata.a549.train2, importance = TRUE)

#rand.test2 <- rand.test[,2:ncol(rand.test)]
#names(rand.test2)[names(rand.test2) == 'SNORD114-5'] <- 'SNORD1145'
#rand.train <- rand.training[,2:ncol(rand.training)]
#names(rand.train)[names(rand.train) == 'SNORD114-5'] <- 'SNORD1145'
#model.rand <- randomForest(General_Grade ~., data = rand.train, importance = TRUE)
#model.rand.rpart <- rpart(General_Grade ~., data = rand.train, method = "class")

#predict.ips <- predict(model.ips, ips.test2)
#predict.rand <- predict(model.rand, rand.test2)
#predict.a549 <- predict(model.a549, a549.test2)
#predict.randomdata.a549 <- predict(model.randomdata.a549, randomdata.a549.test2)
#predict.genes <- predict(model.genes.rf, genes.test2)
#fancyRpartPlot(model.ips.rpart)

#fit <- rpart(General_Stages ~ ., data = stage, method = "class")

# Get importance
importance.mcf7 <- varImp(mcf7.model.rf, scale = FALSE)
#importance.genes.rf <- varImp(model.genes.rf, scale = FALSE)

#importance.ips <- varImp(model.ips, scale=FALSE)
#importance.rand <- varImp(model.rand, scale = FALSE)
#importance.a549 <- varImp(model.a549, scale = FALSE)
#importance.randomdata.a549 <- varImp(model.randomdata.a549, scale = FALSE)
print(importance.mcf7)
plot(importance.mcf7)
#print(importance.randomdata.a549)
#plot(importance.randomdata.a549)
#print(importance.genes.rf)
#plot(importance.genes.rf)

importance.mcf7.2 <- importance(mcf7.model.rf2)
#importance.rand <- importance(model.rand)
#importance.tgenes <- importance(model.genes.rf)
#importance.a549.2 <- importance(model.a549)
#importance.randomdata.a549.2 <- importance(model.randomdata.a549)

# summarize and plot importance
print(importance.mcf7.2)
plot(importance.mcf7.2)
#print(importance.rand)
#plot(importance.rand)
#plot(importance.tgenes)

varImportance.mcf7 <- data.frame(Variables = row.names(importance.mcf7.2), 
                                  Importance = round(importance.mcf7.2[ ,'MeanDecreaseGini'],2))

# Create a rank variable based on importance
rankImportance.mcf7 <- varImportance.mcf7 %>%
  mutate(Rank = paste0('#',dense_rank(desc(Importance))))

# Use ggplot2 to visualize the relative importance of variables
ggplot(rankImportance.mcf7, aes(x = reorder(Variables, Importance), 
                                 y = Importance, fill = Importance)) +
  geom_bar(stat='identity') + 
  geom_text(aes(x = Variables, y = 0.5, label = Rank),
            hjust=0, vjust=0.55, size = 4, colour = 'darkblue') +
  labs(x = 'Variables') +
  theme(axis.text.y = element_text(size = 6)) + #, color = ifelse(rankImportance.mcf7$Variables %in% ipsgene, "red", "black"))) + 
  coord_flip() #+ 
#theme_few()







#Additional Filter Parameters =======================================
grades <- factor(c("1", "2", "3"))
grades <- ordered(grades, levels = c("1", "2", "3"))
stages <- factor(c("1", "2", "3", "4"))
stages <- ordered(stages, levels = c("1", "2", "3", "4"))

stage_diffs <- patientDF
#stage_diffs <- subset(patientDF, GENE %in% totaltestgenes)
stage_diffs <- subset(stage_diffs, GENE %in% mcf7ips_overlap) #Subset rows that are matches to those in the targets
#stage_diffs <- subset(stage_diffs, GENE %in% gene)
#stage_diffs <- subset(stage_diffs, GRADE %in% grades)
#stage_diffs <- subset(stage_diffs, TUMOR_STAGE %in% stages)

#stage_diffs <- filter(stage_diffs, TUMOR_STAGE == 4)
#stage_diffs <- filter(stage_diffs, ((GRADE == 1) | (GRADE == 2)))
stage_diffs <- filter(stage_diffs, GRADE == 3)
stage_diffs <- filter(stage_diffs, ER_STATUS == "+")
stage_diffs <- filter(stage_diffs, PR_STATUS == "+")
stage_diffs <- filter(stage_diffs, HER2_STATUS == "-")

#stage_diffs <- filter(stage_diffs, RATE_OF_PROLIF == "Low")

patient <- dcast(stage_diffs, PATIENT_ID+GRADE+TUMOR_STAGE+ER_STATUS+HER2_STATUS+PR_STATUS ~ GENE, value.var = "EXPRESSION_LEVEL")

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
patient <- subset(patient, select = -c(GRADE, TUMOR_STAGE, PR_STATUS, ER_STATUS, HER2_STATUS))
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
grade4 <- test1
#heatmap(test1, Colv=NA, col=greenred(10),scale="none")
#cor(t(test1))
t <- 1-cor(t(test1))
x <- cor(t(test1))
y <- cor(t(grade1), t(grade4))
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
  geom_text(aes(x=variable, y=Gene, label=round(value, digits = 2)), size=3) + 
  xlab("Correlation between Stage 1 vs Stage 4") + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=8),
        legend.text = element_text(size=8), 
        axis.title.x = element_text(size=8)) 

## For plot to sort numerically
extracted <- arrange(extracted, desc(value))
ggplot(extracted) +
  geom_tile(aes(x=variable, y=Gene, fill=value)) + 
  scale_fill_distiller(palette = "RdBu") + 
  scale_y_discrete(name="", limits = rev(extracted$Gene)) + 
  geom_text(aes(x=variable, y=Gene, label = round(value, digits = 2)), size=3) + 
  xlab("Correlation between Stage 1 vs Stage 4") + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=8),
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

results$editedpvalues <- NA
for (row in 1:nrow(results)) {
  if (results$pvalues[row] < 0.001) {
    results$editedpvalues[row] <- "<0.001"
  } else {
    results$editedpvalues[row] <- round(results$pvalues[row], digits = 4) 
  }
}

ggplot(results) + 
  geom_tile(aes(x=Stage, y=Genes, fill = colorscale), color = "white") + 
  scale_fill_brewer(palette = "RdBu") +
  scale_y_discrete(name="", limits = results$Genes) + 
  geom_text(aes(x=Stage, y=Genes, label = editedpvalues), size = 3) + 
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
#t_patient_means <- t_patient_m1

t_patient_means$mean_diff_1 <- t_patient_means$mean_1-t_patient_means$mean_1
t_patient_means$mean_diff_2 <- t_patient_means$mean_2-t_patient_means$mean_1
t_patient_means$mean_diff_3 <- t_patient_means$mean_3-t_patient_means$mean_1
t_patient_means$mean_diff_4 <- t_patient_means$mean_4-t_patient_means$mean_1

cast_patient <- cbind.data.frame(t_patient_means$Row.names, t_patient_means$mean_diff_1, t_patient_means$mean_diff_2, t_patient_means$mean_diff_3)#, t_patient_means$mean_diff_4)
colnames(cast_patient) <- c("Genes", "Stage1", "Stage2", "Stage3", "Stage4")
#colnames(cast_patient) <- c("Genes", "HER2+", "ER+/HER2-", "ER-/HER2-")
colnames(cast_patient) <- c("Genes", "Grade1", "Grade2", "Grade3")
melt_patient <- melt(cast_patient)
colnames(melt_patient) <- c("Genes", "Stages", "Difference")
melt_patient$Genes <- as.factor(melt_patient$Genes)
#plot_patient <- cbind.data.frame(t_patient$Row.names, t_patient$mean_diff)

ggplot(melt_patient, aes(x=Stages, y=Genes)) + 
  geom_tile(aes(fill=Difference)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="", limits = rev(levels(melt_patient$Genes))) + 
  geom_text(aes(x=Stages, y=Genes, label = round(Difference, digits = 2)), size=3) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8))


# To see descending values
#melt_patient <- filter(melt_patient, Stages == "ER+/HER2-")
melt_patient <- filter(melt_patient, Stages == "Grade3")
melt_patient <- arrange(melt_patient, desc(Difference))
#melt_patient <- arrange(melt_patient, rev(Difference))
ggplot(melt_patient) + 
  geom_tile(aes(x=Stages, y=Genes, fill = Difference)) + 
  scale_fill_distiller(palette = "RdBu") +
  scale_y_discrete(name="", limits = rev(melt_patient$Genes)) + 
  geom_text(aes(x=Stages, y=Genes, label = round(Difference, digits = 2)), size=3) + 
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
