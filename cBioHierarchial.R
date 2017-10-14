library(tidyverse)
library(reshape2)
library(stringr)

file.path <- ("~/Desktop/brca_metabric")
tumordata <- read.csv("data_clinical_supp_sample.csv", sep = ",", stringsAsFactors = FALSE)
patientdata <- read.csv("data_clinical_supp_patient.csv", sep = ",", stringsAsFactors = FALSE)
expression <- read.csv("data_expression.csv", sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

#Merging Data frames
patientDF <- merge(x = tumordata, y = patientdata, by = "PATIENT_ID", all = TRUE)

x <- subset(expression, select = -Entrez_Gene_Id)
x <- melt(x)

#m1 <- as.data.frame(t(expression))
#d2 <- data.frame(r1= row.names(m1), m1, row.names=NULL) 

y <- dcast(x, variable ~ Hugo_Symbol)


#z <- transpose(x)

write.csv(y, "test1.csv")

test <- read.csv("test1.csv", sep = ",", stringsAsFactors = FALSE)

#dcast(melt(as.matrix(expression)), Var2~paste0('r', Var1), value.var='value')


#Sorting patients based on grade

stage4 <- filter(patientDF, TUMOR_STAGE == 3 | TUMOR_STAGE == 4)

df <- filter(y, variable %in% stage4$PATIENT_ID)
