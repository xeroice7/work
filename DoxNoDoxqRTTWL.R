###### R Script to analyze Dox/NoDox qRT data with standard 18S controls

### Packages needed

install.packages("dplyr")
install.packages("ggplot2")
install.packages("stringr")
library(dplyr)
library(ggplot2)
library(stringr)

### Bring in the data

file.path <- ("~/Desktop/Clay/qRT/")
#data <- data.frame(read.csv(paste(file.path, "Plate1090314_data.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
data <- data.frame(read.csv(paste(file.path, "41316twlMYH9DCDGRP78_data copy.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

### Deleting unnecessary rows and columnns from the data

data <- data[-c(1:6),]
data <- data[,-c(11:ncol(data))]
data <- data[,-c(4:9)]
data <- data[-1,]
colnames(data) <- c("Well", "Sample", "Target", "CT")
data$CT <- as.numeric(data$CT)
data <- data[!is.na(data$CT),]
data1 <- data[((data$Sample == 'MDA 231 ND')| (data$Sample == 'MDA 231 Dox')),]
splitcolumns <- str_split_fixed(data1$Sample, " ", 3)
data1 <- cbind.data.frame(data1, splitcolumns)
colnames(data1) <- c("Well", "Sample", "Target", "CT", "Line", "Dox", "Txt")
data1 <- subset(data1, select = -Dox)
data2 <- subset(data, Sample != 'MDA 231 ND')
data2 <- subset(data2, Sample != 'MDA 231 Dox')
splitcolumns <- str_split_fixed(data2$Sample, " ", 2)
data2 <- cbind.data.frame(data2, splitcolumns)
colnames(data2) <- c("Well", "Sample", "Target", "CT", "Line", "Txt")
data <- rbind.data.frame(data1, data2)
data <- subset(data, select = -c(Well, Sample))

### Grouping and making the calculations

data <- data %>%
          group_by(Target, Line, Txt) %>%
          mutate(meanCT = mean(CT))

controlmeanCT <- data[data$Target == "18s",]
controlmeanCT <- controlmeanCT[(!duplicated(controlmeanCT$Txt)),]

doxcontrolmeanCT <- controlmeanCT %>%
                      select(Target, Line, Txt, meanCT) %>%
                      filter(Txt == "Dox")
                      
                      
nodoxcontrolmeanCT <- controlmeanCT %>%
                        select(Target, Txt, meanCT) %>%
                        filter(Txt == "No")
                                
deltaCT <- c()                        
for (i in 1:nrow(data)) {
  
  if (data$Txt[i] == "No") {
    
    deltaCT[i] <- data$CT[i]-nodoxcontrolmeanCT$meanCT
  
  }else if (data$Txt[i] == "Dox") {
      
    deltaCT[i] <- data$CT[i]-doxcontrolmeanCT$meanCT
}
}

data <- cbind.data.frame(data, deltaCT)

data <- data %>%
          group_by(Target, Txt) %>%
          mutate(meandeltaCT = mean(deltaCT))

deltadeltaCT <- c()
for (i in 1:nrow(data)) {
  if (data$Txt[i] == "No") {
    deltadeltaCT[i] <- data$deltaCT[i]-data$meandeltaCT[i]
  } else if (data$Txt[i] == "Dox") {
    deltadeltaCT[i] <- data$deltaCT[i]-data$meandeltaCT[i-3]  #i-length(group_by(data, Target, Txt))
  }
    
}

data <- cbind.data.frame(data, deltadeltaCT)
data <- data %>%
          mutate(normalized = 2^(-deltadeltaCT))

data <- data %>%
          group_by(Target, Txt) %>%
          mutate(meannormalized = mean(normalized))

data <- data %>%
          group_by(Target, Txt) %>%
          mutate(stdev = sd(normalized))

data <- data %>%
          group_by(Target, Txt) %>%
          mutate(sem = sd(normalized)/sqrt(length(normalized)))


levels(data$Txt)[levels(data$Txt)=="No"] <- "No Dox" 
data$Txt <- ordered(data$Txt, levels = c("No Dox", "Dox"))

data1 <- data[!duplicated(data[,c('Target', 'Txt')]),]

ggplot(data1, aes(x=Target, y=meannormalized)) + 
  geom_bar(aes(fill=Txt), stat="identity", position = position_dodge(width = 0.90)) + 
  geom_hline(aes(yintercept=1), linetype="dashed", color="green") + 
  geom_errorbar(aes(ymin=meannormalized-sem, ymax=meannormalized+sem, fill=Txt), color="black", width = 0.2, position = position_dodge(width = 0.90)) +  
  scale_fill_brewer(palette = "Set1") + 
  ylab("Normalized Gene Expression") + 
  theme_gray() + 
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())
