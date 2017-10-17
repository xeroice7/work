install.packages("dplyr")
install.packages("ggplot2")
install.packages("stringr")
install.packages("scales")
install.packages("devtools")
install.packages("cowplot")
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(reshape2)
library(devtools)
library(cowplot)

### Bring in the data

file.path <- ("~/Desktop/Clay/qRT/")
data <- data.frame(read.csv(paste(file.path, "Plate1021417NCL_data.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

### Deleting unnecessary rows and columnns from the data

data <- data[-c(1:6),]
data <- data[,-c(11:ncol(data))]
data <- data[,-c(4:9)]
data <- data[-1,]
colnames(data) <- c("Well", "Sample", "Target", "CT")
data$CT <- as.numeric(data$CT)
#data <- data[!is.na(data$CT),]
data <- data[-c(97:nrow(data)),]
splitcolumns <- str_split_fixed(data$Sample, " ", 4)
data <- cbind.data.frame(data, splitcolumns)
data <- subset(data, select = -c(Well, Sample, `2`, `4`))
colnames(data) <- c("Target", "CT", "Sample", "Txt")

### Grouping and making the calculations

data <- data %>%
  group_by(Target, Sample, Txt) %>%
  mutate(meanCT = mean(CT, na.rm = TRUE))

controlmeanCT <- data[data$Target == "18S",]
controlmeanCT <- controlmeanCT[(!duplicated(subset(controlmeanCT, select = c(Sample, Txt)))),]

deltaCT <- data$CT - controlmeanCT$meanCT[((match(data$Sample, controlmeanCT$Sample)) && (match(data$Txt, controlmeanCT$Txt)))]
# or just subtract elementwise if this doesnt work?

data <- cbind.data.frame(data, deltaCT)

data <- data %>%
  group_by(Sample, Target, Txt) %>%
  mutate(meandeltaCT = mean(deltaCT, na.rm = TRUE))





deltadeltaCT <- c()
    for (i in 1:nrow(data)) {
      data %>%
        group_by(Target, Sample) %>%
          if (data$Txt[i] == "No" || data$Txt[i] == "Non-CSCs" || data$Sample[i] == "IMR90") {
            deltadeltaCT[i] <- data$deltaCT[i]-data$meandeltaCT[i]
          } else if ((data$Txt[i] == "3d") || (data$Txt[i] == "CSCs") || (data$Sample[i] == "FiPS4F5")) {
            deltadeltaCT[i] <- data$deltaCT[i]-data$meandeltaCT[i-1]
          } else {
            deltadeltaCT[i] <- "NA"
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