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
data <- data.frame(read.csv(paste(file.path, "Plate1HCCJ2ZsGRP78CSCs050517_data.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

### Deleting unnecessary rows and columnns from the data

data <- data[-c(1:6),]
data <- data[,-c(11:ncol(data))]
data <- data[,-c(4:9)]
data <- data[-1,]
colnames(data) <- c("Well", "Sample", "Target", "CT")
data$CT <- as.numeric(data$CT)
data <- data[!is.na(data$CT),]
splitcolumns <- str_split_fixed(data$Sample, " ", 2)
data <- cbind.data.frame(data, splitcolumns)
data <- subset(data, select = -c(Well, Sample, `1`))
colnames(data) <- c("Target", "CT", "Population")
data <- filter(data, ((Population == "sGRP78-") | (Population == "sGRP78+")))
  
### Grouping and making the calculations

data <- data %>%
  group_by(Target, Population) %>%
  mutate(meanCT = mean(CT, na.rm = TRUE))

negdata <- filter(data, Population == "sGRP78-")
posdata <- filter(data, Population == "sGRP78+")

negcontrolmeanCT <- negdata[negdata$Target == "18S",]
negcontrolmeanCT <- negcontrolmeanCT$meanCT[1]
negdeltaCT <- negdata$CT-negcontrolmeanCT
negdeltaCT1 <- negdeltaCT[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35)]
negdeltaCT2 <- negdeltaCT[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)]
names <- unique(negdata$Target)
negmeanCT <- negdata$meanCT[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35)]
negDF <- cbind.data.frame(names, negdeltaCT1, negdeltaCT2, negmeanCT)
write.csv(names, "names.csv")

controlmeanCT <- controlmeanCT[(!duplicated(subset(controlmeanCT, select = c(Sample, Txt)))),]




###
file.path <- ("~/Desktop/Clay/qRT/")
#markerDF <- data.frame(read.csv(paste(file.path, "MCF7DoxiPSEMTMET.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
markerDF <- data.frame(read.csv(paste(file.path, "MCF7DoxiPSEMTMETDoxToNoDox.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
#colnames(markerDF) <- c("Gene", "Marker Type", "3d+Dox", "HUViPS4F1")
colnames(markerDF) <- c("Gene", "Marker Type", "3d+Dox")
markerDF <- melt(markerDF)
colnames(markerDF) <- c("Gene", "Marker Type", "Population", "Relative Gene Expression")

epi <- filter(markerDF, `Marker Type` == "Epi")
mes <- filter(markerDF, `Marker Type` == "Mes")
ran <- filter(markerDF, `Marker Type` == "Ran")
ms <- filter(markerDF, `Marker Type` == "MS")

zz <- ggplot(epi) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  #scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkblue", name = "") + 
  scale_fill_continuous(limits=c(0, 3), breaks=seq(0,3,by=1),low='white',high='darkblue') + 
  ggtitle("Epithelial Markers") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) 

yy <- ggplot(mes) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  #scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkblue", name = "") + 
  scale_fill_continuous(limits=c(0, 3), breaks=seq(0,3,by=1),low='white',high='darkblue') + 
  ggtitle("Mesenchymal Markers") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) 

xx <- ggplot(ran) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkblue", name = "") + 
  ggtitle("Random Markers") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

ww <- ggplot(ms) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkblue", name = "") + 
  ggtitle("MS Markers") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

plot_grid(zz,yy,xx,ww)

plot_grid(zz, yy)






























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
