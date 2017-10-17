#### Total Analysis of EMT to MCF7 sGRP78+

# Install the proper packages
install.packages("dplyr")
install.packages("xlsx")
install.packages("ggplot2")
install.packages("plotrix")

# Load the packages
library(ggplot2)
library(xlsx)
library(dplyr)
library(plotrix)

####### Define functions
#Clean the data
tidydata <- function(x) {
  x <- x[,c("X96fast", "X", "X.7")]
  x <- x[-c(1:7),]
  splitcolumns <- str_split_fixed(x$X96fast, " ", 2)
  x <- cbind.data.frame(x, splitcolumns)
  x <- subset(x, select = -c(X96fast, `1`))
  colnames(x) <- c("Target", "CT", "Population")
  x$CT <- as.numeric(x$CT)
  x <- x[!is.na(x$CT),]
}

#Isolate 18s mean
controlmean <- function(x) {
  x <- x %>%
    group_by(Target, Population) %>%
    mutate(meanCT = mean(CT, na.rm = TRUE))
}
  
############ First Data Set
setwd(dir = "Desktop")
setwd(dir = "Clay")
setwd(dir = "qRT")
datadf <- data.frame(read.csv("Plate1071017MCF7sGRP78EMTMShits_data.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

#Tidying the data and calculate 18s mean
datadf <- tidydata(datadf)
datadf <- controlmean(datadf)

#Filter out the different populations
csc <- as.data.frame(filter(datadf, Population == "sGRP78+"))
noncsc <- as.data.frame(filter(datadf, Population == "sGRP78-"))

csccontrolmeanCT <- csc[csc$Target == "18S",]
noncsccontrolmeanCT <- noncsc[noncsc$Target == "18S",]

csccontrolmeanCT <- csccontrolmeanCT$meanCT[1]
noncsccontrolmeanCT <- noncsccontrolmeanCT$meanCT[1]

#Delta CT
csc <- csc %>%
  group_by(Target) %>%
  mutate(deltaCT = CT-csccontrolmeanCT)

noncsc <- noncsc %>%
  group_by(Target) %>%
  mutate(deltaCT = CT-noncsccontrolmeanCT)

#Average Delta CT
csc <- csc %>%
  group_by(Target) %>%
  mutate(avgdeltaCT = mean(deltaCT, na.rm = TRUE))

noncsc <- noncsc %>%
  group_by(Target) %>%
  mutate(avgdeltaCT = mean(deltaCT, na.rm = TRUE))

#Delta Delta CT
noncsc$Target <- as.factor(noncsc$Target)
uniquenoncsc <- filter(noncsc, !duplicated(Target))

noncsc$deltadeltaCT <- noncsc$deltaCT-uniquenoncsc$avgdeltaCT[match(noncsc$Target, uniquenoncsc$Target)]
csc$deltadeltaCT <- csc$deltaCT-uniquenoncsc$avgdeltaCT[match(csc$Target, uniquenoncsc$Target)]

#Combine two populations into 1 data frame
cscDF <- rbind.data.frame(csc, noncsc)

cscDF <- cscDF %>%
  group_by(Population, Target) %>%
  mutate(normalized = 2^(-deltadeltaCT))

############ Second Data Set
datadf1 <- data.frame(read.csv("Plate1HCCJ2ZsGRP78CSCs050517_data.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE))

#Tidying the data and calculate 18s mean
datadf1 <- tidydata(datadf1)
datadf1 <- controlmean(datadf1)

#Filter out the different populations
csc <- as.data.frame(filter(datadf1, Population == "sGRP78+"))
noncsc <- as.data.frame(filter(datadf1, Population == "sGRP78-"))

csccontrolmeanCT <- csc[csc$Target == "18S",]
noncsccontrolmeanCT <- noncsc[noncsc$Target == "18S",]

csccontrolmeanCT <- csccontrolmeanCT$meanCT[1]
noncsccontrolmeanCT <- noncsccontrolmeanCT$meanCT[1]

#Delta CT
csc <- csc %>%
  group_by(Target) %>%
  mutate(deltaCT = CT-csccontrolmeanCT)

noncsc <- noncsc %>%
  group_by(Target) %>%
  mutate(deltaCT = CT-noncsccontrolmeanCT)

#Average Delta CT
csc <- csc %>%
  group_by(Target) %>%
  mutate(avgdeltaCT = mean(deltaCT, na.rm = TRUE))

noncsc <- noncsc %>%
  group_by(Target) %>%
  mutate(avgdeltaCT = mean(deltaCT, na.rm = TRUE))

#Delta Delta CT
noncsc$Target <- as.factor(noncsc$Target)
uniquenoncsc <- filter(noncsc, !duplicated(Target))

noncsc$deltadeltaCT <- noncsc$deltaCT-uniquenoncsc$avgdeltaCT[match(noncsc$Target, uniquenoncsc$Target)]
csc$deltadeltaCT <- csc$deltaCT-uniquenoncsc$avgdeltaCT[match(csc$Target, uniquenoncsc$Target)]

#Combine two populations into 1 data frame
cscDF1 <- rbind.data.frame(csc, noncsc)

cscDF1 <- cscDF1 %>%
  group_by(Population, Target) %>%
  mutate(normalized = 2^(-deltadeltaCT))

############ Third Data Set (contains Keratin 5 from first RNA run, had to be done on second plate because primers came in later)
datadf2 <- data.frame(read.csv("Plate1HCCJ2ZMCF7TotalEMT051017_data.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE))

#Tidying the data and calculate 18s mean
datadf2 <- tidydata(datadf2)
datadf2 <- controlmean(datadf2)

#Filter out the different populations
csc <- as.data.frame(filter(datadf2, Population == "sGRP78+"))
noncsc <- as.data.frame(filter(datadf2, Population == "sGRP78-"))

csccontrolmeanCT <- csc[csc$Target == "18S",]
noncsccontrolmeanCT <- noncsc[noncsc$Target == "18S",]

csccontrolmeanCT <- csccontrolmeanCT$meanCT[1]
noncsccontrolmeanCT <- noncsccontrolmeanCT$meanCT[1]

#Delta CT
csc <- csc %>%
  group_by(Target) %>%
  mutate(deltaCT = CT-csccontrolmeanCT)

noncsc <- noncsc %>%
  group_by(Target) %>%
  mutate(deltaCT = CT-noncsccontrolmeanCT)

#Average Delta CT
csc <- csc %>%
  group_by(Target) %>%
  mutate(avgdeltaCT = mean(deltaCT, na.rm = TRUE))

noncsc <- noncsc %>%
  group_by(Target) %>%
  mutate(avgdeltaCT = mean(deltaCT, na.rm = TRUE))

#Delta Delta CT
noncsc$Target <- as.factor(noncsc$Target)
uniquenoncsc <- filter(noncsc, !duplicated(Target))

noncsc$deltadeltaCT <- noncsc$deltaCT-uniquenoncsc$avgdeltaCT[match(noncsc$Target, uniquenoncsc$Target)]
csc$deltadeltaCT <- csc$deltaCT-uniquenoncsc$avgdeltaCT[match(csc$Target, uniquenoncsc$Target)]

#Combine two populations into 1 data frame
cscDF2 <- rbind.data.frame(csc, noncsc)

cscDF2 <- cscDF2 %>%
  group_by(Population, Target) %>%
  mutate(normalized = 2^(-deltadeltaCT))

#Combine two data sets into data frame
emt <- rbind.data.frame(cscDF, cscDF1, cscDF2)

emt <- emt %>%
  group_by(Population, Target) %>%
  mutate(meannormalized = mean(normalized, na.rm=TRUE))

emt <- emt %>%
  group_by(Population, Target) %>%
  mutate(SDnormalized = sd(normalized, na.rm=TRUE))

emt <- emt %>%
  group_by(Population, Target) %>%
  mutate(SEnormalized = SDnormalized/sqrt(4))

epi <- c("Keratin8", "Keratin5", "Areg", "EpCAM", "E-cadherin", "DSP")
mes <- c("Vimentin", "Twist1", "Snail2", "Snail1", "N-cadherin", "Fibronectin")

emt1 <- filter(emt, Target %in% epi)
emt1$Marker <- c("E")

emt2 <- filter(emt, Target %in% mes)
emt2$Marker <- c("M")

emt <- rbind.data.frame(emt1, emt2)

write.csv(emt, "TotalEMTAnalysis-sGRP78.csv")

#Setting up data to plot
plotemt <- emt[!duplicated(emt[,c("Target", "Population")]),]
emt_names <- c(
  `E` = "Epithelial \nMarkers",
  `M` = "Mesenchymal \nMarkers")

# Plotting for bar graph
ggplot(plotemt) + 
  facet_wrap(~ Target, scales = "free") + 
  geom_bar(aes(x = Population, y = meannormalized, fill = Population), width = 0.5, stat = "identity", position = "dodge") + 
  geom_errorbar(aes(x = Population, ymin = meannormalized - SEnormalized, ymax = meannormalized + SEnormalized, fill = Population), width = .08, position = position_dodge(0.5)) +
  theme(plot.title = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title.y = element_text(size = 7.5))


# Plotting for bar graph by Facets
ggplot(plotemt) + 
  facet_wrap(~ Marker, scales = "free") + 
  geom_bar(aes(x = Target, y = meannormalized, fill = Population), width = 0.5, stat = "identity", position = "dodge") + 
  geom_errorbar(aes(x = Target, ymin = meannormalized - SEnormalized, ymax = meannormalized + SEnormalized, fill = Population), width = .08, position = position_dodge(0.5)) +
  theme(plot.title = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title.y = element_text(size = 7.5))

# Plotting for Heat Map
ggplot(plotemt) + 
  facet_wrap(~ Marker, scales = "free", labeller = as_labeller(emt_names)) + 
  geom_tile(aes(x = Population, y = Target, fill = meannormalized)) +  
  scale_fill_gradient(low = "white", high = "darkgreen", name = "") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 9), 
        strip.text.x = element_text(size = 14, face = "bold"))








###
file.path <- ("~/Desktop/Clay/qRT/")
#markerDF <- data.frame(read.csv(paste(file.path, "MCF7sGRP78EMTMET.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
markerDF <- data.frame(read.csv(paste(file.path, "MCF7sGRP78EMTMETPosToNeg.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

names(markerDF)[names(markerDF) == 'sGRP78Pos'] <- 'sGRP78+'
markerDF <- melt(markerDF)
colnames(markerDF) <- c("Gene", "Marker Type", "Population", "Relative Gene Expression")

epi <- filter(markerDF, `Marker Type` == "E")
mes <- filter(markerDF, `Marker Type` == "M")
ran <- filter(markerDF, `Marker Type` == "R")

zz <- ggplot(epi) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_continuous(limits=c(0, 10), breaks=seq(0,10,by=2),low='white',high='darkgreen') + 
  #scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10), low = "white", high = "darkgreen", name = "") + 
  ggtitle("Epithelial Markers") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 10)) 

yy <- ggplot(mes) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_continuous(limits=c(0, 10), breaks=seq(0,10,by=2),low='white',high='darkgreen') + 
  #scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10), low = "white", high = "darkgreen", name = "") + 
  ggtitle("Mesenchymal Markers") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) 

xx <- ggplot(ran) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_continuous(limits=c(0, 10), breaks=seq(0,10,by=2),low='white',high='darkgreen') + 
  #scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10), low = "white", high = "darkgreen", name = "") + 
  ggtitle("Random Markers") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))


plot_grid(zz,yy,xx)
plot_grid(zz,yy)




marker <- c()
for (i in 1:nrow(emt)) {
  if(any(emt$Target[i]==epi)) {
    marker[i] == "E"
  } else if (any(emt$Target[i]==mes)) {
    marker[i] == "M"
  } else {
    marker[i] == "Other"
  }
  print(i)
}