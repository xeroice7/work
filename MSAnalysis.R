### Packages 
install.packages("dplyr")
install.packages("ggplot2")
install.packages("stringr")
install.packages("scales")
install.packages("devtools")
install.packages("cowplot")
install.packages("data.table")
install.packages("tidyverse")
install.packages("VennDiagram")
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(reshape2)
library(devtools)
library(cowplot)
library(data.table)
library(tidyverse)
library(VennDiagram)


### Bring in the data

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
msdata <- data.table(read.csv(paste(file.path, "TotalSummaryFDROnly.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
names(msdata)[names(msdata) == 'A1111AAAAAAAA'] <- 'Line'  #Change the first column to "Line"

### Analyze "Total Evidence" portion
#Get rid of non-human columns

#human <- msdata[(msdata$Species=='HUMAN'),]
#a549 <- filter(human, grepl('A549', Line))


# FiPS/IMR comparison
fips <- data.frame(read.csv(paste(file.path, "IMRFiPSSummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(fips$X, " ", 4)
fips <- cbind.data.frame(splitcolumns, fips)
names(fips)[names(fips) == '1'] <- 'Line'  #Change the first column to "Line"
names(fips)[names(fips) == '2'] <- 'PassNum'  #Change the second column to "PassNum"
names(fips)[names(fips) == '3'] <- 'Biotin'  #Change the first column to "Line"
fips <- subset(fips, select=-c(`4`, X))
fips <- filter(fips, Species == "HUMAN") #Filtering only human hits
splitcolumns <- colsplit(fips$Gene.Names," ",c("Gene","Synonyms"))
fips <- cbind.data.frame(splitcolumns, fips)
fips <- subset(fips, select=-Gene.Names)
fipsslim <- data.frame(fips$Line, fips$Gene, fips$Synonyms, fips$Biotin, fips$Unused, fips$Total, fips$Gene.Ontology)
colnames(fipsslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

imrslim <- filter(fipsslim, Line == "IMR90")  #Keep only IMR
imrbiotin <- filter(imrslim, Biotin == "+") #Keep only +Biotin for IMR
colnames(imrbiotin)[colnames(imrbiotin) == 'Unused'] <- 'UnusedIMRBiotin'
colnames(imrbiotin)[colnames(imrbiotin) == 'Total'] <- 'TotalIMRBiotin'

imrnobiotin <- filter(imrslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(imrnobiotin)[colnames(imrnobiotin) == 'Unused'] <- 'UnusedIMRNoBiotin'
colnames(imrnobiotin)[colnames(imrnobiotin) == 'Total'] <- 'TotalIMRNoBiotin'

imr <- merge.data.frame(imrbiotin, imrnobiotin, by = "Gene", all = TRUE)
imr$UnusedIMRBiotin <- as.numeric(as.character(imr$UnusedIMRBiotin))
imr$UnusedIMRNoBiotin <- as.numeric(as.character(imr$UnusedIMRNoBiotin))
imr$TotalIMRBiotin <- as.numeric(as.character(imr$TotalIMRBiotin))
imr$TotalIMRNoBiotin <- as.numeric(as.character(imr$TotalIMRNoBiotin))

na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}
imr$UnusedIMRBiotin <- na.zero(imr$UnusedIMRBiotin)
imr$TotalIMRBiotin <- na.zero(imr$TotalIMRBiotin)
imr$UnusedIMRNoBiotin <- na.zero(imr$UnusedIMRNoBiotin)
imr$TotalIMRNoBiotin <- na.zero(imr$TotalIMRNoBiotin)

df1 <- data.frame(imrnobiotin$Gene, imrnobiotin$Line)
colnames(df1) <- c("Gene", "Line")
df2 <- data.frame(imrbiotin$Gene, imrbiotin$Line)
colnames(df2) <- c("Gene", "Line")

df3 <- setdiff(df1, df2)
no1 <- nrow(df1)
df4 <- setdiff(df2, df1)
no2 <- nrow(df2)
df5 <- semi_join(df1, df2)
no5 <- nrow(df5)

grid.newpage()
draw.pairwise.venn(no1, no2, no5, category = c("IMR90 No Biotin", "IMR + Biotin"), lty = rep("blank", 2), 
                  fill = c("blue", "red"), alpha = rep(0.5, 2), cat.pos = c(0, 0),
                  cat.dist = rep(0.025, 2), scaled = TRUE)

imrbiotinunique <- anti_join(df2, df1)
write.csv(imrbiotinunique, "UniqueIMRBiotin.csv")
#imr <- merge.data.frame(imrbiotin, imrnobiotin, by = "Gene")
#imr$synonyms <- coalesce(imr$Synonyms.x, imr$Synonyms.y)
#imr$ontology <- coalesce(imr$Ontology.x, imr$Ontology.y)
#imr$line <- coalesce(imr$Line.x, imr$Line.y)
#imr <- subset(imr, select = -c(Synonyms.x, Biotin.x, Synonyms.y, Biotin.y, `Gene Ontology.x`, `Gene Ontology.y`, Line.x, Line.y))
#imr$UnusedIMRBiotin <- as.numeric(as.character(imr$UnusedIMRBiotin))
#imr$UnusedIMRNoBiotin <- as.numeric(as.character(imr$UnusedIMRNoBiotin))
#imr$TotalIMRBiotin <- as.numeric(as.character(imr$TotalIMRBiotin))
#imr$TotalIMRNoBiotin <- as.numeric(as.character(imr$TotalIMRNoBiotin))

fipsslim <- filter(fipsslim, Line == "FiPS4F5")  #Keep only FiPS4F5
fipsbiotin <- filter(fipsslim, Biotin == "+") #Keep only +Biotin for FiPS4F5
colnames(fipsbiotin)[colnames(fipsbiotin) == 'Unused'] <- 'UnusedFiPSBiotin'
colnames(fipsbiotin)[colnames(fipsbiotin) == 'Total'] <- 'TotalFiPSBiotin'

fipsnobiotin <- filter(fipsslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(fipsnobiotin)[colnames(fipsnobiotin) == 'Unused'] <- 'UnusedFiPSNoBiotin'
colnames(fipsnobiotin)[colnames(fipsnobiotin) == 'Total'] <- 'TotalFiPSNoBiotin'

fips <- merge.data.frame(fipsbiotin, fipsnobiotin, by = "Gene", all = TRUE)
fips$UnusedFiPSBiotin <- as.numeric(as.character(fips$UnusedFiPSBiotin))
fips$UnusedFiPSNoBiotin <- as.numeric(as.character(fips$UnusedFiPSNoBiotin))
fips$TotalFiPSBiotin <- as.numeric(as.character(fips$TotalFiPSBiotin))
fips$TotalFiPSNoBiotin <- as.numeric(as.character(fips$TotalFiPSNoBiotin))

fips$UnusedFiPSBiotin <- na.zero(fips$UnusedFiPSBiotin)
fips$UnusedFiPSNoBiotin <- na.zero(fips$UnusedFiPSNoBiotin)
fips$TotalFiPSBiotin <- na.zero(fips$TotalFiPSBiotin)
fips$TotalFiPSNoBiotin <- na.zero(fips$TotalFiPSNoBiotin)

df6 <- data.frame(fipsnobiotin$Gene, fipsnobiotin$Line)
colnames(df6) <- c("Gene", "Line")
df7 <- data.frame(fipsbiotin$Gene, fipsbiotin$Line)
colnames(df7) <- c("Gene", "Line")

df8 <- setdiff(df6, df7)
no6 <- nrow(df6)
df9 <- setdiff(df7, df6)
no7 <- nrow(df7)
df10 <- semi_join(df6, df6)
no10 <- nrow(df10)

grid.newpage()
draw.pairwise.venn(no6, no7, no10, category = c("FiPS4F5 No Biotin", "FiPS4F5 + Biotin"), lty = rep("blank", 2), 
                   fill = c("blue", "red"), alpha = rep(0.5, 2), cat.pos = c(0, 0),
                   cat.dist = rep(0.025, 2), scaled = TRUE)



fipsbiotinunique <- anti_join(df7, df6)
write.csv(fipsbiotinunique, "UniqueFiPSBiotin.csv")



#fips$synonyms <- coalesce(fips$Synonyms.x, fips$Synonyms.y)
#fips$ontology <- coalesce(fips$Ontology.x, fips$Ontology.y)
#fips$line <- coalesce(fips$Line.x, fips$Line.y)
#fips <- subset(fips, select = -c(Synonyms.x, Biotin.x, Synonyms.y, Biotin.y, `Gene Ontology.x`, `Gene Ontology.y`, Line.x, Line.y))
#fips$UnusedFiPSBiotin <- as.numeric(as.character(fips$UnusedFiPSBiotin))
#fips$UnusedFiPSNoBiotin <- as.numeric(as.character(fips$UnusedFiPSNoBiotin))
#fips$TotalFiPSBiotin <- as.numeric(as.character(fips$TotalFiPSBiotin))
#fips$TotalFiPSNoBiotin <- as.numeric(as.character(fips$TotalFiPSNoBiotin))

imrtotal <- filter(imr, TotalIMRNoBiotin < TotalIMRBiotin)
fipstotal <- filter(fips, TotalFiPSNoBiotin < TotalFiPSBiotin)

df11 <- data.frame(fipstotal$Gene)
colnames(df11) <- c("Gene")
df12 <- data.frame(imrtotal$Gene)
colnames(df12) <- c("Gene")

df13 <- setdiff(df11, df12)
no11 <- nrow(df11)
df14 <- setdiff(df12, df11)
no12 <- nrow(df12)
df15 <- semi_join(df11, df12)
no15 <- nrow(df15)

grid.newpage()
draw.pairwise.venn(no11, no12, no15, category = c("FiPS4F5", "IMR90"), lty = rep("blank", 2), 
                   fill = c("blue", "red"), alpha = rep(0.5, 2), cat.pos = c(0, 0),
                   cat.dist = rep(0.025, 2), scaled = TRUE)

imrunique <- anti_join(df12, df11)
fipsunique <- anti_join(df11, df12)
imrfips <- intersect(df11, df12)
write.csv(imrunique, "UniqueIMR.csv")
write.csv(fipsunique, "UniqueFiPS.csv")
write.csv(imrfips, "FiPSIMROverlap.csv")
write.csv(imrtotal$Gene, "IMRReal.csv")
write.csv(fipstotal$Gene, "FiPSReal.csv")

#To compare which is bigger unused
imrvsfips <- merge.data.frame(fipstotal, imrtotal, by = "Gene", all = TRUE)

imrvsfips$UnusedFiPSBiotin <- na.zero(imrvsfips$UnusedFiPSBiotin)
imrvsfips$UnusedFiPSNoBiotin <- na.zero(imrvsfips$UnusedFiPSNoBiotin)
imrvsfips$TotalFiPSBiotin <- na.zero(imrvsfips$TotalFiPSBiotin)
imrvsfips$TotalFiPSNoBiotin <- na.zero(imrvsfips$TotalFiPSNoBiotin)
imrvsfips$UnusedIMRBiotin <- na.zero(imrvsfips$UnusedIMRBiotin)
imrvsfips$UnusedIMRNoBiotin <- na.zero(imrvsfips$UnusedIMRNoBiotin)
imrvsfips$TotalIMRBiotin <- na.zero(imrvsfips$TotalIMRBiotin)
imrvsfips$TotalIMRNoBiotin <- na.zero(imrvsfips$TotalIMRNoBiotin)

fipsup <- filter(imrvsfips, TotalFiPSBiotin > TotalIMRBiotin)
imrup <- filter(imrvsfips, TotalFiPSBiotin < TotalIMRBiotin)
write.csv(fipsup$Gene, "FiPSUp.csv")
write.csv(imrup$Gene, "IMRUp.csv")










#imrfips <- merge.data.frame(imrtotal, fipstotal, by = "Gene")

















fipstotal <- merge(fipsnobiotin, fipsbiotin, by = "Gene", all = TRUE) #merge multiple data sets based on gene names
fipstotal$synonyms <- coalesce(fipstotal$Synonyms.x, fipstotal$Synonyms.y)
fipstotal$ontology <- coalesce(fipstotal$Ontology.x, fipstotal$Ontology.y)
fipstotal$biotin <- coalesce(fipstotal$Biotin.x, fipstotal$Biotin.y)
fipstotal$line <- coalesce(fipstotal$Line.x, fipstotal$Line.y)
fipstotal$unused <- coalesce(fipstotal$Unused.x, fipstotal$Unused.y)
fipstotal$total <- coalesce(fipstotal$Total.x, fipstotal$Total.y)
fipstotal <- subset(fipstotal, select = c(Gene, line, biotin, unused, total, synonyms, ontology))





#A549 
a549 <- data.frame(read.csv(paste(file.path, "A549SummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(a549$X, "/", 2)
a549 <- cbind.data.frame(splitcolumns, a549)
names(a549)[names(a549) == '1'] <- 'Line'  #Change the first column to "Line"
names(a549)[names(a549) == '2'] <- 'YN'  #Change the second column to "YN"
splitcolumns <- str_split_fixed(a549$Line, " ", 2)
a549 <- cbind.data.frame(splitcolumns, a549)
names(a549)[names(a549) == '1'] <- 'Cell'  #Change the first column to "Cell"
names(a549)[names(a549) == '2'] <- 'Dox'  #Change the first column to "Line"
splitcolumns <- str_split_fixed(a549$YN, "", 2)
splitcolumns <- splitcolumns[,1]
a549 <- cbind.data.frame(splitcolumns, a549)
names(a549)[names(a549) == "splitcolumns"] <- "Biotin"
a549 <- subset(a549, select=-c(Line, YN, X))
a549 <- filter(a549, grepl('HUMAN', Species)) #Filtering only human hits


#a549dox <- filter(a549, Dox == "3d+Dox")

a549slim <- data.frame(a549$Biotin, a549$Dox, a549$Unused, a549$Total, a549$Gene.Names, a549$Gene.Ontology)
splitcolumns <- colsplit(a549slim$a549.Gene.Names," ",c("Gene","Synonyms"))
a549slim <- cbind.data.frame(splitcolumns, a549slim)
a549slim <- subset(a549slim, select=-a549.Gene.Names)

dox <- filter(a549slim, a549slim$a549.Dox == "3d+Dox") #Keep only 3d+Dox for now
doxnobiotin <- filter(dox, dox$a549.Biotin == "N") #subset N only rows
colnames(doxnobiotin) <- c("Gene", "Synonyms", "Biotin", "Dox", "DoxNoBiotinUnused", "DoxNoBiotinTotal", "Ontology")
doxyesbiotin <- filter(dox, dox$a549.Biotin == "+") #subset + only rows
colnames(doxyesbiotin) <- c("Gene", "Synonyms", "Biotin", "Dox", "DoxYesBiotinUnused", "DoxYesBiotinTotal", "Ontology")

nodox <- filter(a549slim, a549slim$a549.Dox == "No Dox") 
nodoxnobiotin <- filter(nodox, nodox$a549.Biotin == "N") #subset N only rows
colnames(nodoxnobiotin) <- c("Gene", "Synonyms", "Biotin", "Dox", "NoDoxNoBiotinUnused", "NoDoxNoBiotinTotal", "Ontology")
nodoxyesbiotin <- filter(nodox, nodox$a549.Biotin == "+") #subset + only rows
colnames(nodoxyesbiotin) <- c("Gene", "Synonyms", "Biotin", "Dox", "NoDoxYesBiotinUnused", "NoDoxYesBiotinTotal", "Ontology")

a549total <- merge(doxnobiotin, doxyesbiotin, by = "Gene", all = TRUE) #merge multiple data sets based on gene names
#a549total$dox <- coalesce(a549total$Dox.x, a549total$Dox.y)
#a549total$biotin <-coalesce(a549total$Biotin.x, a549total$Biotin.y)
a549total$synonyms <- coalesce(a549total$Synonyms.x, a549total$Synonyms.y)
a549total$ontology <- coalesce(a549total$Ontology.x, a549total$Ontology.y)
a549total <- subset(a549total, select = -c(Synonyms.x, Biotin.x, Dox.x, Synonyms.y, Biotin.y, Dox.y, Ontology.x, Ontology.y))

a549total <- merge(nodoxyesbiotin, a549total, by = "Gene", all = TRUE)
a549total$synonyms <- coalesce(a549total$synonyms, a549total$Synonyms)
a549total$ontology <- coalesce(a549total$ontology, a549total$Ontology)
a549total <- subset(a549total, select = -c(Synonyms, Biotin, Dox, Ontology))

a549total <- merge(nodoxnobiotin, a549total, by = "Gene", all = TRUE)
a549total$synonyms <- coalesce(a549total$synonyms, a549total$Synonyms)
a549total$ontology <- coalesce(a549total$ontology, a549total$Ontology)
a549total <- subset(a549total, select = -c(Synonyms, Biotin, Dox, Ontology))
#write.csv(a549total, "MyData1.csv")

#Plotting which are likely real targets

aa <- ggplot(a549total) + 
  geom_point(aes(x=NoDoxNoBiotinUnused, y=NoDoxYesBiotinUnused)) + 
  scale_y_continuous(trans = 'log', breaks = c(0, 0.01, 0.1, 1, 2, 3, 4, 5, 10, 15, 20, 30, 50, 100)) + 
  scale_x_continuous(trans = 'log', breaks = c(0, 0.01, 0.1, 1, 2, 3, 4, 5, 10, 15, 20, 30, 50, 100)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_text(aes(x=NoDoxNoBiotinUnused, y=NoDoxYesBiotinUnused, label=ifelse(NoDoxYesBiotinUnused>NoDoxNoBiotinUnused,as.character(Gene),'')),hjust=0, vjust=0) #label only genes above the 1 line

bb <- ggplot(a549total) + 
  geom_point(aes(x=DoxNoBiotinUnused, y=DoxYesBiotinUnused)) + 
  scale_y_continuous(trans = 'log', breaks = c(0, 0.01, 0.1, 1, 2, 3, 4, 5, 10, 15, 20, 30, 50, 100)) + 
  scale_x_continuous(trans = 'log', breaks = c(0, 0.01, 0.1, 1, 2, 3, 4, 5, 10, 15, 20, 30, 50, 100)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_text(aes(x=DoxNoBiotinUnused, y=DoxYesBiotinUnused, label=ifelse(DoxYesBiotinUnused>DoxNoBiotinUnused,as.character(Gene),'')),hjust=0, vjust=0) #label only genes above the 1 line

plot_grid(aa, bb, ncol=2)

#Gene Ontology
a549go <- a549total
splitcolumns <- str_split_fixed(a549go$ontology, ";", 100)
a549go <- cbind.data.frame(a549go, splitcolumns)
a549go <- subset(a549go, select = -c(ontology))
a549godox <- filter(a549go, DoxNoBiotinTotal< DoxYesBiotinTotal)
a549gonodox <- filter(a549go, NoDoxNoBiotinTotal < NoDoxYesBiotinTotal)


a549meltdox <- subset(a549godox, select = -c(NoDoxNoBiotinUnused, NoDoxNoBiotinTotal, NoDoxYesBiotinUnused, NoDoxYesBiotinTotal, DoxNoBiotinUnused,
                                             DoxNoBiotinTotal, DoxYesBiotinUnused, DoxYesBiotinTotal, synonyms, Gene))
a549meltdox$var <- "NA"
a549meltdox <- melt(a549meltdox, id.var='var')
a549meltdox$value <- as.factor(a549meltdox$value)
a549meltnodox <- subset(a549gonodox, select = -c(NoDoxNoBiotinUnused, NoDoxNoBiotinTotal, NoDoxYesBiotinUnused, NoDoxYesBiotinTotal, DoxNoBiotinUnused,
                                             DoxNoBiotinTotal, DoxYesBiotinUnused, DoxYesBiotinTotal, synonyms, Gene))
a549meltnodox$var <- "NA"
a549meltnodox <- melt(a549meltnodox, id.var='var')
a549meltnodox$value <- as.factor(a549meltnodox$value)
a549meltdox[1501:5900,] <- NA #Adjust the rows in dox to match no dox, so we can bind them together in data frame
a549melttotal <- cbind.data.frame(a549meltnodox$value, a549meltdox$value)
colnames(a549melttotal) <- c("NoDox", "Dox")

d <- count(a549melttotal, Dox)
n <- count(a549melttotal, NoDox)

x <- ggplot(d) + 
      geom_bar(aes(x="", y=n, fill=Dox), stat="identity")

y <- ggplot(n) + 
  geom_bar(aes(x="", y=n, fill=NoDox), stat="identity")

plot_grid(x, y, ncol=2)



write.csv(a549melttotal, "GO.csv")








a549melt <- melt(a549melt, id.var = c('a549.Biotin','a549.Dox', 'a549.Gene.Names'))
#write.csv(test, file = "MyData.csv")
a549melt <- subset(a549melt, select = -c(variable))
colnames(a549melt) <- c("Biotin", "Dox", "Gene", "GO")
a549melt$GO <- as.factor(a549melt$GO)






























## Ideas
test <- data.frame(a549$Biotin, a549$Dox, a549$Total, a549$Gene.Names) #keep only relevant columns
test <- filter(test, test$a549.Dox == "3d+Dox") #Keep only 3d+Dox for now
test <- subset(test, select=-a549.Dox) #Get rid of the Dox column after filter
testN <- filter(test, test$a549.Biotin == "N") #subset N only rows
testY <- filter(test, test$a549.Biotin == "+") #subset + only rows
total <- merge(testN, testY, by = "a549.Gene.Names", all = TRUE) #merge two data sets based on gene names
colnames(total) <- c("Gene1", "No", "NoTotal", "Yes", "PosTotal") #rename columns
total <- subset(total, select=-c(No, Yes)) #Get rid of Yes or No columns
total[is.na(total)] <- 0 #Change all NAs to 0
splitcolumns <- colsplit(total$Gene," ",c("Gene","Synonyms"))
total <- cbind.data.frame(splitcolumns, total)
total <- subset(total, select=-Gene1)

ggplot(total) + 
  geom_point(aes(x=NoTotal, y=PosTotal)) + 
  scale_y_continuous(trans = 'log', breaks = c(0, 0.1, 1, 2, 3, 4, 5, 10, 15, 20, 50, 100)) + 
  scale_x_continuous(trans = 'log', breaks = c(0, 0.1, 1, 2, 3, 4, 5, 10, 15, 20, 50, 100)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_text(aes(x=NoTotal, y=PosTotal, label=ifelse(PosTotal>NoTotal,as.character(Gene),'')),hjust=0, vjust=0) #label only genes above the 1 line



### Other
msdata_anxa <- msdata[msdata$Gene.Names == "ANXA",]
msdata_nuc <- msdata[msdata$Gene.Names == "NCL"]
