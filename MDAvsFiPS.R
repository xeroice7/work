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
install.packages("Vennerable", repos="http://R-Forge.R-project.org") 
install.packages("reshape")
install.packages("graph")
install.packages("RBGL")
install.packages("gtools")
install.packages("xtable")
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
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("graph")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("RBGL")

# Define Function for converting NAs to 0s when comparing value scores
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

#FiPS
fips <- data.frame(read.csv(paste(file.path, "IMRFiPSSummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(fips$X, " ", 4)
fips <- cbind.data.frame(splitcolumns, fips)
names(fips)[names(fips) == '1'] <- 'Line'  #Change the first column to "Line"
names(fips)[names(fips) == '2'] <- 'PassNum'  #Change the second column to "PassNum"
names(fips)[names(fips) == '3'] <- 'Biotin'  #Change the first column to "Line"
fips <- subset(fips, select=-c(`4`, X))
fips <- filter(fips, Species == "HUMAN") #Filtering only human hits
splitcolumns <- str_split_fixed(fips$Gene.Names," ", 2)
fips <- cbind.data.frame(splitcolumns, fips)
names(fips)[names(fips) == '1'] <- 'Gene'  #Change the first column to "Line"
names(fips)[names(fips) == '2'] <- 'Synonyms'  #Change the second column to "PassNum"
fips <- subset(fips, select=-Gene.Names)
fipsslim <- data.frame(fips$Line, fips$Gene, fips$Synonyms, fips$Biotin, fips$Unused, fips$Total, fips$Gene.Ontology)
colnames(fipsslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

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

fipsdf1 <- data.frame(fipsnobiotin$Gene, fipsnobiotin$Line)
colnames(fipsdf1) <- c("Gene", "Line")
fipsno1 <- nrow(fipsdf1)
fipsdf2 <- data.frame(fipsbiotin$Gene, fipsbiotin$Line)
colnames(fipsdf2) <- c("Gene", "Line")
fipsno2 <- nrow(fipsdf2)

fipsdf3 <- setdiff(fipsdf1, fipsdf2)
fipsdf4 <- setdiff(fipsdf2, fipsdf1)
fipsdf5 <- semi_join(fipsdf1, fipsdf2)
fipsno5 <- nrow(fipsdf5)

grid.newpage()
draw.pairwise.venn(fipsno1, fipsno2, fipsno5, category = c("FiPS4F5 No Biotin", "FiPS4F5 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("purple", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

fipsbiotinunique <- anti_join(fipsdf2, fipsdf1)
write.csv(fipsbiotinunique, "UniqueFiPSBiotin.csv")


#MDA
mda <- data.frame(read.csv(paste(file.path, "MDASummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(mda$X, "/", 2)
mda <- cbind.data.frame(splitcolumns, mda)
names(mda)[names(mda) == '1'] <- 'Line'  #Change the first column to "Line"
names(mda)[names(mda) == '2'] <- 'YN'  #Change the second column to "YN"
splitcolumns <- str_split_fixed(mda$Line, " ", 3)
mda <- cbind.data.frame(splitcolumns, mda)
names(mda)[names(mda) == '1'] <- 'Cell'  #Change the first column to "Cell"
names(mda)[names(mda) == '3'] <- 'Dox'  #Change the first column to "Line"
splitcolumns <- str_split_fixed(mda$YN, " ", 2)
splitcolumns <- splitcolumns[,1]
mda <- cbind.data.frame(splitcolumns, mda)
names(mda)[names(mda) == "splitcolumns"] <- "Biotin"
mda <- subset(mda, select=-c(`2`, Line, YN, X))
mda <- filter(mda, Species == "HUMAN") #Filtering only human hits
mda <- filter(mda, Dox == "No Dox") #Filtering only No Dox hits
splitcolumns <- str_split_fixed(mda$Gene.Names," ", 2)
mda <- cbind.data.frame(splitcolumns, mda)
names(mda)[names(mda) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(mda)[names(mda) == '2'] <- 'Synonyms'  #Change the first column to "Line"
mda <- subset(mda, select=-Gene.Names)
mdaslim <- data.frame(mda$Cell, mda$Gene, mda$Synonyms, mda$Biotin, mda$Unused, mda$Total, mda$Gene.Ontology)
colnames(mdaslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

mdabiotin <- filter(mdaslim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(mdabiotin)[colnames(mdabiotin) == 'Unused'] <- 'UnusedMDABiotin'
colnames(mdabiotin)[colnames(mdabiotin) == 'Total'] <- 'TotalMDABiotin'

mdanobiotin <- filter(mdaslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(mdanobiotin)[colnames(mdanobiotin) == 'Unused'] <- 'UnusedMDANoBiotin'
colnames(mdanobiotin)[colnames(mdanobiotin) == 'Total'] <- 'TotalMDANoBiotin'

mda <- merge.data.frame(mdabiotin, mdanobiotin, by = "Gene", all = TRUE)
mda$UnusedMDABiotin <- as.numeric(as.character(mda$UnusedMDABiotin))
mda$UnusedMDANoBiotin <- as.numeric(as.character(mda$UnusedMDANoBiotin))
mda$TotalMDABiotin <- as.numeric(as.character(mda$TotalMDABiotin))
mda$TotalMDANoBiotin <- as.numeric(as.character(mda$TotalMDANoBiotin))

mda$UnusedMDABiotin <- na.zero(mda$UnusedMDABiotin)
mda$UnusedMDANoBiotin <- na.zero(mda$UnusedMDANoBiotin)
mda$TotalMDABiotin <- na.zero(mda$TotalMDABiotin)
mda$TotalMDANoBiotin <- na.zero(mda$TotalMDANoBiotin)

mdadf1 <- data.frame(mdanobiotin$Gene, mdanobiotin$Line)
colnames(mdadf1) <- c("Gene", "Line")
mdano1 <- nrow(mdadf1)
mdadf2 <- data.frame(mdabiotin$Gene, mdabiotin$Line)
colnames(mdadf2) <- c("Gene", "Line")
mdano2 <- nrow(mdadf2)

mdadf3 <- setdiff(mdadf1, mdadf2)
mdadf4 <- setdiff(mdadf2, mdadf2)
mdadf5 <- semi_join(mdadf1, mdadf2)
mdano5 <- nrow(mdadf5)

grid.newpage()
draw.pairwise.venn(mdano1, mdano2, mdano5, category = c("MDA No Biotin", "MDA + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("light green", "magenta"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mdabiotinunique <- anti_join(mdadf2, mdadf1)
write.csv(mdabiotinunique, "UniqueMDABiotin.csv")


# To compare together
mdatotal <- filter(mda, TotalMDANoBiotin < TotalMDABiotin)
fipstotal <- filter(fips, TotalFiPSNoBiotin < TotalFiPSBiotin)

fipsvsmda <- merge.data.frame(mdatotal, fipstotal, by = "Gene", all = TRUE)

fipsvsmda$UnusedMDABiotin <- na.zero(fipsvsmda$UnusedMDABiotin)
fipsvsmda$UnusedMDANoBiotin <- na.zero(fipsvsmda$UnusedMDANoBiotin)
fipsvsmda$TotalMDABiotin <- na.zero(fipsvsmda$TotalMDABiotin)
fipsvsmda$TotalMDANoBiotin <- na.zero(fipsvsmda$TotalMDANoBiotin)
fipsvsmda$UnusedFiPSBiotin <- na.zero(fipsvsmda$UnusedFiPSBiotin)
fipsvsmda$UnusedFiPSNoBiotin <- na.zero(fipsvsmda$UnusedFiPSNoBiotin)
fipsvsmda$TotalFiPSBiotin <- na.zero(fipsvsmda$TotalFiPSBiotin)
fipsvsmda$TotalFiPSNoBiotin <- na.zero(fipsvsmda$TotalFiPSNoBiotin)

mdadf6 <- data.frame(unique(mdatotal$Gene))
colnames(mdadf6) <- c("Gene")
mdano6 <- nrow(mdadf6)
fipsdf6 <- data.frame(unique(fipstotal$Gene))
colnames(fipsdf6) <- c("Gene")
fipsno6 <- nrow(fipsdf6)

fipsmdadf1 <- setdiff(mdadf6, fipsdf6)
fipsmdadf2 <- setdiff(fipsdf6, mdadf6)
fipsmdadf3 <- semi_join(mdadf6, fipsdf6)
fipsmdano3 <- nrow(unique(fipsmdadf3))

grid.newpage()
draw.pairwise.venn(mdano6, fipsno6, fipsmdano3, category = c("MDA-MB-231", "FiPS4F5"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("light green", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mdaunique <- anti_join(mdadf6, fipsdf6)
fipsunique <- anti_join(fipsdf6, mdadf6)
mdafips <- intersect(mdadf6$Gene, fipsdf6$Gene)
write.csv(mdaunique, "UniqueMDAtoFIPS.csv")
write.csv(fipsunique, "UniqueFiPStoMDA.csv")
write.csv(mdafips, "FiPSMDAOverlap.csv")
write.csv(mdatotal$Gene, "MDAReal.csv")
write.csv(fipstotal$Gene, "FiPSReal.csv")
