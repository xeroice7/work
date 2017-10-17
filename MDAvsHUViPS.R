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

#HUViPS
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
huv <- data.frame(read.csv(paste(file.path, "HUVECHUViPSSummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(huv$X, " ", 5)
huv <- cbind.data.frame(splitcolumns, huv)
names(huv)[names(huv) == '1'] <- 'Line'  #Change the first column to "Line"
names(huv)[names(huv) == '3'] <- 'PassNum'  #Change the second column to "PassNum"
names(huv)[names(huv) == '4'] <- 'Biotin'  #Change the first column to "Line"
huv <- subset(huv, select=-c(`2`, `5`, X))
huv <- filter(huv, Species == "HUMAN") #Filtering only human hits
splitcolumns <- str_split_fixed(huv$Gene.Names," ", 2)
huv <- cbind.data.frame(splitcolumns, huv)
names(huv)[names(huv) == '1'] <- 'Gene'  #Change the first column to "Line"
names(huv)[names(huv) == '2'] <- 'Synonyms'  #Change the second column to "PassNum"
huv <- subset(huv, select=-Gene.Names)
huvslim <- data.frame(huv$Line, huv$Gene, huv$Synonyms, huv$Biotin, huv$Unused, huv$Total, huv$Gene.Ontology)
colnames(huvslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

huvipsslim <- filter(huvslim, Line == "HUViPS4F1")  #Keep only HUViPS4F1
huvipsbiotin <- filter(huvipsslim, Biotin == "+") #Keep only +Biotin for HUViPS4F1
colnames(huvipsbiotin)[colnames(huvipsbiotin) == 'Unused'] <- 'UnusedHUViPSBiotin'
colnames(huvipsbiotin)[colnames(huvipsbiotin) == 'Total'] <- 'TotalHUViPSBiotin'

huvipsnobiotin <- filter(huvipsslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(huvipsnobiotin)[colnames(huvipsnobiotin) == 'Unused'] <- 'UnusedHUViPSNoBiotin'
colnames(huvipsnobiotin)[colnames(huvipsnobiotin) == 'Total'] <- 'TotalHUViPSNoBiotin'

huvips <- merge.data.frame(huvipsbiotin, huvipsnobiotin, by = "Gene", all = TRUE)
huvips$UnusedHUViPSBiotin <- as.numeric(as.character(huvips$UnusedHUViPSBiotin))
huvips$UnusedHUViPSNoBiotin <- as.numeric(as.character(huvips$UnusedHUViPSNoBiotin))
huvips$TotalHUViPSBiotin <- as.numeric(as.character(huvips$TotalHUViPSBiotin))
huvips$TotalHUViPSNoBiotin <- as.numeric(as.character(huvips$TotalHUViPSNoBiotin))

huvips$UnusedHUViPSBiotin <- na.zero(huvips$UnusedHUViPSBiotin)
huvips$UnusedHUViPSNoBiotin <- na.zero(huvips$UnusedHUViPSNoBiotin)
huvips$TotalHUViPSBiotin <- na.zero(huvips$TotalHUViPSBiotin)
huvips$TotalHUViPSNoBiotin <- na.zero(huvips$TotalHUViPSNoBiotin)

huvdf1 <- data.frame(huvipsnobiotin$Gene, huvipsnobiotin$Line)
colnames(huvdf1) <- c("Gene", "Line")
huvno1 <- nrow(huvdf1)
huvdf2 <- data.frame(huvipsbiotin$Gene, huvipsbiotin$Line)
colnames(huvdf2) <- c("Gene", "Line")
huvno2 <- nrow(huvdf2)

huvdf3 <- setdiff(huvdf1, huvdf2)
huvdf4 <- setdiff(huvdf2, huvdf1)

huvdf5 <- semi_join(huvdf1, huvdf2)
huvno5 <- nrow(huvdf5)

grid.newpage()
draw.pairwise.venn(huvno1, huvno2, huvno5, category = c("HUViPS4F5 No Biotin", "HUViPS4F5 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("yellow", "dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

huvipsbiotinunique <- anti_join(huvdf2, huvdf1)
write.csv(huvipsbiotinunique, "UniqueHUViPSBiotin.csv")


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
huvipstotal <- filter(huvips, TotalHUViPSNoBiotin < TotalHUViPSBiotin)

huvipsvsmda <- merge.data.frame(mdatotal, huvipstotal, by = "Gene", all = TRUE)

huvipsvsmda$UnusedMDABiotin <- na.zero(huvipsvsmda$UnusedMDABiotin)
huvipsvsmda$UnusedMDANoBiotin <- na.zero(huvipsvsmda$UnusedMDANoBiotin)
huvipsvsmda$TotalMDABiotin <- na.zero(huvipsvsmda$TotalMDABiotin)
huvipsvsmda$TotalMDANoBiotin <- na.zero(huvipsvsmda$TotalMDANoBiotin)
huvipsvsmda$UnusedHUViPSBiotin <- na.zero(huvipsvsmda$UnusedHUViPSBiotin)
huvipsvsmda$UnusedHUViPSNoBiotin <- na.zero(huvipsvsmda$UnusedHUViPSNoBiotin)
huvipsvsmda$TotalHUViPSBiotin <- na.zero(huvipsvsmda$TotalHUViPSBiotin)
huvipsvsmda$TotalHUViPSNoBiotin <- na.zero(huvipsvsmda$TotalHUViPSNoBiotin)

mdadf6 <- data.frame(unique(mdatotal$Gene))
colnames(mdadf6) <- c("Gene")
mdano6 <- nrow(mdadf6)
huvdf6 <- data.frame(unique(huvipstotal$Gene))
colnames(huvdf6) <- c("Gene")
huvno6 <- nrow(huvdf6)

huvmdadf1 <- setdiff(mdadf6, huvdf6)
huvmdadf2 <- setdiff(huvdf6, mdadf6)
huvmdadf3 <- semi_join(mdadf6, huvdf6)
huvmdano3 <- nrow(unique(huvmdadf3))

grid.newpage()
draw.pairwise.venn(mdano6, huvno6, huvmdano3, category = c("MDA-MB-231", "HUViPS4F1"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("light green", "dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mdaunique <- anti_join(mdadf6, huvdf6)
huvipsunique <- anti_join(huvdf6, mdadf6)
mdahuvips <- intersect(mdadf6$Gene, huvdf6$Gene)
write.csv(mdaunique, "UniqueMDAtoHUVIPS.csv")
write.csv(huvipsunique, "UniqueHUViPStoMDA.csv")
write.csv(mdahuvips, "HUViPSMDAOverlap.csv")
write.csv(mdatotal$Gene, "MDAReal.csv")
write.csv(huvipstotal$Gene, "HUViPSReal.csv")
