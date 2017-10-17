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


#NCI-H1299
nci <- data.frame(read.csv(paste(file.path, "NCIH1299SummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(nci$X, "/", 2)
nci <- cbind.data.frame(splitcolumns, nci)
names(nci)[names(nci) == '1'] <- 'Line'  #Change the first column to "Line"
names(nci)[names(nci) == '2'] <- 'YN'  #Change the second column to "YN"
splitcolumns <- str_split_fixed(nci$Line, " ", 3)
nci <- cbind.data.frame(splitcolumns, nci)
names(nci)[names(nci) == '1'] <- 'Cell'  #Change the first column to "Cell"
names(nci)[names(nci) == '3'] <- 'Dox'  #Change the first column to "Line"
splitcolumns <- str_split_fixed(nci$YN, " ", 2)
splitcolumns <- splitcolumns[,1]
nci <- cbind.data.frame(splitcolumns, nci)
names(nci)[names(nci) == "splitcolumns"] <- "Biotin"
nci <- subset(nci, select=-c(`2`, Line, YN, X))
nci <- filter(nci, Species == "HUMAN") #Filtering only human hits
nci <- filter(nci, Dox == "No Dox") #Filtering only No Dox hits
splitcolumns <- str_split_fixed(nci$Gene.Names," ", 2)
nci <- cbind.data.frame(splitcolumns, nci)
names(nci)[names(nci) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(nci)[names(nci) == '2'] <- 'Synonyms'  #Change the first column to "Line"
nci <- subset(nci, select=-Gene.Names)
ncislim <- data.frame(nci$Cell, nci$Gene, nci$Synonyms, nci$Biotin, nci$Unused, nci$Total, nci$Gene.Ontology)
colnames(ncislim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

ncibiotin <- filter(ncislim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(ncibiotin)[colnames(ncibiotin) == 'Unused'] <- 'UnusedNCIBiotin'
colnames(ncibiotin)[colnames(ncibiotin) == 'Total'] <- 'TotalNCIBiotin'

ncinobiotin <- filter(ncislim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(ncinobiotin)[colnames(ncinobiotin) == 'Unused'] <- 'UnusedNCINoBiotin'
colnames(ncinobiotin)[colnames(ncinobiotin) == 'Total'] <- 'TotalNCINoBiotin'

nci <- merge.data.frame(ncibiotin, ncinobiotin, by = "Gene", all = TRUE)
nci$UnusedNCIBiotin <- as.numeric(as.character(nci$UnusedNCIBiotin))
nci$UnusedNCINoBiotin <- as.numeric(as.character(nci$UnusedNCINoBiotin))
nci$TotalNCIBiotin <- as.numeric(as.character(nci$TotalNCIBiotin))
nci$TotalNCINoBiotin <- as.numeric(as.character(nci$TotalNCINoBiotin))

nci$UnusedNCIBiotin <- na.zero(nci$UnusedNCIBiotin)
nci$UnusedNCINoBiotin <- na.zero(nci$UnusedNCINoBiotin)
nci$TotalNCIBiotin <- na.zero(nci$TotalNCIBiotin)
nci$TotalNCINoBiotin <- na.zero(nci$TotalNCINoBiotin)

ncidf1 <- data.frame(ncinobiotin$Gene, ncinobiotin$Line)
colnames(ncidf1) <- c("Gene", "Line")
ncino1 <- nrow(ncidf1)
ncidf2 <- data.frame(ncibiotin$Gene, ncibiotin$Line)
colnames(ncidf2) <- c("Gene", "Line")
ncino2 <- nrow(ncidf2)

ncidf3 <- setdiff(ncidf1, ncidf2)
ncidf4 <- setdiff(ncidf2, ncidf1)
ncidf5 <- semi_join(ncidf1, ncidf2)
ncino5 <- nrow(ncidf5)

grid.newpage()
draw.pairwise.venn(ncino1, ncino2, ncino5, category = c("NCI-H1299 No Biotin", "NCI-H1299 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("turquoise", "orange"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

ncibiotinunique <- anti_join(ncidf2, ncidf1)
write.csv(ncibiotinunique, "UniqueNCIBiotin.csv")


# To compare together
ncitotal <- filter(nci, TotalNCINoBiotin < TotalNCIBiotin)
huvipstotal <- filter(huvips, TotalHUViPSNoBiotin < TotalHUViPSBiotin)

huvipsvsnci <- merge.data.frame(ncitotal, huvipstotal, by = "Gene", all = TRUE)

huvipsvsnci$UnusedNCIBiotin <- na.zero(huvipsvsnci$UnusedNCIBiotin)
huvipsvsnci$UnusedNCINoBiotin <- na.zero(huvipsvsnci$UnusedNCINoBiotin)
huvipsvsnci$TotalNCIBiotin <- na.zero(huvipsvsnci$TotalNCIBiotin)
huvipsvsnci$TotalNCINoBiotin <- na.zero(huvipsvsnci$TotalNCINoBiotin)
huvipsvsnci$UnusedHUViPSBiotin <- na.zero(huvipsvsnci$UnusedHUViPSBiotin)
huvipsvsnci$UnusedHUViPSNoBiotin <- na.zero(huvipsvsnci$UnusedHUViPSNoBiotin)
huvipsvsnci$TotalHUViPSBiotin <- na.zero(huvipsvsnci$TotalHUViPSBiotin)
huvipsvsnci$TotalHUViPSNoBiotin <- na.zero(huvipsvsnci$TotalHUViPSNoBiotin)

ncidf6 <- data.frame(ncitotal$Gene)
colnames(ncidf6) <- c("Gene")
ncino6 <- nrow(ncidf6)
huvdf6 <- data.frame(huvipstotal$Gene)
colnames(huvdf6) <- c("Gene")
huvno6 <- nrow(huvdf6)

huvncidf1 <- setdiff(ncidf6, huvdf6)
huvncidf2 <- setdiff(huvdf6, ncidf6)
huvncidf3 <- semi_join(ncidf6, huvdf6)
huvncino3 <- nrow(huvncidf3)

grid.newpage()
draw.pairwise.venn(ncino6, huvno6, huvncino3, category = c("NCI-H1299", "HUViPS4F1"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("orange", "dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

nciunique <- anti_join(ncidf6, huvdf6)
huvipsunique <- anti_join(huvdf6, ncidf6)
ncihuvips <- intersect(ncidf6$Gene, huvdf6$Gene)
write.csv(nciunique, "UniqueNCItoHUVIPS.csv")
write.csv(huvipsunique, "UniqueHUViPStoNCI.csv")
write.csv(ncihuvips, "HUViPSNCIOverlap.csv")
write.csv(ncitotal$Gene, "NCIReal.csv")
write.csv(huvipstotal$Gene, "HUViPSReal.csv")
