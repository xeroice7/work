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
fipstotal <- filter(fips, TotalFiPSNoBiotin < TotalFiPSBiotin)

fipsvsnci <- merge.data.frame(ncitotal, fipstotal, by = "Gene", all = TRUE)

fipsvsnci$UnusedNCIBiotin <- na.zero(fipsvsnci$UnusedNCIBiotin)
fipsvsnci$UnusedNCINoBiotin <- na.zero(fipsvsnci$UnusedNCINoBiotin)
fipsvsnci$TotalNCIBiotin <- na.zero(fipsvsnci$TotalNCIBiotin)
fipsvsnci$TotalNCINoBiotin <- na.zero(fipsvsnci$TotalNCINoBiotin)
fipsvsnci$UnusedFiPSBiotin <- na.zero(fipsvsnci$UnusedFiPSBiotin)
fipsvsnci$UnusedFiPSNoBiotin <- na.zero(fipsvsnci$UnusedFiPSNoBiotin)
fipsvsnci$TotalFiPSBiotin <- na.zero(fipsvsnci$TotalFiPSBiotin)
fipsvsnci$TotalFiPSNoBiotin <- na.zero(fipsvsnci$TotalFiPSNoBiotin)

ncidf6 <- data.frame(unique(ncitotal$Gene))
colnames(ncidf6) <- c("Gene")
ncino6 <- nrow(ncidf6)
fipsdf6 <- data.frame(unique(fipstotal$Gene))
colnames(fipsdf6) <- c("Gene")
fipsno6 <- nrow(fipsdf6)

fipsncidf1 <- setdiff(ncidf6, fipsdf6)
fipsncidf2 <- setdiff(fipsdf6, ncidf6)
fipsncidf3 <- semi_join(ncidf6, fipsdf6)
fipsncino3 <- nrow(unique(fipsncidf3))

grid.newpage()
draw.pairwise.venn(ncino6, fipsno6, fipsncino3, category = c("NCI-H1299", "FiPS4F5"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("orange", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

nciunique <- anti_join(ncidf6, fipsdf6)
fipsunique <- anti_join(fipsdf6, ncidf6)
ncifips <- intersect(ncidf6$Gene, fipsdf6$Gene)
write.csv(nciunique, "UniqueNCItoFIPS.csv")
write.csv(fipsunique, "UniqueFiPStoNCI.csv")
write.csv(ncifips, "FiPSNCIOverlap.csv")
write.csv(ncitotal$Gene, "NCIReal.csv")
write.csv(fipstotal$Gene, "FiPSReal.csv")
