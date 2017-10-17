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


#MCF7
mcf7 <- data.frame(read.csv(paste(file.path, "MCF7SummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(mcf7$X, "/", 2)
mcf7 <- cbind.data.frame(splitcolumns, mcf7)
names(mcf7)[names(mcf7) == '1'] <- 'Line'  #Change the first column to "Line"
names(mcf7)[names(mcf7) == '2'] <- 'YN'  #Change the second column to "YN"
splitcolumns <- str_split_fixed(mcf7$Line, " ", 3)
mcf7 <- cbind.data.frame(splitcolumns, mcf7)
names(mcf7)[names(mcf7) == '1'] <- 'Cell'  #Change the first column to "Cell"
names(mcf7)[names(mcf7) == '3'] <- 'Dox'  #Change the first column to "Line"
splitcolumns <- str_split_fixed(mcf7$YN, " ", 2)
splitcolumns <- splitcolumns[,1]
mcf7 <- cbind.data.frame(splitcolumns, mcf7)
names(mcf7)[names(mcf7) == "splitcolumns"] <- "Biotin"
mcf7 <- subset(mcf7, select=-c(`2`, Line, YN, X))
mcf7 <- filter(mcf7, Species == "HUMAN") #Filtering only human hits
mcf7 <- filter(mcf7, Dox == "No Dox") #Filtering only No Dox hits
splitcolumns <- str_split_fixed(mcf7$Gene.Names," ", 2)
mcf7 <- cbind.data.frame(splitcolumns, mcf7)
names(mcf7)[names(mcf7) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(mcf7)[names(mcf7) == '2'] <- 'Synonyms'  #Change the first column to "Line"
mcf7 <- subset(mcf7, select=-Gene.Names)
mcf7slim <- data.frame(mcf7$Cell, mcf7$Gene, mcf7$Synonyms, mcf7$Biotin, mcf7$Unused, mcf7$Total, mcf7$Gene.Ontology)
colnames(mcf7slim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

mcf7biotin <- filter(mcf7slim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(mcf7biotin)[colnames(mcf7biotin) == 'Unused'] <- 'UnusedMCF7Biotin'
colnames(mcf7biotin)[colnames(mcf7biotin) == 'Total'] <- 'TotalMCF7Biotin'

mcf7nobiotin <- filter(mcf7slim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(mcf7nobiotin)[colnames(mcf7nobiotin) == 'Unused'] <- 'UnusedMCF7NoBiotin'
colnames(mcf7nobiotin)[colnames(mcf7nobiotin) == 'Total'] <- 'TotalMCF7NoBiotin'

mcf7 <- merge.data.frame(mcf7biotin, mcf7nobiotin, by = "Gene", all = TRUE)
mcf7$UnusedMCF7Biotin <- as.numeric(as.character(mcf7$UnusedMCF7Biotin))
mcf7$UnusedMCF7NoBiotin <- as.numeric(as.character(mcf7$UnusedMCF7NoBiotin))
mcf7$TotalMCF7Biotin <- as.numeric(as.character(mcf7$TotalMCF7Biotin))
mcf7$TotalMCF7NoBiotin <- as.numeric(as.character(mcf7$TotalMCF7NoBiotin))

mcf7$UnusedMCF7Biotin <- na.zero(mcf7$UnusedMCF7Biotin)
mcf7$UnusedMCF7NoBiotin <- na.zero(mcf7$UnusedMCF7NoBiotin)
mcf7$TotalMCF7Biotin <- na.zero(mcf7$TotalMCF7Biotin)
mcf7$TotalMCF7NoBiotin <- na.zero(mcf7$TotalMCF7NoBiotin)

mcf7df1 <- data.frame(mcf7nobiotin$Gene, mcf7nobiotin$Line)
colnames(mcf7df1) <- c("Gene", "Line")
mcf7no1 <- nrow(mcf7df1)
mcf7df2 <- data.frame(mcf7biotin$Gene, mcf7biotin$Line)
colnames(mcf7df2) <- c("Gene", "Line")
mcf7no2 <- nrow(mcf7df2)

mcf7df3 <- setdiff(mcf7df1, mcf7df2)
mcf7df4 <- setdiff(mcf7df2, mcf7df2)
mcf7df5 <- semi_join(mcf7df1, mcf7df2)
mcf7no5 <- nrow(mcf7df5)

grid.newpage()
draw.pairwise.venn(mcf7no1, mcf7no2, mcf7no5, category = c("MCF7 No Biotin", "MCF7 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("deepskyblue", "firebrick1"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mcf7biotinunique <- anti_join(mcf7df2, mcf7df1)
write.csv(mcf7biotinunique, "UniqueMCF7Biotin.csv")


# To compare together
mcf7total <- filter(mcf7, TotalMCF7NoBiotin < TotalMCF7Biotin)
fipstotal <- filter(fips, TotalFiPSNoBiotin < TotalFiPSBiotin)

fipsvsmcf7 <- merge.data.frame(mcf7total, fipstotal, by = "Gene", all = TRUE)

fipsvsmcf7$UnusedMCF7Biotin <- na.zero(fipsvsmcf7$UnusedMCF7Biotin)
fipsvsmcf7$UnusedMCF7NoBiotin <- na.zero(fipsvsmcf7$UnusedMCF7NoBiotin)
fipsvsmcf7$TotalMCF7Biotin <- na.zero(fipsvsmcf7$TotalMCF7Biotin)
fipsvsmcf7$TotalMCF7NoBiotin <- na.zero(fipsvsmcf7$TotalMCF7NoBiotin)
fipsvsmcf7$UnusedFiPSBiotin <- na.zero(fipsvsmcf7$UnusedFiPSBiotin)
fipsvsmcf7$UnusedFiPSNoBiotin <- na.zero(fipsvsmcf7$UnusedFiPSNoBiotin)
fipsvsmcf7$TotalFiPSBiotin <- na.zero(fipsvsmcf7$TotalFiPSBiotin)
fipsvsmcf7$TotalFiPSNoBiotin <- na.zero(fipsvsmcf7$TotalFiPSNoBiotin)

mcf7df6 <- data.frame(unique(mcf7total$Gene))
colnames(mcf7df6) <- c("Gene")
mcf7no6 <- nrow(mcf7df6)
fipsdf6 <- data.frame(unique(fipstotal$Gene))
colnames(fipsdf6) <- c("Gene")
fipsno6 <- nrow(fipsdf6)

fipsmcf7df1 <- setdiff(mcf7df6, fipsdf6)
fipsmcf7df2 <- setdiff(fipsdf6, mcf7df6)
fipsmcf7df3 <- semi_join(mcf7df6, fipsdf6)
fipsmcf7no3 <- nrow(unique(fipsmcf7df3))

grid.newpage()
draw.pairwise.venn(mcf7no6, fipsno6, fipsmcf7no3, category = c("MCF7", "FiPS4F5"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("firebrick1", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mcf7uniquetofips <- anti_join(mcf7df6, fipsdf6)
fipsuniquetomcf7 <- anti_join(fipsdf6, mcf7df6)
mcf7fips <- intersect(mcf7df6$Gene, fipsdf6$Gene)
write.csv(mcf7uniquetofips, "UniqueMCF7toFIPS.csv")
write.csv(fipsuniquetomcf7, "UniqueFiPStoMCF7.csv")
write.csv(mcf7fips, "FiPSMCF7Overlap.csv")
write.csv(mcf7total$Gene, "MCF7Real.csv")
write.csv(fipstotal$Gene, "FiPSReal.csv")

