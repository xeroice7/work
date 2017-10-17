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
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
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
names(fips)[names(fips) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(fips)[names(fips) == '2'] <- 'Synonyms'  #Change the first column to "Line"
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


#A549
a549 <- data.frame(read.csv(paste(file.path, "A549SummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(a549$X, "/", 2)
a549 <- cbind.data.frame(splitcolumns, a549)
names(a549)[names(a549) == '1'] <- 'Line'  #Change the first column to "Line"
names(a549)[names(a549) == '2'] <- 'YN'  #Change the second column to "YN"
splitcolumns <- str_split_fixed(a549$Line, " ", 3)
a549 <- cbind.data.frame(splitcolumns, a549)
names(a549)[names(a549) == '1'] <- 'Cell'  #Change the first column to "Cell"
names(a549)[names(a549) == '2'] <- 'Dox'  #Change the first column to "Line"
splitcolumns <- str_split_fixed(a549$YN, "", 2)
splitcolumns <- splitcolumns[,1]
a549 <- cbind.data.frame(splitcolumns, a549)
names(a549)[names(a549) == "splitcolumns"] <- "Biotin"
a549 <- subset(a549, select=-c(`3`, Line, YN, X))
a549 <- filter(a549, Species == "HUMAN") #Filtering only human hits
a549 <- filter(a549, Dox == "No") #Filtering only No Dox hits
splitcolumns <- str_split_fixed(a549$Gene.Names," ", 2)
a549 <- cbind.data.frame(splitcolumns, a549)
names(a549)[names(a549) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(a549)[names(a549) == '2'] <- 'Synonyms'  #Change the first column to "Line"
a549 <- subset(a549, select=-Gene.Names)
a549slim <- data.frame(a549$Cell, a549$Gene, a549$Synonyms, a549$Biotin, a549$Unused, a549$Total, a549$Gene.Ontology)
colnames(a549slim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

a549biotin <- filter(a549slim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(a549biotin)[colnames(a549biotin) == 'Unused'] <- 'UnusedA549Biotin'
colnames(a549biotin)[colnames(a549biotin) == 'Total'] <- 'TotalA549Biotin'

a549nobiotin <- filter(a549slim, Biotin == "N") #Keep only No Biotin for FiPS
colnames(a549nobiotin)[colnames(a549nobiotin) == 'Unused'] <- 'UnusedA549NoBiotin'
colnames(a549nobiotin)[colnames(a549nobiotin) == 'Total'] <- 'TotalA549NoBiotin'

a549 <- merge.data.frame(a549biotin, a549nobiotin, by = "Gene", all = TRUE)
a549$UnusedA549Biotin <- as.numeric(as.character(a549$UnusedA549Biotin))
a549$UnusedA549NoBiotin <- as.numeric(as.character(a549$UnusedA549NoBiotin))
a549$TotalA549Biotin <- as.numeric(as.character(a549$TotalA549Biotin))
a549$TotalA549NoBiotin <- as.numeric(as.character(a549$TotalA549NoBiotin))

a549$UnusedA549Biotin <- na.zero(a549$UnusedA549Biotin)
a549$UnusedA549NoBiotin <- na.zero(a549$UnusedA549NoBiotin)
a549$TotalA549Biotin <- na.zero(a549$TotalA549Biotin)
a549$TotalA549NoBiotin <- na.zero(a549$TotalA549NoBiotin)

a549df1 <- data.frame(a549nobiotin$Gene, a549nobiotin$Line)
colnames(a549df1) <- c("Gene", "Line")
a549no1 <- nrow(a549df1)
a549df2 <- data.frame(a549biotin$Gene, a549biotin$Line)
colnames(a549df2) <- c("Gene", "Line")
a549no2 <- nrow(a549df2)

a549df3 <- setdiff(a549df1, a549df2)
a549df4 <- setdiff(a549df2, a549df1)
a549df5 <- semi_join(a549df1, a549df2)
a549no5 <- nrow(a549df5)

grid.newpage()
draw.pairwise.venn(a549no1, a549no2, a549no5, category = c("A549 No Biotin", "A549 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark red", "dark blue"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

a549biotinunique <- anti_join(a549df2, a549df1)
write.csv(a549biotinunique, "UniqueA549Biotin.csv")



# To compare together
a549total <- filter(a549, TotalA549NoBiotin < TotalA549Biotin)
fipstotal <- filter(fips, TotalFiPSNoBiotin < TotalFiPSBiotin)

fipsvsa549 <- merge.data.frame(a549total, fipstotal, by = "Gene", all = TRUE)

fipsvsa549$UnusedA549Biotin <- na.zero(fipsvsa549$UnusedA549Biotin)
fipsvsa549$UnusedA549NoBiotin <- na.zero(fipsvsa549$UnusedA549NoBiotin)
fipsvsa549$TotalA549Biotin <- na.zero(fipsvsa549$TotalA549Biotin)
fipsvsa549$TotalA549NoBiotin <- na.zero(fipsvsa549$TotalA549NoBiotin)
fipsvsa549$UnusedFiPSBiotin <- na.zero(fipsvsa549$UnusedFiPSBiotin)
fipsvsa549$UnusedFiPSNoBiotin <- na.zero(fipsvsa549$UnusedFiPSNoBiotin)
fipsvsa549$TotalFiPSBiotin <- na.zero(fipsvsa549$TotalFiPSBiotin)
fipsvsa549$TotalFiPSNoBiotin <- na.zero(fipsvsa549$TotalFiPSNoBiotin)

a549df6 <- data.frame(unique(a549total$Gene))
colnames(a549df6) <- c("Gene")
a549no6 <- nrow(a549df6)
fipsdf6 <- data.frame(unique(fipstotal$Gene))
colnames(fipsdf6) <- c("Gene")
fipsno6 <- nrow(fipsdf6)

fipsa549df1 <- setdiff(a549df6, fipsdf6)
fipsa549df2 <- setdiff(fipsdf6, a549df6)
fipsa549df3 <- semi_join(a549df6, fipsdf6)
fipsa549no3 <- nrow(unique(fipsa549df3))

grid.newpage()
draw.pairwise.venn(a549no6, fipsno6, fipsa549no3, category = c("A549", "FiPS4F5"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark blue", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

a549unique <- anti_join(a549df6, fipsdf6)
fipsunique <- anti_join(fipsdf6, a549df6)
a549fips <- intersect(a549df6$Gene, fipsdf6$Gene)
write.csv(a549unique, "UniqueA549toFIPS.csv")
write.csv(fipsunique, "UniqueFiPStoA549.csv")
write.csv(a549fips, "FiPSA549Overlap.csv")
write.csv(a549total$Gene, "A549Real.csv")
write.csv(fipstotal$Gene, "FiPSReal.csv")

