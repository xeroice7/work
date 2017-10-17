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
names(huv)[names(huv) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(huv)[names(huv) == '2'] <- 'Synonyms'  #Change the first column to "Line"
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
huvipstotal <- filter(huvips, TotalHUViPSNoBiotin < TotalHUViPSBiotin)

huvipsvsa549 <- merge.data.frame(a549total, huvipstotal, by = "Gene", all = TRUE)

huvipsvsa549$UnusedA549Biotin <- na.zero(huvipsvsa549$UnusedA549Biotin)
huvipsvsa549$UnusedA549NoBiotin <- na.zero(huvipsvsa549$UnusedA549NoBiotin)
huvipsvsa549$TotalA549Biotin <- na.zero(huvipsvsa549$TotalA549Biotin)
huvipsvsa549$TotalA549NoBiotin <- na.zero(huvipsvsa549$TotalA549NoBiotin)
huvipsvsa549$UnusedHUViPSBiotin <- na.zero(huvipsvsa549$UnusedHUViPSBiotin)
huvipsvsa549$UnusedHUViPSNoBiotin <- na.zero(huvipsvsa549$UnusedHUViPSNoBiotin)
huvipsvsa549$TotalHUViPSBiotin <- na.zero(huvipsvsa549$TotalHUViPSBiotin)
huvipsvsa549$TotalHUViPSNoBiotin <- na.zero(huvipsvsa549$TotalHUViPSNoBiotin)

a549df6 <- data.frame(a549total$Gene)
colnames(a549df6) <- c("Gene")
a549no6 <- nrow(a549df6)
huvdf6 <- data.frame(huvipstotal$Gene)
colnames(huvdf6) <- c("Gene")
huvno6 <- nrow(huvdf6)

huva549df1 <- setdiff(a549df6, huvdf6)
huva549df2 <- setdiff(huvdf6, a549df6)
huva549df3 <- semi_join(a549df6, huvdf6)
huva549no3 <- nrow(unique(huva549df3))

grid.newpage()
draw.pairwise.venn(a549no6, huvno6, huva549no3, category = c("A549", "HUViPS4F1"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark blue", "dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

a549unique <- anti_join(a549df6, huvdf6)
huvipsunique <- anti_join(huvdf6, a549df6)
a549huvips <- intersect(a549df6$Gene, huvdf6$Gene)
write.csv(a549unique, "UniqueA549toHUVIPS.csv")
write.csv(huvipsunique, "UniqueHUViPStoA549.csv")
write.csv(huva549df3, "HUViPSA549Overlap.csv")
write.csv(a549total$Gene, "A549Real.csv")
write.csv(huvipstotal$Gene, "HUViPSReal.csv")

