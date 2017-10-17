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
splitcolumns <- str_split_fixed(mda$Gene.Names," ", 2)
mda <- cbind.data.frame(splitcolumns, mda)
names(mda)[names(mda) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(mda)[names(mda) == '2'] <- 'Synonyms'  #Change the first column to "Line"
mda <- subset(mda, select=-Gene.Names)

#Sorting NO DOX in different Biotin categories
mdaNDslim <- filter(mda, Dox == "No Dox") #Filtering only No Dox hits
mdaNDslim <- data.frame(mdaNDslim$Cell, mdaNDslim$Gene, mdaNDslim$Synonyms, mdaNDslim$Biotin, mdaNDslim$Unused, mdaNDslim$Total, mdaNDslim$Gene.Ontology)
colnames(mdaNDslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

mdaNDbiotin <- filter(mdaNDslim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(mdaNDbiotin)[colnames(mdaNDbiotin) == 'Unused'] <- 'UnusedMDANDBiotin'
colnames(mdaNDbiotin)[colnames(mdaNDbiotin) == 'Total'] <- 'TotalMDANDBiotin'

mdaNDnobiotin <- filter(mdaNDslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(mdaNDnobiotin)[colnames(mdaNDnobiotin) == 'Unused'] <- 'UnusedMDANDNoBiotin'
colnames(mdaNDnobiotin)[colnames(mdaNDnobiotin) == 'Total'] <- 'TotalMDANDNoBiotin'

mdaND <- merge.data.frame(mdaNDbiotin, mdaNDnobiotin, by = "Gene", all = TRUE)
mdaND$UnusedMDANDBiotin <- as.numeric(as.character(mdaND$UnusedMDANDBiotin))
mdaND$UnusedMDANDNoBiotin <- as.numeric(as.character(mdaND$UnusedMDANDNoBiotin))
mdaND$TotalMDANDBiotin <- as.numeric(as.character(mdaND$TotalMDANDBiotin))
mdaND$TotalMDANDNoBiotin <- as.numeric(as.character(mdaND$TotalMDANDNoBiotin))

mdaND$UnusedMDANDBiotin <- na.zero(mdaND$UnusedMDANDBiotin)
mdaND$UnusedMDANDNoBiotin <- na.zero(mdaND$UnusedMDANDNoBiotin)
mdaND$TotalMDANDBiotin <- na.zero(mdaND$TotalMDANDBiotin)
mdaND$TotalMDANDNoBiotin <- na.zero(mdaND$TotalMDANDNoBiotin)

mdaNDdf1 <- data.frame(mdaNDnobiotin$Gene, mdaNDnobiotin$Line)
colnames(mdaNDdf1) <- c("Gene", "Line")
mdaNDno1 <- nrow(mdaNDdf1)
mdaNDdf2 <- data.frame(mdaNDbiotin$Gene, mdaNDbiotin$Line)
colnames(mdaNDdf2) <- c("Gene", "Line")
mdaNDno2 <- nrow(mdaNDdf2)

mdaNDdf3 <- setdiff(mdaNDdf1, mdaNDdf2)
mdaNDdf4 <- setdiff(mdaNDdf2, mdaNDdf1)
mdaNDdf5 <- semi_join(mdaNDdf1, mdaNDdf2)
mdaNDno5 <- nrow(mdaNDdf5)

grid.newpage()
draw.pairwise.venn(mdaNDno1, mdaNDno2, mdaNDno5, category = c("MDA No Dox/No Biotin", "MDA No Dox/+ Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("brown", "tan"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mdaNDbiotinunique <- anti_join(mdaNDdf2, mdaNDdf1)
write.csv(mdaNDbiotinunique, "UniqueMDANDBiotin.csv")

#Sorting 3D+DOX in different Biotin categories
mda3dslim <- filter(mda, Dox == "3d+Dox") #Filtering only No Dox hits
mda3dslim <- data.frame(mda3dslim$Cell, mda3dslim$Gene, mda3dslim$Synonyms, mda3dslim$Biotin, mda3dslim$Unused, mda3dslim$Total, mda3dslim$Gene.Ontology)
colnames(mda3dslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

mda3dbiotin <- filter(mda3dslim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(mda3dbiotin)[colnames(mda3dbiotin) == 'Unused'] <- 'UnusedMDA3dBiotin'
colnames(mda3dbiotin)[colnames(mda3dbiotin) == 'Total'] <- 'TotalMDA3dBiotin'

mda3dnobiotin <- filter(mda3dslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(mda3dnobiotin)[colnames(mda3dnobiotin) == 'Unused'] <- 'UnusedMDA3dNoBiotin'
colnames(mda3dnobiotin)[colnames(mda3dnobiotin) == 'Total'] <- 'TotalMDA3dNoBiotin'

mda3d <- merge.data.frame(mda3dbiotin, mda3dnobiotin, by = "Gene", all = TRUE)
mda3d$UnusedMDA3dBiotin <- as.numeric(as.character(mda3d$UnusedMDA3dBiotin))
mda3d$UnusedMDA3dNoBiotin <- as.numeric(as.character(mda3d$UnusedMDA3dNoBiotin))
mda3d$TotalMDA3dBiotin <- as.numeric(as.character(mda3d$TotalMDA3dBiotin))
mda3d$TotalMDA3dNoBiotin <- as.numeric(as.character(mda3d$TotalMDA3dNoBiotin))

mda3d$UnusedMDA3dBiotin <- na.zero(mda3d$UnusedMDA3dBiotin)
mda3d$UnusedMDA3dNoBiotin <- na.zero(mda3d$UnusedMDA3dNoBiotin)
mda3d$TotalMDA3dBiotin <- na.zero(mda3d$TotalMDA3dBiotin)
mda3d$TotalMDA3dNoBiotin <- na.zero(mda3d$TotalMDA3dNoBiotin)

mda3ddf1 <- data.frame(mda3dnobiotin$Gene, mda3dnobiotin$Line)
colnames(mda3ddf1) <- c("Gene", "Line")
mda3dno1 <- nrow(mda3ddf1)
mda3ddf2 <- data.frame(mda3dbiotin$Gene, mda3dbiotin$Line)
colnames(mda3ddf2) <- c("Gene", "Line")
mda3dno2 <- nrow(mda3ddf2)

mda3ddf3 <- setdiff(mda3ddf1, mda3ddf2)
mda3ddf4 <- setdiff(mda3ddf2, mda3ddf1)
mda3ddf5 <- semi_join(mda3ddf1, mda3ddf2)
mda3dno5 <- nrow(mda3ddf5)

grid.newpage()
draw.pairwise.venn(mda3dno1, mda3dno2, mda3dno5, category = c("MDA 3d+Dox/No Biotin", "MDA 3d+Dox/+ Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("wheat", "snow"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mda3dbiotinunique <- anti_join(mda3ddf2, mda3ddf1)
write.csv(mda3dbiotinunique, "UniqueMDA3dBiotin.csv")


# To compare DOX to NO DOX
mdaNDtotal <- filter(mdaND, TotalMDANDNoBiotin < TotalMDANDBiotin)
mda3dtotal <- filter(mda3d, TotalMDA3dNoBiotin < TotalMDA3dBiotin)

mdaNDvs3d <- merge.data.frame(mdaNDtotal, mda3dtotal, by = "Gene", all = TRUE)

mdaNDvs3d$UnusedMDANDBiotin <- na.zero(mdaNDvs3d$UnusedMDANDBiotin)
mdaNDvs3d$UnusedMDANDNoBiotin <- na.zero(mdaNDvs3d$UnusedMDANDNoBiotin)
mdaNDvs3d$TotalMDANDBiotin <- na.zero(mdaNDvs3d$TotalMDANDBiotin)
mdaNDvs3d$TotalMDANDNoBiotin <- na.zero(mdaNDvs3d$TotalMDANDNoBiotin)
mdaNDvs3d$UnusedMDA3dBiotin <- na.zero(mdaNDvs3d$UnusedMDA3dBiotin)
mdaNDvs3d$UnusedMDA3dNoBiotin <- na.zero(mdaNDvs3d$UnusedMDA3dNoBiotin)
mdaNDvs3d$TotalMDA3dBiotin <- na.zero(mdaNDvs3d$TotalMDA3dBiotin)
mdaNDvs3d$TotalMDA3dNoBiotin <- na.zero(mdaNDvs3d$TotalMDA3dNoBiotin)

mdaNDdf6 <- data.frame(mdaNDtotal$Gene)
colnames(mdaNDdf6) <- c("Gene")
mdaNDno6 <- nrow(mdaNDdf6)
mda3ddf6 <- data.frame(mda3dtotal$Gene)
colnames(mda3ddf6) <- c("Gene")
mda3dno6 <- nrow(mda3ddf6)

mdaNDvs3ddf1 <- setdiff(mdaNDdf6, mda3ddf6)
mdaNDvs3ddf2 <- setdiff(mda3ddf6, mdaNDdf6)
mdaNDvs3ddf3 <- semi_join(mdaNDdf6, mda3ddf6)
mdaNDvs3dno3 <- nrow(mdaNDvs3ddf3)

grid.newpage()
draw.pairwise.venn(mdaNDno6, mda3dno6, mdaNDvs3dno3, category = c("MDA No Dox", "MDA 3d+Dox"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("tan", "wheat"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mdaNDunique <- anti_join(mdaNDdf6, mda3ddf6)
mda3dunique <- anti_join(mda3ddf6, mdaNDdf6)
mdaoverlap <- intersect(mdaNDdf6$Gene, mda3ddf6$Gene)
write.csv(mdaNDunique, "UniqueMDANDtoDox.csv")
write.csv(mda3dunique, "UniqueMDADoxtoND.csv")
write.csv(mdaoverlap, "MDAND3dOverlap.csv")
write.csv(mdaNDtotal$Gene, "MDANDReal.csv")
write.csv(mda3dtotal$Gene, "MDA3dReal.csv")
