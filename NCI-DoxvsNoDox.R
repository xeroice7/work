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
splitcolumns <- str_split_fixed(nci$Gene.Names," ", 2)
nci <- cbind.data.frame(splitcolumns, nci)
names(nci)[names(nci) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(nci)[names(nci) == '2'] <- 'Synonyms'  #Change the first column to "Line"
nci <- subset(nci, select=-Gene.Names)


#Sorting NO DOX in different Biotin categories
nciNDslim <- filter(nci, Dox == "No Dox") #Filtering only No Dox hits
nciNDslim <- data.frame(nciNDslim$Cell, nciNDslim$Gene, nciNDslim$Synonyms, nciNDslim$Biotin, nciNDslim$Unused, nciNDslim$Total, nciNDslim$Gene.Ontology)
colnames(nciNDslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

nciNDbiotin <- filter(nciNDslim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(nciNDbiotin)[colnames(nciNDbiotin) == 'Unused'] <- 'UnusedNCINDBiotin'
colnames(nciNDbiotin)[colnames(nciNDbiotin) == 'Total'] <- 'TotalNCINDBiotin'

nciNDnobiotin <- filter(nciNDslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(nciNDnobiotin)[colnames(nciNDnobiotin) == 'Unused'] <- 'UnusedNCINDNoBiotin'
colnames(nciNDnobiotin)[colnames(nciNDnobiotin) == 'Total'] <- 'TotalNCINDNoBiotin'

nciND <- merge.data.frame(nciNDbiotin, nciNDnobiotin, by = "Gene", all = TRUE)
nciND$UnusedNCINDBiotin <- as.numeric(as.character(nciND$UnusedNCINDBiotin))
nciND$UnusedNCINDNoBiotin <- as.numeric(as.character(nciND$UnusedNCINDNoBiotin))
nciND$TotalNCINDBiotin <- as.numeric(as.character(nciND$TotalNCINDBiotin))
nciND$TotalNCINDNoBiotin <- as.numeric(as.character(nciND$TotalNCINDNoBiotin))

nciND$UnusedNCINDBiotin <- na.zero(nciND$UnusedNCINDBiotin)
nciND$UnusedNCINDNoBiotin <- na.zero(nciND$UnusedNCINDNoBiotin)
nciND$TotalNCINDBiotin <- na.zero(nciND$TotalNCINDBiotin)
nciND$TotalNCINDNoBiotin <- na.zero(nciND$TotalNCINDNoBiotin)

nciNDdf1 <- data.frame(nciNDnobiotin$Gene, nciNDnobiotin$Line)
colnames(nciNDdf1) <- c("Gene", "Line")
nciNDno1 <- nrow(nciNDdf1)
nciNDdf2 <- data.frame(nciNDbiotin$Gene, nciNDbiotin$Line)
colnames(nciNDdf2) <- c("Gene", "Line")
nciNDno2 <- nrow(nciNDdf2)

nciNDdf3 <- setdiff(nciNDdf1, nciNDdf2)
nciNDdf4 <- setdiff(nciNDdf2, nciNDdf1)
nciNDdf5 <- semi_join(nciNDdf1, nciNDdf2)
nciNDno5 <- nrow(nciNDdf5)

grid.newpage()
draw.pairwise.venn(nciNDno1, nciNDno2, nciNDno5, category = c("NCI-H1299 No Dox/No Biotin", "NCI-H1299 No Dox/+ Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark orange", "orange"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

nciNDbiotinunique <- anti_join(nciNDdf2, nciNDdf1)
write.csv(nciNDbiotinunique, "UniqueNCINDBiotin.csv")

#Sorting 3D+DOX in different Biotin categories
nci3dslim <- filter(nci, Dox == "3d+Dox") #Filtering only No Dox hits
nci3dslim <- data.frame(nci3dslim$Cell, nci3dslim$Gene, nci3dslim$Synonyms, nci3dslim$Biotin, nci3dslim$Unused, nci3dslim$Total, nci3dslim$Gene.Ontology)
colnames(nci3dslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

nci3dbiotin <- filter(nci3dslim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(nci3dbiotin)[colnames(nci3dbiotin) == 'Unused'] <- 'UnusedNCI3dBiotin'
colnames(nci3dbiotin)[colnames(nci3dbiotin) == 'Total'] <- 'TotalNCI3dBiotin'

nci3dnobiotin <- filter(nci3dslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(nci3dnobiotin)[colnames(nci3dnobiotin) == 'Unused'] <- 'UnusedNCI3dNoBiotin'
colnames(nci3dnobiotin)[colnames(nci3dnobiotin) == 'Total'] <- 'TotalNCI3dNoBiotin'

nci3d <- merge.data.frame(nci3dbiotin, nci3dnobiotin, by = "Gene", all = TRUE)
nci3d$UnusedNCI3dBiotin <- as.numeric(as.character(nci3d$UnusedNCI3dBiotin))
nci3d$UnusedNCI3dNoBiotin <- as.numeric(as.character(nci3d$UnusedNCI3dNoBiotin))
nci3d$TotalNCI3dBiotin <- as.numeric(as.character(nci3d$TotalNCI3dBiotin))
nci3d$TotalNCI3dNoBiotin <- as.numeric(as.character(nci3d$TotalNCI3dNoBiotin))

nci3d$UnusedNCI3dBiotin <- na.zero(nci3d$UnusedNCI3dBiotin)
nci3d$UnusedNCI3dNoBiotin <- na.zero(nci3d$UnusedNCI3dNoBiotin)
nci3d$TotalNCI3dBiotin <- na.zero(nci3d$TotalNCI3dBiotin)
nci3d$TotalNCI3dNoBiotin <- na.zero(nci3d$TotalNCI3dNoBiotin)

nci3ddf1 <- data.frame(nci3dnobiotin$Gene, nci3dnobiotin$Line)
colnames(nci3ddf1) <- c("Gene", "Line")
nci3dno1 <- nrow(nci3ddf1)
nci3ddf2 <- data.frame(nci3dbiotin$Gene, nci3dbiotin$Line)
colnames(nci3ddf2) <- c("Gene", "Line")
nci3dno2 <- nrow(nci3ddf2)

nci3ddf3 <- setdiff(nci3ddf1, nci3ddf2)
nci3ddf4 <- setdiff(nci3ddf2, nci3ddf1)
nci3ddf5 <- semi_join(nci3ddf1, nci3ddf2)
nci3dno5 <- nrow(nci3ddf5)

grid.newpage()
draw.pairwise.venn(nci3dno1, nci3dno2, nci3dno5, category = c("NCI-H1299 3d+Dox/No Biotin", "NCI-H1299 3d+Dox/+ Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark green", "green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

nci3dbiotinunique <- anti_join(nci3ddf2, nci3ddf1)
write.csv(nci3dbiotinunique, "UniqueNCI3dBiotin.csv")


# To compare DOX to NO DOX
nciNDtotal <- filter(nciND, TotalNCINDNoBiotin < TotalNCINDBiotin)
nci3dtotal <- filter(nci3d, TotalNCI3dNoBiotin < TotalNCI3dBiotin)

nciNDvs3d <- merge.data.frame(nciNDtotal, nci3dtotal, by = "Gene", all = TRUE)

nciNDvs3d$UnusedNCINDBiotin <- na.zero(nciNDvs3d$UnusedNCINDBiotin)
nciNDvs3d$UnusedNCINDNoBiotin <- na.zero(nciNDvs3d$UnusedNCINDNoBiotin)
nciNDvs3d$TotalNCINDBiotin <- na.zero(nciNDvs3d$TotalNCINDBiotin)
nciNDvs3d$TotalNCINDNoBiotin <- na.zero(nciNDvs3d$TotalNCINDNoBiotin)
nciNDvs3d$UnusedNCI3dBiotin <- na.zero(nciNDvs3d$UnusedNCI3dBiotin)
nciNDvs3d$UnusedNCI3dNoBiotin <- na.zero(nciNDvs3d$UnusedNCI3dNoBiotin)
nciNDvs3d$TotalNCI3dBiotin <- na.zero(nciNDvs3d$TotalNCI3dBiotin)
nciNDvs3d$TotalNCI3dNoBiotin <- na.zero(nciNDvs3d$TotalNCI3dNoBiotin)

nciNDdf6 <- data.frame(nciNDtotal$Gene)
colnames(nciNDdf6) <- c("Gene")
nciNDno6 <- nrow(nciNDdf6)
nci3ddf6 <- data.frame(nci3dtotal$Gene)
colnames(nci3ddf6) <- c("Gene")
nci3dno6 <- nrow(nci3ddf6)

nciNDvs3ddf1 <- setdiff(nciNDdf6, nci3ddf6)
nciNDvs3ddf2 <- setdiff(nci3ddf6, nciNDdf6)
nciNDvs3ddf3 <- semi_join(nciNDdf6, nci3ddf6)
nciNDvs3dno3 <- nrow(nciNDvs3ddf3)

grid.newpage()
draw.pairwise.venn(nciNDno6, nci3dno6, nciNDvs3dno3, category = c("NCI No Dox", "NCI 3d+Dox"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("orange", "green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

nciNDunique <- anti_join(nciNDdf6, nci3ddf6)
nci3dunique <- anti_join(nci3ddf6, nciNDdf6)
ncioverlap <- intersect(nciNDdf6$Gene, nci3ddf6$Gene)
write.csv(nciNDunique, "UniqueNCINDtoDox.csv")
write.csv(nci3dunique, "UniqueNCIDoxtoND.csv")
write.csv(ncioverlap, "NCIND3dOverlap.csv")
write.csv(nciNDtotal$Gene, "NCINDReal.csv")
write.csv(nci3dtotal$Gene, "NCI3dReal.csv")

