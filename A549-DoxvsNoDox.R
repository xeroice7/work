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
splitcolumns <- str_split_fixed(a549$Gene.Names," ", 2)
a549 <- cbind.data.frame(splitcolumns, a549)
names(a549)[names(a549) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(a549)[names(a549) == '2'] <- 'Synonyms'  #Change the first column to "Line"
a549 <- subset(a549, select=-Gene.Names)

#Sorting NO DOX in different Biotin categories
a549NDslim <- filter(a549, Dox == "No") #Filtering only No Dox hits
a549NDslim <- data.frame(a549NDslim$Cell, a549NDslim$Gene, a549NDslim$Synonyms, a549NDslim$Biotin, a549NDslim$Unused, a549NDslim$Total, a549NDslim$Gene.Ontology)
colnames(a549NDslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

a549NDbiotin <- filter(a549NDslim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(a549NDbiotin)[colnames(a549NDbiotin) == 'Unused'] <- 'UnusedA549NDBiotin'
colnames(a549NDbiotin)[colnames(a549NDbiotin) == 'Total'] <- 'TotalA549NDBiotin'

a549NDnobiotin <- filter(a549NDslim, Biotin == "N") #Keep only No Biotin for FiPS
colnames(a549NDnobiotin)[colnames(a549NDnobiotin) == 'Unused'] <- 'UnusedA549NDNoBiotin'
colnames(a549NDnobiotin)[colnames(a549NDnobiotin) == 'Total'] <- 'TotalA549NDNoBiotin'

a549ND <- merge.data.frame(a549NDbiotin, a549NDnobiotin, by = "Gene", all = TRUE)
a549ND$UnusedA549NDBiotin <- as.numeric(as.character(a549ND$UnusedA549NDBiotin))
a549ND$UnusedA549NDNoBiotin <- as.numeric(as.character(a549ND$UnusedA549NDNoBiotin))
a549ND$TotalA549NDBiotin <- as.numeric(as.character(a549ND$TotalA549NDBiotin))
a549ND$TotalA549NDNoBiotin <- as.numeric(as.character(a549ND$TotalA549NDNoBiotin))

a549ND$UnusedA549NDBiotin <- na.zero(a549ND$UnusedA549NDBiotin)
a549ND$UnusedA549NDNoBiotin <- na.zero(a549ND$UnusedA549NDNoBiotin)
a549ND$TotalA549NDBiotin <- na.zero(a549ND$TotalA549NDBiotin)
a549ND$TotalA549NDNoBiotin <- na.zero(a549ND$TotalA549NDNoBiotin)

a549NDdf1 <- data.frame(a549NDnobiotin$Gene, a549NDnobiotin$Line)
colnames(a549NDdf1) <- c("Gene", "Line")
a549NDno1 <- nrow(a549NDdf1)
a549NDdf2 <- data.frame(a549NDbiotin$Gene, a549NDbiotin$Line)
colnames(a549NDdf2) <- c("Gene", "Line")
a549NDno2 <- nrow(a549NDdf2)

a549NDdf3 <- setdiff(a549NDdf1, a549NDdf2)
a549NDdf4 <- setdiff(a549NDdf2, a549NDdf1)
a549NDdf5 <- semi_join(a549NDdf1, a549NDdf2)
a549NDno5 <- nrow(a549NDdf5)

grid.newpage()
draw.pairwise.venn(a549NDno1, a549NDno2, a549NDno5, category = c("A549 No Dox/No Biotin", "A549 No Dox/+ Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark red", "red"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

a549NDbiotinunique <- anti_join(a549NDdf2, a549NDdf1)
write.csv(a549NDbiotinunique, "UniqueA549NDBiotin.csv")

#Sorting DOX in different Biotin categories
a5493dslim <- filter(a549, Dox == "3d+Dox") #Filtering only No Dox hits
a5493dslim <- data.frame(a5493dslim$Cell, a5493dslim$Gene, a5493dslim$Synonyms, a5493dslim$Biotin, a5493dslim$Unused, a5493dslim$Total, a5493dslim$Gene.Ontology)
colnames(a5493dslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

a5493dbiotin <- filter(a5493dslim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(a5493dbiotin)[colnames(a5493dbiotin) == 'Unused'] <- 'UnusedA5493dBiotin'
colnames(a5493dbiotin)[colnames(a5493dbiotin) == 'Total'] <- 'TotalA5493dBiotin'

a5493dnobiotin <- filter(a5493dslim, Biotin == "N") #Keep only No Biotin for FiPS
colnames(a5493dnobiotin)[colnames(a5493dnobiotin) == 'Unused'] <- 'UnusedA5493dNoBiotin'
colnames(a5493dnobiotin)[colnames(a5493dnobiotin) == 'Total'] <- 'TotalA5493dNoBiotin'

a5493d <- merge.data.frame(a5493dbiotin, a5493dnobiotin, by = "Gene", all = TRUE)
a5493d$UnusedA5493dBiotin <- as.numeric(as.character(a5493d$UnusedA5493dBiotin))
a5493d$UnusedA5493dNoBiotin <- as.numeric(as.character(a5493d$UnusedA5493dNoBiotin))
a5493d$TotalA5493dBiotin <- as.numeric(as.character(a5493d$TotalA5493dBiotin))
a5493d$TotalA5493dNoBiotin <- as.numeric(as.character(a5493d$TotalA5493dNoBiotin))

a5493d$UnusedA5493dBiotin <- na.zero(a5493d$UnusedA5493dBiotin)
a5493d$UnusedA5493dNoBiotin <- na.zero(a5493d$UnusedA5493dNoBiotin)
a5493d$TotalA5493dBiotin <- na.zero(a5493d$TotalA5493dBiotin)
a5493d$TotalA5493dNoBiotin <- na.zero(a5493d$TotalA5493dNoBiotin)

a5493ddf1 <- data.frame(a5493dnobiotin$Gene, a5493dnobiotin$Line)
colnames(a5493ddf1) <- c("Gene", "Line")
a5493dno1 <- nrow(a5493ddf1)
a5493ddf2 <- data.frame(a5493dbiotin$Gene, a5493dbiotin$Line)
colnames(a5493ddf2) <- c("Gene", "Line")
a5493dno2 <- nrow(a5493ddf2)

a5493ddf3 <- setdiff(a5493ddf1, a5493ddf2)
a5493ddf4 <- setdiff(a5493ddf2, a5493ddf1)
a5493ddf5 <- semi_join(a5493ddf1, a5493ddf2)
a5493dno5 <- nrow(a5493ddf5)

grid.newpage()
draw.pairwise.venn(a5493dno1, a5493dno2, a5493dno5, category = c("A549 3d+Dox/No Biotin", "A549 3d+Dox/+ Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark blue", "blue"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

a5493dbiotinunique <- anti_join(a5493ddf2, a5493ddf1)
write.csv(a5493dbiotinunique, "UniqueA5493dBiotin.csv")


# To compare DOX to NO DOX
a549NDtotal <- filter(a549ND, TotalA549NDNoBiotin < TotalA549NDBiotin)
a5493dtotal <- filter(a5493d, TotalA5493dNoBiotin < TotalA5493dBiotin)

a549NDvs3d <- merge.data.frame(a549NDtotal, a5493dtotal, by = "Gene", all = TRUE)

a549NDvs3d$UnusedA549NDBiotin <- na.zero(a549NDvs3d$UnusedA549NDBiotin)
a549NDvs3d$UnusedA549NDNoBiotin <- na.zero(a549NDvs3d$UnusedA549NDNoBiotin)
a549NDvs3d$TotalA549NDBiotin <- na.zero(a549NDvs3d$TotalA549NDBiotin)
a549NDvs3d$TotalA549NDNoBiotin <- na.zero(a549NDvs3d$TotalA549NDNoBiotin)
a549NDvs3d$UnusedA5493dBiotin <- na.zero(a549NDvs3d$UnusedA5493dBiotin)
a549NDvs3d$UnusedA5493dNoBiotin <- na.zero(a549NDvs3d$UnusedA5493dNoBiotin)
a549NDvs3d$TotalA5493dBiotin <- na.zero(a549NDvs3d$TotalA5493dBiotin)
a549NDvs3d$TotalA5493dNoBiotin <- na.zero(a549NDvs3d$TotalA5493dNoBiotin)

a549NDdf6 <- data.frame(a549NDtotal$Gene)
colnames(a549NDdf6) <- c("Gene")
a549NDno6 <- nrow(a549NDdf6)
a5493ddf6 <- data.frame(a5493dtotal$Gene)
colnames(a5493ddf6) <- c("Gene")
a5493dno6 <- nrow(a5493ddf6)

a549NDvs3ddf1 <- setdiff(a549NDdf6, a5493ddf6)
a549NDvs3ddf2 <- setdiff(a5493ddf6, a549NDdf6)
a549NDvs3ddf3 <- semi_join(a549NDdf6, a5493ddf6)
a549NDvs3dno3 <- nrow(a549NDvs3ddf3)

grid.newpage()
draw.pairwise.venn(a549NDno6, a5493dno6, a549NDvs3dno3, category = c("A549 No Dox", "A549 3d+Dox"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("red", "blue"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

a549NDunique <- anti_join(a549NDdf6, a5493ddf6)
a5493dunique <- anti_join(a5493ddf6, a549NDdf6)
a549overlap <- intersect(a549NDdf6$Gene, a5493ddf6$Gene)
write.csv(a549NDunique, "UniqueA549NDtoDox.csv")
write.csv(a5493dunique, "UniqueA549DoxtoND.csv")
write.csv(a549overlap, "A549ND3dOverlap.csv")
write.csv(a549NDtotal$Gene, "A549NDReal.csv")
write.csv(a5493dtotal$Gene, "A5493dReal.csv")

