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
splitcolumns <- str_split_fixed(mcf7$Gene.Names," ", 2)
mcf7 <- cbind.data.frame(splitcolumns, mcf7)
names(mcf7)[names(mcf7) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(mcf7)[names(mcf7) == '2'] <- 'Synonyms'  #Change the first column to "Line"
mcf7 <- subset(mcf7, select=-Gene.Names)

#Sorting NO DOX in different Biotin categories
mcf7NDslim <- filter(mcf7, Dox == "No Dox") #Filtering only No Dox hits
mcf7NDslim <- data.frame(mcf7NDslim$Cell, mcf7NDslim$Gene, mcf7NDslim$Synonyms, mcf7NDslim$Biotin, mcf7NDslim$Unused, mcf7NDslim$Total, mcf7NDslim$Gene.Ontology)
colnames(mcf7NDslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

mcf7NDbiotin <- filter(mcf7NDslim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(mcf7NDbiotin)[colnames(mcf7NDbiotin) == 'Unused'] <- 'UnusedMCF7NDBiotin'
colnames(mcf7NDbiotin)[colnames(mcf7NDbiotin) == 'Total'] <- 'TotalMCF7NDBiotin'

mcf7NDnobiotin <- filter(mcf7NDslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(mcf7NDnobiotin)[colnames(mcf7NDnobiotin) == 'Unused'] <- 'UnusedMCF7NDNoBiotin'
colnames(mcf7NDnobiotin)[colnames(mcf7NDnobiotin) == 'Total'] <- 'TotalMCF7NDNoBiotin'

mcf7ND <- merge.data.frame(mcf7NDbiotin, mcf7NDnobiotin, by = "Gene", all = TRUE)
mcf7ND$UnusedMCF7NDBiotin <- as.numeric(as.character(mcf7ND$UnusedMCF7NDBiotin))
mcf7ND$UnusedMCF7NDNoBiotin <- as.numeric(as.character(mcf7ND$UnusedMCF7NDNoBiotin))
mcf7ND$TotalMCF7NDBiotin <- as.numeric(as.character(mcf7ND$TotalMCF7NDBiotin))
mcf7ND$TotalMCF7NDNoBiotin <- as.numeric(as.character(mcf7ND$TotalMCF7NDNoBiotin))

mcf7ND$UnusedMCF7NDBiotin <- na.zero(mcf7ND$UnusedMCF7NDBiotin)
mcf7ND$UnusedMCF7NDNoBiotin <- na.zero(mcf7ND$UnusedMCF7NDNoBiotin)
mcf7ND$TotalMCF7NDBiotin <- na.zero(mcf7ND$TotalMCF7NDBiotin)
mcf7ND$TotalMCF7NDNoBiotin <- na.zero(mcf7ND$TotalMCF7NDNoBiotin)

mcf7NDdf1 <- data.frame(mcf7NDnobiotin$Gene, mcf7NDnobiotin$Line)
colnames(mcf7NDdf1) <- c("Gene", "Line")
mcf7NDno1 <- nrow(mcf7NDdf1)
mcf7NDdf2 <- data.frame(mcf7NDbiotin$Gene, mcf7NDbiotin$Line)
colnames(mcf7NDdf2) <- c("Gene", "Line")
mcf7NDno2 <- nrow(mcf7NDdf2)

mcf7NDdf3 <- setdiff(mcf7NDdf1, mcf7NDdf2)
mcf7NDdf4 <- setdiff(mcf7NDdf2, mcf7NDdf1)
mcf7NDdf5 <- semi_join(mcf7NDdf1, mcf7NDdf2)
mcf7NDno5 <- nrow(mcf7NDdf5)

grid.newpage()
draw.pairwise.venn(mcf7NDno1, mcf7NDno2, mcf7NDno5, category = c("MCF7 No Dox/No Biotin", "MCF7 No Dox/+ Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("magenta", "purple"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mcf7NDbiotinunique <- anti_join(mcf7NDdf2, mcf7NDdf1)
write.csv(mcf7NDbiotinunique, "UniqueMCF7NDBiotin.csv")

#Sorting 3D+DOX in different Biotin categories
mcf73dslim <- filter(mcf7, Dox == "3d+Dox") #Filtering only No Dox hits
mcf73dslim <- data.frame(mcf73dslim$Cell, mcf73dslim$Gene, mcf73dslim$Synonyms, mcf73dslim$Biotin, mcf73dslim$Unused, mcf73dslim$Total, mcf73dslim$Gene.Ontology)
colnames(mcf73dslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

mcf73dbiotin <- filter(mcf73dslim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(mcf73dbiotin)[colnames(mcf73dbiotin) == 'Unused'] <- 'UnusedMCF73dBiotin'
colnames(mcf73dbiotin)[colnames(mcf73dbiotin) == 'Total'] <- 'TotalMCF73dBiotin'

mcf73dnobiotin <- filter(mcf73dslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(mcf73dnobiotin)[colnames(mcf73dnobiotin) == 'Unused'] <- 'UnusedMCF73dNoBiotin'
colnames(mcf73dnobiotin)[colnames(mcf73dnobiotin) == 'Total'] <- 'TotalMCF73dNoBiotin'

mcf73d <- merge.data.frame(mcf73dbiotin, mcf73dnobiotin, by = "Gene", all = TRUE)
mcf73d$UnusedMCF73dBiotin <- as.numeric(as.character(mcf73d$UnusedMCF73dBiotin))
mcf73d$UnusedMCF73dNoBiotin <- as.numeric(as.character(mcf73d$UnusedMCF73dNoBiotin))
mcf73d$TotalMCF73dBiotin <- as.numeric(as.character(mcf73d$TotalMCF73dBiotin))
mcf73d$TotalMCF73dNoBiotin <- as.numeric(as.character(mcf73d$TotalMCF73dNoBiotin))

mcf73d$UnusedMCF73dBiotin <- na.zero(mcf73d$UnusedMCF73dBiotin)
mcf73d$UnusedMCF73dNoBiotin <- na.zero(mcf73d$UnusedMCF73dNoBiotin)
mcf73d$TotalMCF73dBiotin <- na.zero(mcf73d$TotalMCF73dBiotin)
mcf73d$TotalMCF73dNoBiotin <- na.zero(mcf73d$TotalMCF73dNoBiotin)

mcf73ddf1 <- data.frame(mcf73dnobiotin$Gene, mcf73dnobiotin$Line)
colnames(mcf73ddf1) <- c("Gene", "Line")
mcf73dno1 <- nrow(mcf73ddf1)
mcf73ddf2 <- data.frame(mcf73dbiotin$Gene, mcf73dbiotin$Line)
colnames(mcf73ddf2) <- c("Gene", "Line")
mcf73dno2 <- nrow(mcf73ddf2)

mcf73ddf3 <- setdiff(mcf73ddf1, mcf73ddf2)
mcf73ddf4 <- setdiff(mcf73ddf2, mcf73ddf1)
mcf73ddf5 <- semi_join(mcf73ddf1, mcf73ddf2)
mcf73dno5 <- nrow(mcf73ddf5)

grid.newpage()
draw.pairwise.venn(mcf73dno1, mcf73dno2, mcf73dno5, category = c("MCF7 No Dox/No Biotin", "MCF7 3d+Dox/+ Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("yellow", "gold"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mcf73dbiotinunique <- anti_join(mcf73ddf2, mcf73ddf1)
write.csv(mcf73dbiotinunique, "UniqueMCF73dBiotin.csv")


# To compare DOX to NO DOX
mcf7NDtotal <- filter(mcf7ND, TotalMCF7NDNoBiotin < TotalMCF7NDBiotin)
mcf73dtotal <- filter(mcf73d, TotalMCF73dNoBiotin < TotalMCF73dBiotin)

mcf7NDvs3d <- merge.data.frame(mcf7NDtotal, mcf73dtotal, by = "Gene", all = TRUE)

mcf7NDvs3d$UnusedMCF7NDBiotin <- na.zero(mcf7NDvs3d$UnusedMCF7NDBiotin)
mcf7NDvs3d$UnusedMCF7NDNoBiotin <- na.zero(mcf7NDvs3d$UnusedMCF7NDNoBiotin)
mcf7NDvs3d$TotalMCF7NDBiotin <- na.zero(mcf7NDvs3d$TotalMCF7NDBiotin)
mcf7NDvs3d$TotalMCF7NDNoBiotin <- na.zero(mcf7NDvs3d$TotalMCF7NDNoBiotin)
mcf7NDvs3d$UnusedMCF73dBiotin <- na.zero(mcf7NDvs3d$UnusedMCF73dBiotin)
mcf7NDvs3d$UnusedMCF73dNoBiotin <- na.zero(mcf7NDvs3d$UnusedMCF73dNoBiotin)
mcf7NDvs3d$TotalMCF73dBiotin <- na.zero(mcf7NDvs3d$TotalMCF73dBiotin)
mcf7NDvs3d$TotalMCF73dNoBiotin <- na.zero(mcf7NDvs3d$TotalMCF73dNoBiotin)

mcf7NDdf6 <- data.frame(mcf7NDtotal$Gene)
colnames(mcf7NDdf6) <- c("Gene")
mcf7NDno6 <- nrow(mcf7NDdf6)
mcf73ddf6 <- data.frame(mcf73dtotal$Gene)
colnames(mcf73ddf6) <- c("Gene")
mcf73dno6 <- nrow(mcf73ddf6)

mcf7NDvs3ddf1 <- setdiff(mcf7NDdf6, mcf73ddf6)
mcf7NDvs3ddf2 <- setdiff(mcf73ddf6, mcf7NDdf6)
mcf7NDvs3ddf3 <- semi_join(mcf7NDdf6, mcf73ddf6)
mcf7NDvs3dno3 <- nrow(mcf7NDvs3ddf3)

grid.newpage()
draw.pairwise.venn(mcf7NDno6, mcf73dno6, mcf7NDvs3dno3, category = c("MCF7 No Dox", "MCF7 3d+Dox"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("purple", "gold"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mcf7NDunique <- anti_join(mcf7NDdf6, mcf73ddf6)
mcf73dunique <- anti_join(mcf73ddf6, mcf7NDdf6)
mcf7overlap <- intersect(mcf7NDdf6$Gene, mcf73ddf6$Gene)
write.csv(mcf7NDunique, "UniqueMCF7NDtoDox.csv")
write.csv(mcf73dunique, "UniqueMCF7DoxtoND.csv")
write.csv(mcf7overlap, "MCF7ND3dOverlap.csv")
write.csv(mcf7NDtotal$Gene, "MCF7NDReal.csv")
write.csv(mcf73dtotal$Gene, "MCF73dReal.csv")

