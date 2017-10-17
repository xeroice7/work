### Packages needed

install.packages("dplyr")
install.packages("ggplot2")
install.packages("stringr")
install.packages("scales")
install.packages("devtools")
install.packages("cowplot")
install.packages("matrixStats")
install.packages("tidyverse")
install.packages("forcats")
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(reshape2)
library(devtools)
library(cowplot)
library(matrixStats)
library(tidyverse)
library(forcats)
### Bring in the data

file.path <- ("~/Desktop/Clay/qRT/SCAssay/")
scassay <- data.frame(read.csv(paste(file.path, "MCF7IPSCSCssGRP78FinalAnalysis-R.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
scassay1 <- data.frame(read.csv(paste(file.path, "MCF7CSCssGRP78toMCF7total.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

###### SCASSAY1 - TO MCF7 Heatmap

colnames(scassay1) <- c("Gene", "CD24-/CD44+", "CD24+/CD44+", "sGRP78+", "sGRP78-")
scassay1 <- melt(scassay1)
colnames(scassay1) <- c("Gene", "Population", "Relative Gene Expression")
scassay1$Gene <- as.factor(scassay1$Gene)
scassay1$Gene <- fct_rev(scassay1$Gene)


ggplot(scassay1) +
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) + 
  scale_fill_gradient(trans = 'log', breaks = c(0.1, 1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9, size = 10), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8))




#Data Manipulation - TO Non-CSCs and sGRP78-

scassay$iPS..FiPS.HUViPS. <- as.numeric(scassay$iPS..FiPS.HUViPS.)
scassay$MCF7.CD24.CD44.NonCSCs <- as.numeric(scassay$MCF7.CD24.CD44.NonCSCs)
scassay$MCF7.CD24.CD44.CSCs <- as.numeric(scassay$MCF7.CD24.CD44.CSCs)
scassay$MCF7.sGRP78.Neg <- as.numeric(scassay$MCF7.sGRP78.Neg)
scassay$MCF7.sGRP78.Pos <- as.numeric(scassay$MCF7.sGRP78.Pos)

sGRP78PtosGRP78N = 2^(-(scassay$MCF7.sGRP78.Pos - scassay$MCF7.sGRP78.Neg)) 
sGRP78PtoiPS = 2^(-(scassay$MCF7.sGRP78.Pos-scassay$iPS..FiPS.HUViPS.))
CSCtoNonCSC = 2^(-(scassay$MCF7.CD24.CD44.CSCs-scassay$MCF7.CD24.CD44.NonCSCs))
CSCtoiPS = 2^(-(scassay$MCF7.CD24.CD44.CSCs-scassay$iPS..FiPS.HUViPS.))

#Heat map generation of CSC (to Non-CSC) and sGRP78+ (to sGRP78-)
map <- cbind.data.frame(scassay$Gene, sGRP78PtosGRP78N, CSCtoNonCSC)
colnames(map) <- c("gene", "sGRP78", "CSC")
map1 <- melt(map)
colnames(map1) <- c("Gene", "Population", "Relative Gene Expression")

map1$Gene <- factor(map1$Gene, levels = rev(levels(map1$Gene)))

ggplot(map1) +
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) + 
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500), low = "white", high = "darkred", name = "") + 
  theme_classic()


#Correlations
lm_eqn = function(m) {
  
  l <- list(#a = format(coef(m)[1], digits = 2),
            #b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 4));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(~~italic(R)^2~"="~r2,l)
  } else {
    eq <- substitute(~~italic(R)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}
#If I want the equation too, add these before the "~~italic" for the if, and then the else
#italic(y) == a + b %.% italic(x)*","
#italic(y) == a - b %.% italic(x)*","

#sGRP78+ to CSCs
o <- ggplot(scassay, aes(x=MCF7.CD24.CD44.CSCs, y=MCF7.sGRP78.Pos)) +
  geom_point() + 
  geom_smooth(method = 'lm', level = 0) + #substitute level= 0.95 for 95% CI
  geom_text(aes(x = 20, y = 7, label = lm_eqn(lm(MCF7.sGRP78.Pos ~ MCF7.CD24.CD44.CSCs, scassay))), parse = TRUE, size = 8) + 
  xlim(5,25) + 
  ylim(5,25) + 
  ylab("Relative CTs (sGRP78+)") + 
  xlab("Relative CTs (CD24-/CD44+)") +
  ggtitle("sGRP78+ vs CD24-/CD44+") + 
  theme(axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6), 
        plot.title = element_text(hjust = 0.5)) + 
  theme_classic()

#sGRP78+ to Non-CSCs
p <- ggplot(scassay, aes(x=MCF7.CD24.CD44.NonCSCs, y=MCF7.sGRP78.Pos)) +
  geom_point() + 
  geom_smooth(method = 'lm', level = 0) + #substitute level= 0.95 for 95% CI 
  geom_text(aes(x = 20, y = 5, label = lm_eqn(lm(MCF7.sGRP78.Pos ~ MCF7.CD24.CD44.NonCSCs, scassay))), parse = TRUE) + 
  xlim(5,25) + 
  ylim(5,25) + 
  ylab("Relative CTs (sGRP78+)") + 
  xlab("Relative CTs (CD24+/CD44+)") +
  ggtitle("sGRP78+ vs CD24+/CD44+") + 
  theme(axis.title.x = element_text(size=6), 
        axis.title.y = element_text(size=6)) + 
  theme_classic()

#sGRP78- to CSCs
q <- ggplot(scassay, aes(x=MCF7.CD24.CD44.CSCs, y=MCF7.sGRP78.Neg)) +
  geom_point() + 
  geom_smooth(method = 'lm', level = 0) + #substitute level= 0.95 for 95% CI 
  geom_text(aes(x = 20, y = 7, label = lm_eqn(lm(MCF7.sGRP78.Neg ~ MCF7.CD24.CD44.CSCs, scassay))), parse = TRUE, size = 8) + 
  xlim(5,25) + 
  ylim(5,25) + 
  ylab("Relative CTs (sGRP78-)") + 
  xlab("Relative CTs (CD24-/CD44+)") +
  ggtitle("sGRP78- vs CD24-/CD44+") +
  theme(axis.title.x = element_text(size=6), 
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust = 0.5)) + 
  theme_classic()

#sGRP78- to Non-CSCs
r <- ggplot(scassay, aes(x=MCF7.CD24.CD44.NonCSCs, y=MCF7.sGRP78.Neg)) +
  geom_point() + 
  geom_smooth(method = 'lm', level=0.95) + 
  geom_text(aes(x = 20, y = 5, label = lm_eqn(lm(MCF7.sGRP78.Neg ~ MCF7.CD24.CD44.NonCSCs, scassay))), parse = TRUE) + 
  xlim(5,25) + 
  ylim(5,25) + 
  ylab("Relative CTs (sGRP78-)") + 
  xlab("Relative CTs (CD24+/CD44+)") +
  ggtitle("sGRP78- vs CD24+/CD44+") +
  theme(axis.title.x = element_text(size=6), 
        axis.title.y = element_text(size=6)) + 
  theme_classic()

#sGRP78+ to sGRP78-
s <- ggplot(scassay, aes(x=MCF7.sGRP78.Pos, y=MCF7.sGRP78.Neg)) +
  geom_point() + 
  geom_smooth(method = 'lm', level = 0) + #substitute level= 0.95 for 95% CI 
  geom_text(aes(x = 20, y = 7, label = lm_eqn(lm(MCF7.sGRP78.Neg ~ MCF7.sGRP78.Pos, scassay))), parse = TRUE, size = 8) + 
  xlim(5,25) + 
  ylim(5,25) + 
  ylab("Relative CTs (sGRP78-)") + 
  xlab("Relative CTs (sGRP78+)") +
  ggtitle("sGRP78- vs sGRP78+") +
  theme(axis.title.x = element_text(size=6), 
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust = 0.5)) + 
  theme_classic()

#CSCs to Non-CSCs
t <- ggplot(scassay, aes(x=MCF7.CD24.CD44.NonCSCs, y=MCF7.CD24.CD44.CSCs)) +
  geom_point() + 
  geom_smooth(method = 'lm', level = 0) + #substitute level= 0.95 for 95% CI 
  geom_text(aes(x = 20, y = 5, label = lm_eqn(lm(MCF7.CD24.CD44.CSCs ~ MCF7.CD24.CD44.NonCSCs, scassay))), parse = TRUE) + 
  xlim(5,25) + 
  ylim(5,25) + 
  ylab("Relative CTs (CD24-/CD44+)") + 
  xlab("Relative CTs (CD24+/CD44+)") +
  ggtitle("CD24-/CD44+ vs CD24+/CD44+") +
  theme(axis.title.x = element_text(size=6), 
        axis.title.y = element_text(size=6)) + 
  theme_classic()

plot_grid(o, p, q, r, s, t, ncol=2)
plot_grid(o,q,s, ncol = 1)
#sGRP78+ to iPS
u <- ggplot(scassay, aes(x=iPS..FiPS.HUViPS., y=MCF7.sGRP78.Pos)) +
  geom_point() + 
  geom_smooth(method = 'lm', level=0.95) + 
  geom_text(aes(x = 25, y = 5, label = lm_eqn(lm(MCF7.sGRP78.Pos ~ iPS..FiPS.HUViPS., scassay))), parse = TRUE) + 
  xlim(5,30) + 
  ylim(5,30) + 
  ylab("Relative CTs (sGRP78+)") + 
  xlab("Relative CTs (iPS)") +
  ggtitle("sGRP78+ to iPS") +
  theme(axis.title.x = element_text(size=6), 
        axis.title.y = element_text(size=6)) + 
  theme_classic()

#sGRP78- to iPS
v <- ggplot(scassay, aes(x=iPS..FiPS.HUViPS., y=MCF7.sGRP78.Neg)) +
  geom_point() + 
  geom_smooth(method = 'lm', level=0.95) + 
  geom_text(aes(x = 25, y = 5, label = lm_eqn(lm(MCF7.sGRP78.Neg ~ iPS..FiPS.HUViPS., scassay))), parse = TRUE) + 
  xlim(5,30) + 
  ylim(5,30) + 
  ylab("Relative CTs (sGRP78-)") + 
  xlab("Relative CTs (iPS)") +
  ggtitle("sGRP78- to iPS") +
  theme(axis.title.x = element_text(size=6), 
        axis.title.y = element_text(size=6)) + 
  theme_classic()

#CSCs to iPS
w <- ggplot(scassay, aes(x=iPS..FiPS.HUViPS., y=MCF7.CD24.CD44.CSCs)) +
  geom_point() + 
  geom_smooth(method = 'lm', level=0.95) + 
  geom_text(aes(x = 25, y = 5, label = lm_eqn(lm(MCF7.CD24.CD44.CSCs ~ iPS..FiPS.HUViPS., scassay))), parse = TRUE) + 
  xlim(5,30) + 
  ylim(5,30) + 
  ylab("Relative CTs (CD24-/CD44+)") + 
  xlab("Relative CTs (iPS)") +
  ggtitle("CD24-/CD44+ to iPS") +
  theme(axis.title.x = element_text(size=6), 
        axis.title.y = element_text(size=6)) + 
  theme_classic()

#Non-CSCs to iPS
x <- ggplot(scassay, aes(x=iPS..FiPS.HUViPS., y=MCF7.CD24.CD44.NonCSCs)) +
  geom_point() + 
  geom_smooth(method = 'lm', level=0.95) + 
  geom_text(aes(x = 25, y = 5, label = lm_eqn(lm(MCF7.CD24.CD44.NonCSCs ~ iPS..FiPS.HUViPS., scassay))), parse = TRUE) + 
  xlim(5,30) + 
  ylim(5,30) + 
  ylab("Relative CTs (CD24+/CD44+)") + 
  xlab("Relative CTs (iPS)") +
  ggtitle("CD24+/CD44+ to iPS") +
  theme(axis.title.x = element_text(size=6), 
        axis.title.y = element_text(size=6)) + 
  theme_classic()


plot_grid(o, p, q, r, s, t, u, v, w, x, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"), ncol=3)

#Faceting based on gene ontology

facetDF <- data.frame(read.csv(paste(file.path, "FacetsMCF7R.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

# Regulation of cell cycle
regulationDF <- subset(facetDF, facetDF$Regulation.of.Cell.Cycle == 'Y')
regulationDF <- cbind.data.frame(regulationDF$Gene, regulationDF$Non.CSCs, regulationDF$CSCs, regulationDF$sGRP78Neg, regulationDF$sGRP78Pos)
colnames(regulationDF) <- c("Gene", "Non-CSCs", "CSCs", "sGRP78-", "sGRP78+")
regulationDF <- melt(regulationDF)
colnames(regulationDF) <- c("Gene", "Population", "Relative Gene Expression")
regulationDF$Gene <- factor(regulationDF$Gene, levels = rev(levels(regulationDF$Gene)))
regulationDF <- filter(regulationDF, Population != "Non-CSCs")


a <- ggplot(regulationDF) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  ggtitle("Regulation of Cell Cycle") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

# Chromatin Modification and Remodeling

chromDF <- subset(facetDF, facetDF$Chromatin.Modification.and.Remodeling.Factors == 'Y')
chromDF <- cbind.data.frame(chromDF$Gene, chromDF$Non.CSCs, chromDF$CSCs, chromDF$sGRP78Neg, chromDF$sGRP78Pos)
colnames(chromDF) <- c("Gene", "Non-CSCs", "CSCs", "sGRP78-", "sGRP78+")
chromDF <- melt(chromDF)
colnames(chromDF) <- c("Gene", "Population", "Relative Gene Expression")
chromDF$Gene <- factor(chromDF$Gene, levels = rev(levels(chromDF$Gene)))
chromDF <- filter(chromDF, Population != "Non-CSCs")

b <- ggplot(chromDF) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  ggtitle("Chromatin Remodeling and Modification") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

# Symmetric and Asymmetric Cell Division

symDF <- subset(facetDF, facetDF$Symmetric.and.Asymmetric.Cell.Division == 'Y')
symDF <- cbind.data.frame(symDF$Gene, symDF$Non.CSCs, symDF$CSCs, symDF$sGRP78Neg, symDF$sGRP78Pos)
colnames(symDF) <- c("Gene", "Non-CSCs", "CSCs", "sGRP78-", "sGRP78+")
symDF <- melt(symDF)
colnames(symDF) <- c("Gene", "Population", "Relative Gene Expression")
symDF$Gene <- factor(symDF$Gene, levels = rev(levels(symDF$Gene)))
symDF <- filter(symDF, Population != "Non-CSCs")

c <- ggplot(symDF) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  ggtitle("Symmetric and Asymmetric Cell Division") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

# Pluripotency (formerly Self-Renewal)

selfDF <- subset(facetDF, facetDF$Self.Renewal.Markers == 'Y')
selfDF <- cbind.data.frame(selfDF$Gene, selfDF$Non.CSCs, selfDF$CSCs, selfDF$sGRP78Neg, selfDF$sGRP78Pos)
colnames(selfDF) <- c("Gene", "Non-CSCs", "CSCs", "sGRP78-", "sGRP78+")
selfDF <- melt(selfDF)
colnames(selfDF) <- c("Gene", "Population", "Relative Gene Expression")
selfDF$Gene <- factor(selfDF$Gene, levels = rev(levels(selfDF$Gene)))
selfDF <- filter(selfDF, Population != "Non-CSCs")

d <- ggplot(selfDF) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  ggtitle("Pluripotency") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

# Cytokines/Growth Factors

cytoDF <- subset(facetDF, facetDF$Cytokines.Growth.Factors == 'Y')
cytoDF <- cbind.data.frame(cytoDF$Gene, cytoDF$Non.CSCs, cytoDF$CSCs, cytoDF$sGRP78Neg, cytoDF$sGRP78Pos)
colnames(cytoDF) <- c("Gene", "Non-CSCs", "CSCs", "sGRP78-", "sGRP78+")
cytoDF <- melt(cytoDF)
colnames(cytoDF) <- c("Gene", "Population", "Relative Gene Expression")
cytoDF$Gene <- factor(cytoDF$Gene, levels = rev(levels(cytoDF$Gene)))
cytoDF <- filter(cytoDF, Population != "Non-CSCs")

e <- ggplot(cytoDF) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  ggtitle("Cytokines/Growth Factors") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

# Cell-Cell Communication

cellcomDF <- subset(facetDF, facetDF$Cell.Cell.Communication == 'Y')
cellcomDF <- cbind.data.frame(cellcomDF$Gene, cellcomDF$Non.CSCs, cellcomDF$CSCs, cellcomDF$sGRP78Neg, cellcomDF$sGRP78Pos)
colnames(cellcomDF) <- c("Gene", "Non-CSCs", "CSCs", "sGRP78-", "sGRP78+")
cellcomDF <- melt(cellcomDF)
colnames(cellcomDF) <- c("Gene", "Population", "Relative Gene Expression")
cellcomDF$Gene <- factor(cellcomDF$Gene, levels = rev(levels(cellcomDF$Gene)))
cellcomDF <- filter(cellcomDF, Population != "Non-CSCs")

f <- ggplot(cellcomDF) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  ggtitle("Cell-Cell Communication") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

# Cell Adhesion 

celladDF <- subset(facetDF, facetDF$Cell.Adhesion.Molecules == 'Y')
celladDF <- cbind.data.frame(celladDF$Gene, celladDF$Non.CSCs, celladDF$CSCs, celladDF$sGRP78Neg, celladDF$sGRP78Pos)
colnames(celladDF) <- c("Gene", "Non-CSCs", "CSCs", "sGRP78-", "sGRP78+")
celladDF <- melt(celladDF)
colnames(celladDF) <- c("Gene", "Population", "Relative Gene Expression")
celladDF$Gene <- factor(celladDF$Gene, levels = rev(levels(celladDF$Gene)))
celladDF <- filter(celladDF, Population != "Non-CSCs")

g <- ggplot(celladDF) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  ggtitle("Cell Adhesion") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

# Metabolic Markers

metaDF <- subset(facetDF, facetDF$Metabolic.Markers == 'Y')
metaDF <- cbind.data.frame(metaDF$Gene, metaDF$Non.CSCs, metaDF$CSCs, metaDF$sGRP78Neg, metaDF$sGRP78Pos)
colnames(metaDF) <- c("Gene", "Non-CSCs", "CSCs", "sGRP78-", "sGRP78+")
metaDF <- melt(metaDF)
colnames(metaDF) <- c("Gene", "Population", "Relative Gene Expression")
metaDF$Gene <- factor(metaDF$Gene, levels = rev(levels(metaDF$Gene)))
metaDF <- filter(metaDF, Population != "Non-CSCs")

h <- ggplot(metaDF) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  ggtitle("Metabolism") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

# Embryonic Cell Lineage Markers

#embDF <- subset(facetDF, facetDF$Embryonic.Cell.Lineage.Markers == 'Y')
#embDF <- cbind.data.frame(embDF$Gene, embDF$CSC, embDF$sGRP78Pos)
#colnames(embDF) <- c("Gene", "CSC", "sGRP78+")
#embDF <- melt(embDF)
#colnames(embDF) <- c("Gene", "Population", "Relative Gene Expression")
#embDF$Gene <- factor(embDF$Gene, levels = rev(levels(embDF$Gene)))

#i <- ggplot(embDF) + 
  #geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  #scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500), low = "white", high = "darkred", name = "") + 
  #ggtitle("Embryonic Cell Lineage") + 
  #xlab("") + 
  #ylab("") + 
  #theme_classic()

# Hematopoietic Cell Lineage Markers

#hemaDF <- subset(facetDF, facetDF$Hematopoietic.Cell.Lineage.Markers== 'Y')
#hemaDF <- cbind.data.frame(hemaDF$Gene, hemaDF$CSC, hemaDF$sGRP78Pos)
#colnames(hemaDF) <- c("Gene", "CSC", "sGRP78+")
#hemaDF <- melt(hemaDF)
#colnames(hemaDF) <- c("Gene", "Population", "Relative Gene Expression")
#hemaDF$Gene <- factor(hemaDF$Gene, levels = rev(levels(hemaDF$Gene)))

#j <- ggplot(hemaDF) + 
  #geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  #scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500), low = "white", high = "darkred", name = "") + 
  #ggtitle("Hematopoietic Cell Lineage") + 
  #xlab("") + 
  #ylab("") + 
  #theme_classic()

# Mesenchymal Cell Lineage Markers

#mesDF <- subset(facetDF, facetDF$Mesenchymal.Cell.Lineage.Markers == 'Y')
#mesDF <- cbind.data.frame(mesDF$Gene, mesDF$CSC, mesDF$sGRP78Pos)
#colnames(mesDF) <- c("Gene", "CSC", "sGRP78+")
#mesDF <- melt(mesDF)
#colnames(mesDF) <- c("Gene", "Population", "Relative Gene Expression")
#cellcomDF$Gene <- factor(mesDF$Gene, levels = rev(levels(mesDF$Gene)))

#k <- ggplot(mesDF) + 
  #geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  #scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500), low = "white", high = "darkred", name = "") + 
  #ggtitle("Mesenchymal Cell Lineage") + 
  #xlab("") + 
  #ylab("") + 
  #theme_classic()

# Neural Cell Lineage Markers

#neuralDF <- subset(facetDF, facetDF$Neural.Cell.Lineage.Markers== 'Y')
#neuralDF <- cbind.data.frame(neuralDF$Gene, neuralDF$CSC, neuralDF$sGRP78Pos)
#colnames(neuralDF) <- c("Gene", "CSC", "sGRP78+")
#neuralDF <- melt(neuralDF)
#colnames(neuralDF) <- c("Gene", "Population", "Relative Gene Expression")
#neuralDF$Gene <- factor(neuralDF$Gene, levels = rev(levels(neuralDF$Gene)))

#l <- ggplot(neuralDF) + 
  #geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  #scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500), low = "white", high = "darkred", name = "") + 
  #ggtitle("Neural Cell Lineage") + 
  #xlab("") + 
  #ylab("") + 
  #theme_classic()

# Notch Signaling

notchDF <- subset(facetDF, facetDF$Notch.Signaling == 'Y')
notchDF <- cbind.data.frame(notchDF$Gene, notchDF$Non.CSCs, notchDF$CSCs, notchDF$sGRP78Neg, notchDF$sGRP78Pos)
colnames(notchDF) <- c("Gene", "Non-CSCs", "CSCs", "sGRP78-", "sGRP78+")
notchDF <- melt(notchDF)
colnames(notchDF) <- c("Gene", "Population", "Relative Gene Expression")
notchDF$Gene <- factor(notchDF$Gene, levels = rev(levels(notchDF$Gene)))
notchDF <- filter(notchDF, Population != "Non-CSCs")

m <- ggplot(notchDF) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  ggtitle("Notch Signaling") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

# Wnt Signaling

wntDF <- subset(facetDF, facetDF$Wnt.Signaling == 'Y')
wntDF <- cbind.data.frame(wntDF$Gene, wntDF$Non.CSCs, wntDF$CSCs, wntDF$sGRP78Neg, wntDF$sGRP78Pos)
colnames(wntDF) <- c("Gene", "Non-CSCs", "CSCs", "sGRP78-", "sGRP78+")
wntDF <- melt(wntDF)
colnames(wntDF) <- c("Gene", "Population", "Relative Gene Expression")
wntDF$Gene <- factor(wntDF$Gene, levels = rev(levels(wntDF$Gene)))
wntDF <- filter(wntDF, Population != "Non-CSCs")

n <- ggplot(wntDF) + 
  geom_tile(aes(x = Population, y = Gene, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks = c(1, 2, 5, 10, 20, 50, 100, 300, 500, 1000), low = "white", high = "darkred", name = "") + 
  ggtitle("Wnt Signaling") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

#plot_grid(a, b, c, d, e, f, g, h, i, j, k, l, m, n, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"), ncol=3)
plot_grid(a, b, c, d, e, f, g, h, m, n, ncol=2)


#sGRP78+ expression in CSCs

grp <- data.frame(read.csv(paste(file.path, "MCF7CSCsSCAssay.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

grp$Mean <- rowMeans(grp[,2:3])
grp <- transform(grp, SD=apply(grp[,2:3], 1, sd, na.rm=TRUE))
grp$X <- as.factor(grp$X)
grp$X <- ordered(grp$X, levels = c("MCF7", "CD24+/CD44+", "CD24-/CD44+"))

ggplot(grp, aes(x=X, y=Mean)) + 
  geom_bar(aes(fill = X), stat = "identity", color = "black") + 
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), width = 0.3) + 
  scale_fill_brewer(palette = "Set1") + 
  ggtitle("Normalized GRP78 Expression") + 
  ylab("Relative Gene Expression") + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        legend.position = "none")
