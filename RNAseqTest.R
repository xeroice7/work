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


file.path <- ("~/Desktop/")
#rna <- data.frame(read.csv(paste(file.path, "HeatMapCSV.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = TRUE))
#rna <- data.frame(read.csv(paste(file.path, "HeatMap10000.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = TRUE))
rna <- data.frame(read.csv(paste(file.path, "HeatMap1000.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = TRUE))
rna$Cell.Line <- factor(rna$Cell.Line, levels = rev(as.character(rna[[1]])))
rnamelt <- melt(rna)
colnames(rnamelt) <- c("Cell Line", "X", "Gene", "Relative Gene Expression")
rnamelt$`Cell Line` <- factor(rnamelt$`Cell Line`, levels = as.character(rna[[1]]))
rnamelt$`Cell Line` <- with(rnamelt, factor(`Cell Line`, levels = rev(levels(`Cell Line`))))

rnamelt1 <- subset(rnamelt, rnamelt$Gene == "end.OCT4")
rnamelt2 <- subset(rnamelt, rnamelt$Gene == "end.SOX2")
rnamelt3 <- subset(rnamelt, rnamelt$Gene == "NANOG")

# Version 1 - Plot_Grid
aa <- ggplot(rnamelt1) + 
  geom_tile(aes(x = Gene, y = `Cell Line`, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', low = "white", high = "firebrick1", name = "") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 8, angle=90), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

bb <- ggplot(rnamelt2) + 
  geom_tile(aes(x = Gene, y = `Cell Line`, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', low = "white", high = "firebrick1", name = "") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 8, angle=90), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

cc <- ggplot(rnamelt3) + 
  geom_tile(aes(x = Gene, y = `Cell Line`, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', low = "white", high = "firebrick1", name = "") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 8, angle=90), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

plot_grid(aa, bb, cc, ncol=3)

# Version 2 - Facets
ggplot(rnamelt) + 
  facet_wrap(~ Gene, scales = "free") +
  geom_tile(aes(x = Gene, y = `Cell Line`, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', low = "white", high = "firebrick1", name = "") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 8, angle=90), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 

#Version 3 - Same Log Scale
ggplot(rnamelt) + 
  geom_tile(aes(x = Gene, y = `Cell Line`, fill = `Relative Gene Expression`)) +  
  scale_fill_gradient(trans = 'log', breaks=c(0,1,10,100,1000,10000), low = "white", high = "firebrick1", name = "") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8)) 
