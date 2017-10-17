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

file.path <- ("~/Desktop/Clay/In Vivo/")
lung <- data.frame(read.csv(paste(file.path, "InVivo2daysLung.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
brain <- data.frame(read.csv(paste(file.path, "InVivo2daysBrain.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

lung <- lung[-c(4:nrow(lung)),]
lung$Mean <- rowMeans(lung[,2:5], na.rm = TRUE)
lung <- transform(lung, SD=apply(lung[,2:5], 1, sd, na.rm=TRUE))
lung$X <- factor(lung$X, levels = c("Total MCF7", "CD24-/CD44+", "sGRP78+"))
lung$X[3] <- "sGRP78+"

z <- ggplot(lung, aes(x=X, y=Mean)) +
  geom_bar(aes(fill=X), stat = "identity", color="black") + 
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), width = 0.3) + 
  scale_fill_brewer(palette = "Set1") + 
  ggtitle("Total Lung Spots per Mouse") + 
  ylab("Average # Lung Spots/Mouse") + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        legend.position = "none")

brain <- brain[-c(4:nrow(brain)),]
brain$Mean <- rowMeans(brain[,2:5], na.rm = TRUE)
brain <- transform(brain, SD=apply(brain[,2:5], 1, sd, na.rm=TRUE))
brain$X <- factor(brain$X, levels = c("Total MCF7", "CD24-/CD44+", "sGRP78+"))
brain$X[3] <- "sGRP78+"

y <- ggplot(brain, aes(x=X, y=Mean)) +
  geom_bar(aes(fill=X), stat = "identity", color="black") + 
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), width = 0.3) + 
  scale_fill_brewer(palette = "Set1") + 
  ggtitle("Total Brain Spots per Mouse") + 
  ylab("Average # Brain Spots/Mouse") + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        legend.position = "none")


plot_grid(z, y, ncol=2)
