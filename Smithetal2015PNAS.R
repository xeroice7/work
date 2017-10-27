library(tidyverse)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ISLR)
library(gplots)

##################
# ADDITIONAL GENE CRITERIA WE WANT TO FILTER AND COMPARE OUR DATASET TO
##################

#Total iPS overlaps (including somatic targets - 93 genes)
huvfips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/HUViPSvsFiPS/HUViPSFiPSOverlap.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
huvfips_overlap <- as.factor(huvfips_overlap$x)

#All Human Cell surface-specific targets indicated by GO analysis on GO site 
data <- read.csv("~/Desktop/data.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) # Import CSV of surface expression
surface_genes <- data$V1 
surface_genes <- unique(surface_genes) 

#Unique iPS targets(no somatic source - 34 genes)
gene <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/CSVs/TotalUniqueiPSGENEandPID.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

#iPS to MCF7 Comparisons (No Somatic Filters - can include somatic hits - 29 genes) - Jun's list
mcf7huvips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsHUViPS/HUViPSMCF7OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
mcf7fips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsFiPS/FiPSMCF7OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
mcf7ips_overlap <- intersect(mcf7fips_overlap$V1, mcf7huvips_overlap$V1)

#A smaller subset of genes pulled from FiPS and HUViPS MS sets to test clustering - 390 genes
testgenes <- read.csv("~/Desktop/TestGeneSetUnique.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) # Import CSV of surface expression
testgenes <- as.factor(testgenes$x)
c_testgenes <- as.character(testgenes)
c_gene <- as.character(gene)
totaltestgenes <- c(c_testgenes, c_gene)
totaltestgenes <- unique(totaltestgenes)
totaltestgenes <- as.factor(totaltestgenes)

#Combination of Jun's list (29 genes) and my unique list (34 genes) - total of 50 genes
ipsgene <- c(c_gene, mcf7ips_overlap)
ipsgene <- unique(ipsgene)

##################
# READ IN PATIENT EXPRESSION DATA AND TIDY
##################

expressiondata <- data.frame(read.csv("~/Desktop/pnas_1518007112_sd01.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

expressiondata <- filter(expressiondata, gene_id != "?") #Getting rid of rows with "?" as a gene ID

expressiondata <- melt(expressiondata)

split <- str_split_fixed(expressiondata$variable, "_", 2)

expressiondata <- cbind.data.frame(split, expressiondata)

colnames(expressiondata) <- c("Patient", "Diagnosis", "Gene", "P_D", "Value")

#### FILTERING DATA WITH OUR CRITERIA OF INTEREST
expressiondata$Diagnosis <- ordered(expressiondata$Diagnosis, levels = c("BenignLo", "BenignHi", "CancerLo", "CancerHi"))
expressiondata <- arrange(expressiondata, Diagnosis)
#stage_diffs <- subset(patientDF, GENE %in% huvfips_overlap) #Subset rows that are matches to those in the targets
stage_diffs <- subset(expressiondata, Gene %in% totaltestgenes)
#stage_diffs <- subset(stage_diffs, GENE %in% surface_genes)



patient <- dcast(stage_diffs, P_D+Patient+Diagnosis ~ Gene, value.var = "Value") #Make sure there are no duplicate rows
patient <- arrange(patient, Diagnosis)
patient <- subset(patient, select = -c(Patient, Diagnosis))
patient[,1] <- as.character(patient[,1])

patient[,2:ncol(patient)] <- sapply(patient[,2:ncol(patient)],as.numeric)

sapply(patient, class)  #to check classes

patient <- patient %>% remove_rownames %>% column_to_rownames(var="P_D") # Making the patient IDs the rownames, not the first column

redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)  # Making a color palette for the heatmap 

#Setting up the data for a heatmap/cluster analysis
patient <- patient[, apply(patient, 2, sum)!=0] ### Required to remove all the columns with 0 in them to get a distance matrix
test <- as.matrix(patient)
test1 <- t(test)
heatmap.2(test, 
          #distfun = dist_cor, hclustfun = clus_wd2, 
          # clustering
          distfun = dist, 
          #hclust = clus,
          # scaling (genes are in rows)
          scale = "row",
          # color
          Rowv = FALSE,
          col = redblackgreen, 
          cexRow = 0.5,
          # labels
          #ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none")


heatmap(test1, Colv=NA, col=greenred(10),scale="none")
cor(t(test1))
t <- 1-cor(t(test1))
t <- as.dist(1-cor(t(test1)))
hc <- hclust(as.dist(1-cor(t(test1))))
plot(hc)
heatmap(test1, Rowv=as.dendrogram(hc) , Colv=NA, col=greenred(10), cexRow = 0.2)





h1 <- c(10,20,10,20,10,20,10,20)
h2 <- c(20,10,20,10,20,10,20,10)

l1 <- c(1,3,1,3,1,3,1,3)
l2 <- c(3,1,3,1,3,1,3,1)

mat <- rbind(h1,h2,l1,l2)

par(mar=c(4,4,1,1))
plot(1:8,rep(0,8), ylim=c(0,35), pch="", xlab="Time", ylab="Gene Expression")

for (i in 1:nrow(mat)) {
  lines(1:8,mat[i,], lwd=3, col=i)
}

legend(1,35,rownames(mat), 1:4, cex=0.7)

dist(mat)

hc <- hclust(dist(mat))
plot(hc)

heatmap(mat, Colv=NA, col=greenred(10))

heatmap(mat, Colv=NA, col=greenred(10), scale="none")

cor(t(mat))

1-cor(t(mat))

hc <- hclust(as.dist(1-cor(t(mat))))
plot(hc)

heatmap(mat, Rowv=as.dendrogram(hc), Colv=NA, col=greenred(10))



library(gplots)

# read the data in from URL
bots <- read.table(url("http://genome-www.stanford.edu/cellcycle/data/rawdata/combined.txt"), sep="t", header=TRUE)

# get just the alpha data
abot <- bots[,c(8:25)]
rownames(abot) <- bots[,1]
abot[1:7,]

# get rid of NAs
abot[is.na(abot)] <- 0

# we need to find a way of reducing the data.
# Sort on max difference and take first 1000
min <-apply(abot, 1, min)
max <- apply(abot, 1, max)
sabot <- abot[order(max – min, decreasing=TRUE),][1:1000,]

# cluster on correlation
hc <- hclust(as.dist(1 – cor(t(sabot))), method=”average”)

# draw a heatmap
heatmap(as.matrix(sabot),
        Rowv=as.dendrogram(hc),
        Colv=NA,
        col=greenred(10),
        labRow=””)



X = t(scale(t(test),center=TRUE,scale=FALSE))
sv = t(X)
w <- which(is.na(as.matrix(sv)))
#sv <- sv[which(!is.finite(as.matrix(sv)))] <- 0
sv[w] <- 0
sv = svd(sv)
U = sv$u
V = sv$v
D = sv$d

aa<- grep("grey",colors())
bb<- grep("green",colors())
cc<-  grep("red",colors())
gcol2<- colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

## use the genes that drive the first PC1. This is the first major patter in the data
k=1
ord1<- order(abs(V[,k]),decreasing=TRUE)
x1<- as.matrix(X[ord1[1:nrow(X)],])
heatmap(x1,col=gcol2)
heatmap.2(x1, col=gcol2)

pheatmap(x1,col = gcol2,
         # clustering
         #distfun = dist, 
         #hclust = clus,
         # scaling (genes are in rows)
         #scale = "row",
         # color
         #col = redblackgreen, 
         Rowv=FALSE, 
         # labels
         #labRow = "", 
         #ColSideColors = class_labels, 
         # tweaking
         #trace = "none",
         fontsize_row = 2,
         #density.info = "none")
)

# PVClust analysis using bootstrapping analysis:
library(pvclust)
result3 <- pvclust(test, nboot = 1000)
plot(result3, cex = 0.35, print.pv = FALSE)
pvrect(result3, alpha = 0.95)

##Site 2 for PCA 
#pca <- prcomp(t(gene_matrix))
pca <- prcomp(t(test)) # To see how genes may be clustered (genes become rows)
pca <- prcomp(test) #To see how Patients/Stage/Grades Cluster
#plot(pca$x[,1:2])
ipsgene <- c(c_gene, mcf7ips_overlap)
ipsgene <- unique(ipsgene)

#Adding back parameter of interest
pca$x <- rownames_to_column(as.data.frame(pca$x))
names(pca$x)[names(pca$x) == 'rowname'] <- 'PATIENT_ID'
pca$x <- merge(pca$x, extras, by="PATIENT_ID")

ggplot(as.data.frame(pca$x[,1:10])) + 
  geom_point(aes(x=PC1, y=PC2, color = pca$x$METASTASIS), size = 2) #+ 
#geom_text(aes(x=PC1, y=PC2, label = rownames(pca$x)), size = 4, hjust=-0.25, vjust=-0.5, color = ifelse(rownames(pca$x) %in% ipsgene, "red", "black"))

#Notes:

#  1.) In general, matrices of gene data are usually samples in columns
#and genes in rows, which is the transpose of what prcomp() expects, so you
#have to use t().

#2.) Usually when I plot the results, I also use pch, col, xlab, ylab,
#main, etc. to make the plotting symbols for each group different
#shapes and colors, add reasonable axis labels, a main title, etc. See ?par
#for other variables you can pass to plot.

#3.) It is nice to follow up with legend() to add a nice legend to the
#plot. People seem to like that stuff ;-D.

#4.) The general recommendation is to set scale. = TRUE in the call to
#prcomp. I'm not sure if it will make much difference with microarray
#data because it is usually normalized anyway, so I usually don't
#bother.
#However, it is worth a try.
