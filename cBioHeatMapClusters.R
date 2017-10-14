library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(ISLR)
library(gplots)
library("Biobase")

nan.na <- function (x) {
  x[is.nan(x)] <- NA
  return(x)
}

dist_cor <- function(x) {
  as.dist(1 - cor(t(x), method = "pearson"))
}

clus_wd2 <- function(x) {
  hclust(x, method = "ward.D2")
}

clus_wd <- function(x) {
  hclust(x, method = "ward.D")
}

file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/brca_tcga_pub/brca_tcga_pub/")
patientdata <- data.frame(read.csv(paste(file.path, "data_clinical.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
expressiondata <- read.csv(paste(file.path, "data_expression_median.csv", sep=""), sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

expressiondata[expressiondata == "null"] <- NA

colnames(patientdata) <- lapply(patientdata[5,], as.character)
patientdata <- patientdata[-(1:5),]
#patientdata <- patientdata[1:100,]

patientdata <- subset(patientdata, select = c(PATIENT_ID, Tumor, `AJCC Stage`)) 

expressiondata <- subset(expressiondata, select = -Entrez_Gene_Id)
expressiondata <- expressiondata %>% drop_na()
#x <- expressiondata[complete.cases(expressiondata), ]
expressionDF <- melt(expressiondata, id.vars = "Hugo_Symbol")

colnames(expressionDF) <- c("GENE", "PATIENT_ID", "EXPRESSION_LEVEL")

#Merging Data frames
patientDF <- merge(x = expressionDF, y = patientdata, by = "PATIENT_ID") # all = TRUE)

patientDF <- dcast(patientDF, PATIENT_ID+Tumor+`AJCC Stage`~ GENE, value.var = "EXPRESSION_LEVEL")

patientDF <- arrange(patientDF, `AJCC Stage`, Tumor)

patientDF <- patientDF %>% 
  unite(ID, PATIENT_ID, Tumor, 'AJCC Stage', sep = "_", remove = TRUE)

patientDF[,2:ncol(patientDF)] <- sapply(patientDF[,2:ncol(patientDF)],as.numeric)
sapply(patientDF, class)  #to check classes
patientDF <- patientDF %>% remove_rownames %>% column_to_rownames(var="ID")

#patientDF <- t(patientDF)




redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)

#class_labels <- ifelse(all_var$mol.biol == "NEG", "grey80", "grey20")

#x <- patientDF %>% drop_na()

#row.has.na <- apply(expressiondata, 1, function(x){any(is.na(x))})
#row.has.nan <- apply(expressiondata, 1, function(x){any(is.nan(x))})
#row.has.inf <- apply(expressiondata, 1, function(x){any(is.infinite(x))})
#row.has.null <- apply(expressiondata, 1, function(x){any(is.null(x))})
#sum(row.has.na)

#final.filtered <- final[!row.has.na,]



test <- as.matrix(patientDF)

heatmap.2(test, distfun = dist_cor, hclustfun = clus_wd2, 
          # clustering
          #distfun = dist, 
          #hclust = clus,
          # scaling (genes are in rows)
          scale = "row",
          # color
          col = redblackgreen, 
          # labels
          labRow = "", 
          #ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none")







set <- expressiondata[, ALL$mol.biol %in% c("BCR/ABL", "ALL1/AF4")]

heatmapdata = t(scale(t(expressiondata),center=TRUE,scale=FALSE))

X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
heatmapdata <- t(expressiondata)
heatmapdata <- t(scale(heatmapdata, center = TRUE, scale = TRUE))

colnames(expressionDF) <- c("GENE", "PATIENT_ID", "EXPRESSION_LEVEL")

#Merging Data frames
patientDF <- merge(x = expressionDF, y = patientdata, by = "PATIENT_ID") # all = TRUE)

patientDF <- arrange(patientDF, PATIENT_ID, `AJCC Stage`, Tumor, GENE)

######Hierarchial Clustering Analysis
zzz <- patientDF
zzz <- zzz[,-(1:2)]
colnames(zzz) <- NULL
#ncidat1 = t(NCI60$data)
ncidat <- t(zzz)
#colnames(ncidat1) = NCI60$labs
#unique(colnames(ncidat1))
#dim(ncidat)

#colnames(ncidat) <- lapply(ncidat[1,], as.character)
colnames(ncidat) <- expressiondata$Hugo_Symbol
#ncidat <- ncidat[-1,]  
#ncidat <- ncidat[,-1]  
X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
sv = t(X)
#is.na(X) <- do.call(cbind,lapply(X, is.infinite))
#length(which(!is.finite(as.matrix(X))))
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


expression$tumor_grade <- patientDF$GRADE[match(expression$variable, patientDF$PATIENT_ID)]
expression$stage <- patientDF$TUMOR_STAGE[match(expression$variable, patientDF$PATIENT_ID)]
#expression$value <- nan.na(ips.hits$value)

x <- subset(expression, select = -variable)

y <- dcast(x, Hugo_Symbol ~ stage, mean)
z <- dcast(x, Hugo_Symbol ~ stage, median)
#write.csv(y, "test.csv")

dy <- dist(as.matrix(y))   # find distance matrix 
dz <- dist(as.matrix(z))

hcy <- hclust(dy)                # apply hirarchical clustering 
hcz <- hclust(dz)

plot(hcy)




file.path <- ("~/Desktop/")
gene <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

split.names <- str_split_fixed(proteinDF$Composite.Element.REF,"[|]", 2) #special character, so needed [ ]
proteinDF <- data.frame(proteinDF, split.names)
proteinDF <- subset(proteinDF, select = -c(Composite.Element.REF, `X2`))
names(proteinDF)[names(proteinDF) == 'X1'] <- 'Gene'

ips.targets <- proteinDF[proteinDF$Gene %in% gene,]
ips.targets <- melt(ips.targets)
ips.targets <- subset(ips.targets, select = -variable)
ips.targets$value <- nan.na(ips.targets$value)
ips.targets <- ips.targets %>%
  group_by(Gene) %>%
  mutate(mean = mean(value, na.rm = TRUE))

ips.hits <- rnaDF[rnaDF$Hugo_Symbol %in% gene,]
ips.hits <- subset(ips.hits, select = -Entrez_Gene_Id)
ips.hits <- melt(ips.hits)
ips.hits <- subset(ips.hits, select = -variable)
ips.hits$value <- nan.na(ips.hits$value)

ggplot(ips.hits) + 
  geom_point(aes(x=Hugo_Symbol, y=value)) + 
  ggtitle("Median Z-Scores (mRNA) vs Gene by Tumor Grade") + 
  xlab("Median Z-score") + 
  ylab("Gene") + 
  theme(axis.text.x = element_text(angle = 90))

ips.hits <- ips.hits %>%
  group_by(Hugo_Symbol) %>%
  mutate(mean = mean(value, na.rm = TRUE))




######Hierarchial Clustering Analysis
install.packages("ISLR")
library(ISLR)
zzz <- expressiondata
zzz <- zzz[,-(1:2)]
colnames(zzz) <- NULL
#ncidat1 = t(NCI60$data)
ncidat <- t(zzz)
#colnames(ncidat1) = NCI60$labs
#unique(colnames(ncidat1))
#dim(ncidat)

#colnames(ncidat) <- lapply(ncidat[1,], as.character)
colnames(ncidat) <- expressiondata$Hugo_Symbol
#ncidat <- ncidat[-1,]  
#ncidat <- ncidat[,-1]  
X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
sv = t(X)
#is.na(X) <- do.call(cbind,lapply(X, is.infinite))
#length(which(!is.finite(as.matrix(X))))
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


expression$tumor_grade <- patientDF$GRADE[match(expression$variable, patientDF$PATIENT_ID)]
expression$stage <- patientDF$TUMOR_STAGE[match(expression$variable, patientDF$PATIENT_ID)]
#expression$value <- nan.na(ips.hits$value)

x <- subset(expression, select = -variable)

y <- dcast(x, Hugo_Symbol ~ stage, mean)
z <- dcast(x, Hugo_Symbol ~ stage, median)
#write.csv(y, "test.csv")

dy <- dist(as.matrix(y))   # find distance matrix 
dz <- dist(as.matrix(z))

hcy <- hclust(dy)                # apply hirarchical clustering 
hcz <- hclust(dz)

plot(hcy)

install.packages("factoextra")
install.packages("cluster")
install.packages("magrittr")
