# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)
library(biomaRt)
library(tidyverse)
library(ggpubr)
library(MCPcounter)
library(data.table)

# load series and platform data from GEO

#treatment information
#------------
#Neoadjuvant chemoradiation consisted of gemcitabine, erlotinib and radiation therapy in 5 patients and 
#gemcitabine and radiation therapy in one patient. Radiation was given as either 30 Gy or 36 Gy in 15 fractions

gset <- getGEO("GSE68172", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

# box-and-whisker plot
par(mar=c(7,4,2,1))
title <- paste ("GSE68172", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE68172", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE68172")

annotation <- fread("C:/data/AML/Microarray - GPL570/GPL570-probe-annotation.txt", sep="\t", header=TRUE)

count.matrix <- gset@assayData[["exprs"]]
count.matrix.df <- data.frame(ID=rownames(count.matrix), count.matrix)
count.matrix.annotated <- merge(annotation, count.matrix.df, by="ID")


#str_extract(count.matrix.annotated$gene_assignment,
            
tt <- str_match(count.matrix.annotated$gene_id,"// (.*?) //") 
count.matrix.annotated$gene_id <- tt[,2]
count.matrix.2 <- count.matrix.annotated[,c(11,17:93)]
count.matrix.2 <- data.frame(filter(count.matrix.2,!count.matrix.2=="" ))
rownames(count.matrix.2)<- make.unique(count.matrix.2[,1])


MetaData <- gset@phenoData@data
MetaData$Sorting <- paste(MetaData$characteristics_ch1.3, MetaData$characteristics_ch1.4)
#MetaData <- dplyr::rename(MetaData, Group = source_name_ch1)  
  
df.joined <- merge(MetaData, t(count.matrix.2[,-1]), by=0)

p <- ggplot(df.joined %>% filter(Sorting == "cd34 status: + plus cd38 status: - minus" & `patient id:ch1` != "patient 2"), aes(x=`disease status:ch1`, y=IL3RA, fill = `disease status:ch1`))+geom_boxplot()+theme_bw()+geom_point()+
  theme(panel.spacing = unit(0.75, "lines"),axis.line = element_line(colour = "black"),text = element_text(size=20),
        axis.text.x = element_text(size=20), axis.text.y = element_text(size=24),  legend.text=element_text(size=20))
#p+stat_compare_means(ref.group = "untreated")
p+stat_compare_means()

write.table(MetaData,"C:/data/AML/Microarray - GPL570/GSE68172-annotation.txt", sep="\t", row.names=FALSE)
write.table(count.matrix.2,"C:/data/AML/Microarray - GPL570/GSE68172-normalized.counts.txt", sep="\t", row.names=FALSE)

#------------------------
#run MCP counter
#--------------------
mcp.output <- MCPcounter.estimate(count.matrix.2[,-1],featuresType = "HUGO_symbols",genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE))
df.mcp.output <- data.frame(SAMPLE_ID=colnames(mcp.output), t(mcp.output))
df.mcp <- merge(df.mcp.output)


