#load library
library(adegenet)
library(tidyverse)
library(ggplot2)
###set your working directory in your local computer to the folder where you genotype structure file is
setwd("")


#read the structure file an save it as a genind object called data
data <- read.structure("sal_rel.stru", n.ind=84, n.loc=3603, onerowperind=FALSE, col.lab=1, col.pop=2, row.marknames=2, NA.char=0, ask=FALSE)

#names for the different populations included in the structure file
levels(data@pop) <- c("M_seg_KOL", "M_sal_SR5", 'M_sal_KER', "M_sal_BBay", "M_sal_PAA", "M_sal_IRE", "M_sal_MJO", "M_rel_TOP", "M_rel_TVA", "M_rel_BBay")

#choose one colour for each population
col <- c("darkgreen","darkgreen","darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "blue","blue", "blue" )
sum(is.na(data$tab))

#fills missing data with average
X <- scaleGen(data, NA.method="mean")

#performs PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)

# estimates contribuion of the first PC's to print them in the final plot
PC1 <- (100*pca1$eig/sum(pca1$eig))[1]
PC2 <- (100*pca1$eig/sum(pca1$eig))[2]
PC3 <- (100*pca1$eig/sum(pca1$eig))[3]

#plots PCA
s.class(pca1$li, pop(data),xax=1,yax=2, col=(col),addaxes=TRUE, axesell=TRUE, cstar=1, cpoint=5, grid=FALSE, sub=paste("PC1=", round(PC1, digits=2), "%, PC2=", round(PC2, digits=2), "%", sep=""), csub=1.5)

#Exact results from dudi.pca(df = x.apo, center = TRUE, scale = FALSE, scannf = FALSE, nf = 3)
#Calculate % of ech component
eig.perc <- 100*pca.apo$eig/sum(pca.apo$eig)
eig.perc
head(eig.perc)

#plot pca
p<-ggplot(pca1$li,aes(x=Axis1,y=Axis2,color=data$pop)) 
p<-p+scale_color_manual(labels = c("M_seg_KOL","M_sal_SR5", 'M_sal_KER',"M_sal_BBay", "M_sal_PAA", "M_sal_IRE", "M_sal_MJO", "M_rel_TOP", "M_rel_TVA", "M_rel_BBay"), values = c("yellowgreen", "yellowgreen","yellowgreen", "yellowgreen", "yellowgreen", "yellowgreen", "yellowgreen", "deepskyblue","deepskyblue", "deepskyblue"))
p<-p + geom_point(aes(shape=data$pop), size=5) + scale_shape_manual(values =c(0,1,2,15,3,17,6,19,8,5,1,2)) 
p<-p + coord_equal() +  theme_classic()
p + xlab("PCA1=23.4")+ylab("PCA2=2.0")


