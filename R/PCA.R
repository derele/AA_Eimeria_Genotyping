library(MultiAmplicon)
library(pheatmap)
library(data.table)
library(DECIPHER)
library(ape)
library(pegas)
library(sidier)
library(adegenet)
library(wordcloud)
library(systemPipeR)
library(Biostrings)
library(ade4)
library(hierfstat)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(poppr)


Nuclear35<- TRUE
Apicoplast5 <- FALSE


if(Nuclear35){
  oneAln <- readDNAStringSet(file="~/GitProjects/AA_Eimeria_Genotyping/Alignments/combinedEfalNucEver_all.fasta", format="fasta")
}

if(Apicoplast5){
  oneAln <- readDNAStringSet(file="~/GitProjects/AA_Eimeria_Genotyping/Alignments/combinedEfalApiEver_all.fasta", format="fasta")
}


## create a Network representation despite missing data

## dmat <- as.matrix(oneAln)
dbin <- as.DNAbin(oneAln)
gind <- DNAbin2genind(dbin) ##Conserves SNPs only 

gind@type <- "PA" ##Presence/Absence

gind@pop <- as.factor(ifelse(rownames(gind@tab)%in%c("SK_2808", "SK_2809", "E_vermiformis", "E_falciformis", 
                                                     "AA_0079_CE", "AA_0100_CE", "SK_3019","AA_0366", "AA_0204_IL",
                                                     "AA_0080_IL", "AA_0054_IL", "AA_0111_IL", "AA_0112_CE","AA_0140_CE"),
                             "Mus", "Apodemus"))

#barplot(table(pop(gind)), col=funky(17), las=3,
#      xlab="Population", ylab="Sample size")

#Y <-summary(gind)

##Assesing population structure
##Fst for population differentation


X <- tab(gind, NA.method="zero")

pca1 <- dudi.pca(X, scannf=FALSE, scale=FALSE)
temp <- as.integer(pop(gind))
myCol <- transp(c("blue","red"),.7)[temp]
myPch <- c(15,17)[temp]


if(Nuclear35){
  pdf("~/GitProjects/AA_Eimeria_Genotyping/Figures/dudiPCA_Nuclear.pdf",  width=15, height=15, onefile=FALSE)
  plot(pca1$li, col=myCol, cex=3, pch=myPch)
  textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=1.4, new=FALSE)
  abline(h=0,v=0,col="grey",lty=2)
  ## here it makes more senst to add up the contribution of different markers
  ## s.arrow(pca1$c1*.5, add.plot=TRUE)
  legend("topright", pch=c(15,17), col=transp(c("blue","red"),.7),
         leg=c("Apodemus","Mus"), pt.cex=2)
  u <- par("usr")
  v <- c(grconvertX(u[1:2], "user", "ndc"),
         grconvertY(u[3:4], "user", "ndc"))
  v <- c( (v[1]+v[2])/2, v[2], (v[3]+v[4])/2, v[4] )
  par(fig=v, new=TRUE, mar=c(0,0,0,0))
  barplot(pca1$eig[1:60],main="PCA eigenvalues", col=heat.colors(50))
  dev.off()
}

if(Apicoplast5){
  pdf("~/GitProjects/AA_Eimeria_Genotyping/Figures/dudiPCA_Api.pdf",  width=15, height=15, onefile=FALSE)
  plot(pca1$li, col=myCol, cex=3, pch=myPch)
  textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=1.4, new=FALSE)
  abline(h=0,v=0,col="grey",lty=2)
  ## here it makes more senst to add up the contribution of different markers
  ## s.arrow(pca1$c1*.5, add.plot=TRUE)
  legend("topright", pch=c(15,17), col=transp(c("blue","red"),.7),
         leg=c("Apodemus","Mus"), pt.cex=2)
  u <- par("usr")
  v <- c(grconvertX(u[1:2], "user", "ndc"),
         grconvertY(u[3:4], "user", "ndc"))
  v <- c( (v[1]+v[2])/2, v[2], (v[3]+v[4])/2, v[4] )
  par(fig=v, new=TRUE, mar=c(0,0,0,0))
  barplot(pca1$eig[1:60],main="PCA eigenvalues", col=heat.colors(50))
  dev.off()
}

###PCA scaled
pca2<- dudi.pca(X, scannf=F, scale=T)

s.class(pca2$li, fac = pop(gind), col = transp(c("#03BA0F","#1703F9"), .6))
add.scatter.eig(pca2$eig,NULL,1,2, posi = "bottomright", ratio = .3)

s.class(pca2$li, fac = pop(gind), 
        col = transp(c("#03BA0F","#1703F9"), .6),
        axesell = T, cstar = 0, cpoint = 3)
add.scatter.eig(pca1$eig,NULL,1,2, posi = "bottomleft", ratio = .3)
loadingplot(pca2$c1^2, threshold = 0.0008)


##FactoMineR
res.pca <- PCA(X, scale.unit = T, ncp= 10, graph = F)
eigenvalues <- res.pca$eig
head(eigenvalues[,1:2])
plot(res.pca,  ellipse.level=0.95)


##Squared cosine (cos2)
if(Nuclear){
  pdf("~/GitProjects/AA_Eimeria_Genotyping/Figures/Squared_cosine_Nuc.pdf",  width=15, height=15, onefile=FALSE)
  fviz_pca_ind(res.pca,col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = T)
  dev.off()
}

if(Apicoplast){
  pdf("~/GitProjects/AA_Eimeria_Genotyping/Figures/Squared_cosine_Api.pdf",  width=15, height=15, onefile=FALSE)
  fviz_pca_ind(res.pca,col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = T)
  dev.off()
}

