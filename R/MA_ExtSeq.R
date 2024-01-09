library(MultiAmplicon)
library(pheatmap)
library(data.table)
library(DECIPHER)
library(ape)
library(pegas)
library(sidier)
library(adegenet)
#library(wordcloud)
#library(systemPipeR)

##Change these parameter according to the steps to re-run and dataset to analyse 

## ## We don't re-run any analyses but use the dat as published in
## "Generalist Eimeria species in rodents: Multilocus analyses
## indicate inadequate resolution of established markers" (2020)
## Víctor Hugo Jarquín-Díaz, Alice Balard, Anna Mácová, Jenny Jost,
## Tabea Roth von Szepesbéla, Karin Berktold, Steffen Tank, Jana
## Kvičerová, Emanuel Heitlinger https://doi.org/10.1002/ece3.5992

Filtering <- FALSE
newMA <- FALSE
newDeDa <- FALSE

Mus <- TRUE ##Database from Mus musculus
Apo <- FALSE ##Database from Apodemus 
Nuc <- FALSE ##Nuclear genes analysis
Api <- FALSE ##Apicoplast regions analysis


if (Mus){files <- list.files(path="/SAN/Victors_playground/mouseGT/raw_fastq/Mus/",
                             pattern=".fastq.gz", full.names=TRUE)}


if (Apo){files <- list.files(path="/SAN/Victors_playground/mouseGT/raw_fastq/Apodemus/",
                             pattern=".fastq.gz", full.names=TRUE)}


if(Mus){
  ## first pool file names
  Ffq.file <- files[grepl("R1_001", files)]
  Rfq.file <- files[grepl("R2_001", files)]
  
  ## first pool sample names
  samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(Ffq.file))
  samples <- gsub("S\\d+-", "\\1", basename(samples))
  samples <- gsub("-", "_", basename(samples))
  
  ## plotQualityProfile(Ffq.file[[14]])
  
  
  filt_path <- "/SAN/Victors_playground/mouseGT/filt/Mus/"
  if(!file_test("-d", filt_path)) dir.create(filt_path)
  
  filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
  names(filtFs) <- samples
  filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
  names(filtRs) <- samples
  
  if (Filtering){
    filter.track <- lapply(seq_along(Ffq.file),  function (i) {
      filterAndTrim(Ffq.file[i], filtFs[i], Rfq.file[i], filtRs[i],
                    truncLen=c(245,245), minLen=c(245,245), 
                    maxN=0, maxEE=4, truncQ=4, 
                    compress=TRUE, verbose=TRUE)
    })
    saveRDS(filter.track, file="/SAN/Victors_playground/mouseGT/filt/Mus/filter.Rds")
  } else {
    filter.track <- readRDS(file="/SAN/Victors_playground/mouseGT/filt/Mus/filter.Rds")
  }

  filter <- do.call(rbind, filter.track)
  colSums(filter)[2]/colSums(filter)[1]
  ## low proportion of data overall cept (but okay for this purpose!)

##
  file.ptable <- "data/primer.file.csv"
  
  ptable <- read.csv(file.ptable, sep=",", header=TRUE)
  
  primerF <- gsub(" ", "", ptable[,2])
  primerR <- gsub(" ", "", ptable[,4])
  
  names(primerF) <- as.character(ptable[,1])
  names(primerR) <- as.character(ptable[,3])
  
  files <- PairedReadFileSet(filtFs, filtRs)
  
  primers <- PrimerPairsSet(primerF, primerR)
  
  rownames(ptable) <- names(primers)
  
  MA <- MultiAmplicon(primers, files)
  
  
  stratfiles <- "/SAN/Victors_playground/mouseGT/Stratified_files/Mus"
  
  if(newMA) {
    if(dir.exists(stratfiles)){
      unlink(stratfiles,
             recursive=TRUE)
    }
    MA1 <- sortAmplicons(MA, filedir=stratfiles,
                         starting.at=1, max.mismatch=4)
    saveRDS(MA1, file="/SAN/Victors_playground/mouseGT/MA_files/Mus/MA1.Rds")
  } else {
    if(!newMA){
      MA1 <- readRDS(file="/SAN/Victors_playground/mouseGT/MA_files/Mus/MA1.Rds")
    } else {stop("Want new sorting or not? Set newMA to TRUE")}
  } 
}

if(Apo){
  ## second pool file names
  Ffq.file <- files[grepl("R1_001", files)]
  Rfq.file <- files[grepl("R2_001", files)]
  
  samples <- gsub("/SAN/Victors_playground/mouseGT/raw_fastq/Apodemus/(.*?)_S\\d{,3}_L001_R1_001.fastq\\.gz", "\\1", Ffq.file)
  samples <- gsub("/", "", basename(samples))
  
  filt_path <- "/SAN/Victors_playground/mouseGT/filt/Apodemus/"
  if(!file_test("-d", filt_path)) dir.create(filt_path)
  
  filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
  names(filtFs) <- samples
  filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
  names(filtRs) <- samples
  
  if (Filtering){
    filter.track <- lapply(seq_along(Ffq.file),  function (i) {
      filterAndTrim(Ffq.file[i], filtFs[i], Rfq.file[i], filtRs[i],
                    truncLen=c(285,285), minLen=c(285,285), 
                    maxN=0, maxEE=2, truncQ=2, 
                    compress=TRUE, verbose=TRUE)
    })
    saveRDS(filter.track, file="/SAN/Victors_playground/mouseGT/filt/Apodemus/filter.Rds")
  } else {
    filter.track <- readRDS(file="/SAN/Victors_playground/mouseGT/filt/Apodemus/filter.Rds")
  }
  
  filter <- do.call(rbind, filter.track)
  colSums(filter)[2]/colSums(filter)[1]
  ## low proportion of data overall (but okay for this purpose!)
  
  file.ptable <- "/SAN/Victors_playground/mouseGT/primer.file.csv"
  
  ptable <- read.csv(file.ptable, sep=",", header=TRUE)
  
  primerF <- gsub(" ", "", ptable[,2])
  primerR <- gsub(" ", "", ptable[,4])
  
  names(primerF) <- as.character(ptable[,1])
  names(primerR) <- as.character(ptable[,3])
  
  files <- PairedReadFileSet(filtFs, filtRs)
  
  primers <- PrimerPairsSet(primerF, primerR)
  
  rownames(ptable) <- names(primers)
  
  MA <- MultiAmplicon(primers, files)
  
  stratfiles <- "/SAN/Victors_playground/mouseGT/Stratified_files/Apodemus/"
  
  if(newMA) {
    if(dir.exists(stratfiles)){
      unlink(stratfiles,
             recursive=TRUE)
    }
    MA1 <- sortAmplicons(MA, filedir=stratfiles,
                         starting.at=1, max.mismatch=4)
    saveRDS(MA1, file="/SAN/Victors_playground/mouseGT/MA_files/Apodemus/MA1.Rds")
    
    pdf("/SAN/Victors_playground/mouseGT/Figures/primers_MA_sorted_Apo.pdf", width=46)
    plotAmpliconNumbers(MA1)
    dev.off()
  } else {
    if(!newMA){
      MA1 <- readRDS(file="/SAN/Victors_playground/mouseGT/MA_files/Apodemus/MA1.Rds")
    } else {stop("Want new sorting or not? Set newMA to TRUE")}
  } 
}

###Reduce the data set just for sequences with high number of amplicons generated

if(Mus){
MA1 <- readRDS(file="/SAN/Victors_playground/mouseGT/MA_files/Mus/MA1.Rds")
clust <- plotAmpliconNumbers(MA1)
two.clusters.col <- cutree(clust$tree_col, k=2)
keep.sample <- names(two.clusters.col)[two.clusters.col==1]
MA1 <- MA1[, which(colnames(MA1)%in%keep.sample)]
}

if(Apo){
MA1 <- readRDS(file="/SAN/Victors_playground/mouseGT/MA_files/Apodemus/MA1.Rds")
clust <- plotAmpliconNumbers(MA1)
two.clusters.col <- cutree(clust$tree_col, k=2)
keep.sample <- names(two.clusters.col)[two.clusters.col==2]
MA1 <- MA1[, which(colnames(MA1)%in%keep.sample)] 
}

####Separate Nuclear from Apicoplast sequences 

if(Nuc){
MA1 <- MA1[which(grepl("^Eim_Nuc", rownames(MA1))), ] ##Nuclear
}

if(Api){
MA1 <- MA1[which(grepl("^Eim_Ap", rownames(MA1))), ] ##Apicoplast 
}

if(newDeDa){
  MA2 <- derepMulti(MA1, mc.cores=20)
  MA3 <- dadaMulti(MA2, Ferr=NULL, Rerr=NULL, selfConsist=TRUE,
                   multithread=FALSE, mc.cores=20, verbose=0, MAX_CONSIST=20)
  MA4 <- mergeMulti(MA3, mc.cores=20, minOverlap= 10)
  MA5 <- makeSequenceTableMulti(MA4, mc.cores=20, orderBy="nsamples")
  MA6 <- removeChimeraMulti(MA5)
  if(Nuc){ 
    if(Mus){ 
  saveRDS(MA6, file="/SAN/Victors_playground/mouseGT/MA_files/Mus/MA6_Nuc.Rds")}
    if(Apo){
  saveRDS(MA6, file="/SAN/Victors_playground/mouseGT/MA_files/Apodemus/MA6_Nuc.Rds")}}
  if(Api){ 
    if(Mus){ 
      saveRDS(MA6, file="/SAN/Victors_playground/mouseGT/MA_files/Mus/MA6_Api.Rds")}
  if(Apo){
    saveRDS(MA6, file="/SAN/Victors_playground/mouseGT/MA_files/Apodemus/MA6_Api.Rds")}
} else {
  if(Nuc){ 
    if(Mus){ 
  MA6 <- readRDS(file="/SAN/Victors_playground/mouseGT/MA_files/Mus/MA6_Nuc.Rds")}
    if(Apo){
  MA6 <- readRDS(file="/SAN/Victors_playground/mouseGT/MA_files/Apodemus/MA6_Nuc.Rds")}
  }
  if(Api){ 
    if(Mus){ 
      MA6 <- readRDS(file="/SAN/Victors_playground/mouseGT/MA_files/Mus/MA6_Api.Rds")}
    if(Apo){
      MA6 <- readRDS(file="/SAN/Victors_playground/mouseGT/MA_files/Apodemus/MA6_Api.Rds")}}
  }
}

## plotAmpliconNumbers(MA6)

### Remove chimeras 
STnoC <- getSequenceTableNoChime(MA6) 

#hist(log10(unlist(STnoC)), breaks = 100) 

##Discard empty samples and ASVs for each amplicon 
STnoC <- lapply(STnoC, function(x) x[rowSums(x)>0, colSums(x)>0])

if(Mus){
  if(Nuc){
    STnoC$Eim_Nuc9409_F.Eim_Nuc9410_R <- NULL ##This amplicon is empty and bug
    STnoC$Eim_Nuc3331_F.Eim_Nuc3332_R <- NULL ##This amplicon is not empty but bug 
    STnoC <- lapply(STnoC, function(x) x[,colSums(x)>400]) 
    STnoC <- lapply(STnoC, function(x) x[rowSums(x)>0,]) #just for Mus Nuclear to avoid problems in the next function
  }
}

STnoC <- STnoC[unlist(lapply(STnoC, function (x) all(dim(x) >1)))]


pdf("/SAN/Victors_playground/mouseGT/Figures/all_amps_heated.pdf") ## Supplementary file 
lapply(STnoC, function (x) {
  pheatmap(x, show_colnames=FALSE, scale="none")
})
dev.off()

## get the most represented sequence for each amplicon and sample
DNA <- lapply(STnoC, function (amp){
  mostAbSeq <- apply(amp, 1, function (x){
    mostAbundant <- which(x==max(x))
    names(x)[mostAbundant]
  })
  DNAStringSet(mostAbSeq)
})

## now we can remove the Ns and just let the alignment put gaps
alnDNA1 <- lapply(DNA, function (x) AlignSeqs(x))

## now exclude too distant sequences (few seqs are from wrong targets e.i. Bacterial DNA)

DNA2 <- lapply(alnDNA1, function(x){
  distDNA <- DistanceMatrix(x, penalizeGapLetterMatches = FALSE)
  d <- colSums(distDNA>0.2, na.rm=TRUE) 
  ## so we can exclude sequences which have with more than the are in
  ## the minority as having huge differences 
  toKeep <- which (d < (length(x)/2)) 
  x[toKeep]
})

if(Api){
  if(Apo){
    DNA2$Eim_Ap9_F.Eim_Ap9_R <- NULL ##Empty
  }
}

## realign removeing the gaps before
alnDNA2 <- lapply(DNA2, function (x) AlignSeqs(DNAStringSet(gsub("-", "", x))))


## checking that huge distances are gone
distDNA2 <- lapply(alnDNA2, DistanceMatrix, penalizeGapLetterMatches = FALSE)

## to see maximal distances
summary(unlist(lapply(distDNA2, max, na.rm=TRUE)))

## to check alignment length vs longest sequence
longSeq <- unlist(lapply(DNA, function (x) max(width(x))))
alnSeq1 <- unlist(lapply(alnDNA1, function (x) max(width(x))))
alnSeq2 <- unlist(lapply(alnDNA2, function (x) max(width(x))))
cbind(longSeq, alnSeq1, alnSeq2)


## remove the worst markers producing hardly any sequence for lots of
## samples

## remove duplicate technical replicate sequences for all samples
## return NNN sequences to have all alignments have the same sequences

if(Mus){
all.samples <- as.factor(unique(unlist(lapply(alnDNA2, names))))
}

if(Apo){
all.reps <- unique(unlist(lapply(alnDNA2, names)))
all.samples <- as.factor(gsub("-chip\\d", "", all.reps))
}


Aln <- lapply(alnDNA2, function (aln){
  alnWidth <- unique(width(aln))
  sampleNames <- gsub("-chip\\d", "", names(aln))
  ## could probably also work with an empty DNAstring("")
  badString <- DNAString(paste(rep("N", times=alnWidth), collapse=""))
  marker <- sapply(levels(all.samples), function (x){
    pair <- which(sampleNames%in%x)
    DNA <- aln[pair]
    ## ## if it's only one sequence return it
    if(length(DNA)<1){
      return(badString)
    }
    if(all(!grepl("[A|G|C|T]", DNA))){return(badString)}
    if(length(DNA)<2){return(DNA[[1]])}
    ## ## if both are gap only return one of them (the first)        
    ## PROBLEMATIC: What if the two technical reps don't agree?
    ## repD <- DistanceMatrix(DNA)
    ## if(any(repD>0)){
    ##     warning(paste("replicate dist >0:", repD[repD>0][1]))
    ##     return(DNAString(paste(rep("N", times=alnWidth),
    ##                            collapse="")))
    ## }
    else {
      return(DNA[[1]])
    }
  })
  DNAStringSet(marker) 
})

## throw out markers with less than optimal number of sequences
#table(unlist(lapply(Aln, function (x) length(grep("[A|G|C|T]", x)))))

#Aln <- Aln[unlist(lapply(Aln, function (x) length(grep("[A|G|C|T]", x))))>13]

## to get the alingment partitions
cumsum(unlist(lapply(Aln, function (x) max(width(x)))))

## paste the string sets togehter seq by seq
pasteDNAstrings <- function (x, y){
  DNAStringSet(paste0(as.character(x), as.character(y)))
}

oneAln <- Reduce(pasteDNAstrings, Aln)
names(oneAln) <- names(Aln[[1]])

## Frame analysis:
inFrameAln <- lapply(Aln, function (x) {
  DNA <- list()
  ## remove gaps and prefix G's which don't produce stop codons to
  ## see what brings the alignment in frame
  DNA[[1]] <- DNAStringSet(gsub("-|N", "", x))
  DNA[[2]] <- DNAStringSet(paste0("G", DNA[[1]]))
  DNA[[3]] <- DNAStringSet(paste0("GG", DNA[[1]]))
  DNA[[4]] <- reverseComplement(DNAStringSet(gsub("-|N", "", x)))
  DNA[[5]] <- reverseComplement(DNAStringSet(paste0(DNA[[1]], "C")))
  DNA[[6]] <- reverseComplement(DNAStringSet(paste0(DNA[[1]], "CC")))
  options(warn=-1)
  ## max because we hope that one sequence at least hd no gaps at
  ## all and we keep everthing at this length adding gaps back in
  ## when aligning at the end
  none <- max(nchar(gsub("\\*.*", "",
                         translate(DNA[[1]][width(DNA[[1]])>0]))))
  one <- max(nchar(gsub("\\*.*", "",
                        translate(DNA[[2]][width(DNA[[2]])>1]))))
  two <- max(nchar(gsub("\\*.*", "",
                        translate(DNA[[3]][width(DNA[[3]])>2]))))
  ## seems to be unnecessary, as all on + strand
  noneR <- max(nchar(gsub("\\*.*", "",
                          translate(DNA[[4]][width(DNA[[4]])>0]))))
  oneR <- max(nchar(gsub("\\*.*", "",
                         translate(DNA[[5]][width(DNA[[5]])>1]))))
  twoR <- max(nchar(gsub("\\*.*", "",
                         translate(DNA[[6]][width(DNA[[6]])>2]))))
  options(warn=1)
  ## hoping that one sequence has no gaps at all
  orf.to <- c(none, one, two, noneR, oneR, twoR)
  ## dodgy: having to avoid that in some cases (amplicon 30) tow
  ## sequences have maximal orf length
  this <- which(orf.to==max(orf.to))[[1]]
  this.cor <- (this+2)%%3
  cat("THIS:", this, "\n")
  if(this.cor>0){
    rm.pattern <- paste0(rep("G", times=this.cor), collapse="")
    rm.pattern <- paste0("^", rm.pattern)
    DNA[[this]] <- DNAStringSet(gsub(rm.pattern, "",  DNA[[this]]))
  }
  list(AlignSeqs(DNA[[this]]), orf.to[[this]]*3, this.cor)
})

## yeah baby, now correct
lapply(inFrameAln, function (x) cbind(unique(width(x[[1]])), x[[2]]))

frameAln <- lapply(inFrameAln, "[[", 1)

toORF <- lapply(inFrameAln, "[[", 2)
toFrame <- unlist(lapply(inFrameAln, "[[", 3))+1

oneFrameAln <- Reduce(pasteDNAstrings, frameAln)
names(oneFrameAln) <- names(Aln[[1]])


partition.to <- cumsum(unlist(lapply(frameAln, function (x) max(width(x)))))


if(Nuc){
  if(Mus){
    writeXStringSet(oneFrameAln,
                    file= "/SAN/Victors_playground/mouseGT/Alignments/combinedEfalNucFrame_Mus.fasta")
    
    write.csv(cbind(partition.to, toFrame),
              file="/SAN/Victors_playground/mouseGT/Alignments/combinedEfalNucFramePartitions_Mus.csv")
  }
  if(Apo){
    writeXStringSet(oneAln,
                    file="/SAN/Victors_playground/mouseGT/Alignments/combinedEfalNucFrame_Apo.fasta")
    write.csv(cbind(partition.to, toFrame),
              file="/SAN/Victors_playground/mouseGT/Alignments/combinedEfalNucFramePartitions_Apo.csv")
    }
}

if(Api){
  if(Mus){
    writeXStringSet(oneFrameAln,
                    file= "/SAN/Victors_playground/mouseGT/Alignments/combinedEfalApiFrame_Mus.fasta")
    
    write.csv(cbind(partition.to, toFrame),
              file="/SAN/Victors_playground/mouseGT/Alignments/combinedEfalApiFramePartitions_Mus.csv")
  }
  if(Apo){
    writeXStringSet(oneFrameAln,
                    file="/SAN/Victors_playground/mouseGT/Alignments/combinedEfalApiFrame_Apo.fasta")
    write.csv(cbind(partition.to, toFrame),
              file="/SAN/Victors_playground/mouseGT/Alignments/combinedEfalApiFramePartitions_Apo.csv")
  }
}
