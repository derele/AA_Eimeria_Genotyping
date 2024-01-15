## Rewriting this completely for the host genotyping
library(MultiAmplicon)
library(dada2)
library(phyloseq)
library(parallel)
library(pheatmap)

## Here, for complete reproducibility, we start from scratch
## downloading SRA data

SRA <- read.csv("data/SraRunTable.txt")
rownames(SRA) <- SRA$Run

## change as required on your system!
dataDir <- "/SAN/MouseAAGT2024/raw"

## Should a new download be performed
DoDownload <- FALSE
## using NCBI toolkit to download (see https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/)
commands <- paste("/usr/local/ncbi/sra-tools/bin/fastq-dump --split-3 --outdir",
                  dataDir, SRA$Run)

if(DoDownload){
    mclapply(commands, system, mc.cores=24)
} else {
    message("using raw data: ",
            length(list.files(dataDir)), " files in ", dataDir)
}

## 
Filtering <- FALSE
newMA <- TRUE
newDeDa <- TRUE

files <- list.files(path=dataDir,
                    pattern=".fastq", full.names=TRUE)

## first pool file names
Ffq.file <- files[grepl("1\\.fastq", files)]
Rfq.file <- files[grepl("2\\.fastq", files)]
  
## first pool sample names
samples <- gsub("_1\\.fastq", "", basename(Ffq.file))
  
filt_path <- "/SAN/MouseAAGT2024/filtered"
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
    saveRDS(filter.track, file="data/filter.Rds")
} else {
    filter.track <- readRDS(file="data/filter.Rds")
}

filter <- do.call(rbind, filter.track)
colSums(filter)[2]/colSums(filter)[1]
## low proportion of data overall cept (but okay for this purpose!)

PF <- read.csv("data/Primers_F.csv")
PR <- read.csv("data/Primers_R.csv")

## the same
## altPrimers <- read.delim("data/primer_table.csv", sep=";", header=FALSE)

ptable <- cbind(PF, PR)
  
primerF <- gsub(" ", "", ptable[,2])
primerR <- gsub(" ", "", ptable[,4])
  
names(primerF) <- as.character(ptable[,1])
names(primerR) <- as.character(ptable[,3])
  
filtfiles <- PairedReadFileSet(filtFs, filtRs)
  
primers <- PrimerPairsSet(primerF, primerR)
  
rownames(ptable) <- names(primers)
  
MA <- MultiAmplicon(primers, filtfiles)
  
  
stratfiles <- "/SAN/MouseAAGT2024/Stratified_files/"
  
if(newMA) {
    if(dir.exists(stratfiles)){
        stop("stratified files", stratfiles, "exists")
    }
    MA1 <- sortAmplicons(MA, filedir=stratfiles, starting.at=1, max.mismatch=4)
    saveRDS(MA1, file="data/MA1GT.Rds")
  } else {
    if(!newMA){
      MA1 <- readRDS(file="data/MA1GT.Rds")
    } else {stop("Want new sorting or not? Set newMA to TRUE")}
  } 

## Not reduce the data set just for sequences with high number of
## amplicons generated... just look at the statified data

clust <- plotAmpliconNumbers(MA1)
## two.clusters.col <- cutree(clust$tree_col, k=2)
## keep.sample <- names(two.clusters.col)[two.clusters.col==1]
## MA1 <- MA1[, which(colnames(MA1)%in%keep.sample)]

if(newDeDa){
  MA2 <- derepMulti(MA1, mc.cores=20)
  MA3 <- dadaMulti(MA2, Ferr=NULL, Rerr=NULL, selfConsist=TRUE,
                   multithread=FALSE, mc.cores=20, verbose=0, MAX_CONSIST=20)
  MA4 <- mergeMulti(MA3, mc.cores=20, minOverlap= 10)
  MA5 <- makeSequenceTableMulti(MA4, mc.cores=20, orderBy="nsamples")
  MA6 <- removeChimeraMulti(MA5)
  saveRDS(MA6, file="data/MA6.Rds")
} else {
      MA6 <- readRDS(file="data/MA6.Rds")}
}

pdf("Figures/HeatAllamp.pdf", width=30, height=7)
plotAmpliconNumbers(MA6)
dev.off()

## Adding SOTA data to SRA
SOTA <- read.csv("../Mouse_Eimeria_Field/data_products/SOTA_Data_Product.csv")
SOTA$Host <-  NULL ## strange strang strange column in SOTA
SRA$SK_ID <- gsub("-chip\\d", "", SRA$Sample_ID)
Sdat <- merge(SRA, SOTA, by.x="SK_ID", by.y="Mouse_ID")
rownames(Sdat) <- Sdat$Run

## adding SOTA also reducing to only house mouse samples the automatic
## sample data addig BUGGs here failing to remove all aspects of the
## non-matching samples (e.g. from rawCounts, but likely others as
## well).
MA7 <- MA6[, which(colnames(MA6)%in%Sdat$Run)]

MA8 <- addSampleData(MA7, Sdat)
table(getSampleData(MA8)$Host%in%"Mus musculus")
### All 90 TRUE

pdf("Figures/HeatAllampMus.pdf", width=15, height=7)
plotAmpliconNumbers(MA8)
dev.off()

MA9 <- MA8[which(grepl("Mus_", rownames(MA8))), ]

## remove non-amplified amplcions and samples
colSums(getRawCounts(MA9))[order(colSums(getRawCounts(MA9)))]
## under 1000 reads are bad for a sample

rowSums(getRawCounts(MA9))[order(rowSums(getRawCounts(MA9)))]
## under 1000 reads are bad for a primer pair

MA9 <- MA9[which(rowSums(getRawCounts(MA9))>1000),
           which(colSums(getRawCounts(MA9))>1000)]

RC9 <- getRawCounts(MA9)

annDat <- as.data.frame(unclass(getSampleData(MA9)))
rownames(annDat) <- annDat$Run
annDat <- annDat[, c("HI", "geo_loc_name")]

pdf("Figures/HeatMusampMus.pdf", width=15, height=7)
plotAmpliconNumbers(MA9, annotation_col=annDat)
dev.off()

lapply(getSequenceTable(MA9), function(x) nrow(x))
lapply(getSequenceTable(MA9), function(x) sum(x))
### Why $Mus_mtDNA0_F.Mus_mtDNA0_R are all zero?!?

lapply(getSequenceTableNoChime(MA9), function(x) nrow(x))
lapply(getSequenceTableNoChime(MA9), function(x) sum(x))
### $Mus_mtDNA0_F.Mus_mtDNA0_R are all excluded as chimeras?!?
MA9  <-  MA9[which(lapply(getSequenceTableNoChime(MA9), function(x) sum(x))>0), ]

## annotate taxonomically to make sure they are mouse sequences
MA9 <- blastTaxAnnot(MA9, negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                     db = "/SAN/db/blastdb/nt/nt",
                     infasta = "data/allMouse.fasta",
                     outblast = "data/allMouse.blt",
                     num_threads = 96)

lapply(getTaxonTable(MA9), function(x) dim(x))

## BUG alert! BUG ALERT!  MultiAmplicon::to Phyloseq
## multi2Single=FALSE does not transfer the tax table, if there is one
## amplicon without tax table! ## I solved this above removeing the
## amplicon with only chimears.. but woth to solve!
PSall <- toPhyloseq(MA9, samples=colnames(MA9), multi2Single=FALSE)

PScom <- toPhyloseq(MA9, samples=colnames(MA9), multi2Single=TRUE)


lapply(PSall, function (x) table(tax_table(x)[,"species"], useNA="always"))       
## only 4 ASVs NA all other Mus musculus!!

## how many reads for the top ASVs per amplicson
lapply(PSall, function (foo){       
    tail(unname(taxa_sums(foo)[order(taxa_sums(foo))]))
})


### Let's classify the mice for each amplicon

getDiffPlot <- function(ps, name){
    foo <- otu_table(ps) ### <- 
    foo <- foo[rowSums(foo) > median(rowSums(foo)), ]
    foo <- foo[, colSums(foo) > sort(colSums(foo), decreasing=TRUE)[5]]
    if (all(dim(foo)>2)){
        annDat <- as.data.frame(unclass(sample_data(ps))) ### <- 
        rownames(annDat) <- annDat$Run
        annDat <- annDat[, c("HI", "geo_loc_name", "Latitude")]

        colnames(foo) <- paste(1:ncol(foo))
        pheatmap(log(foo+1), annotation_row=annDat, filename=name)
    } else {message("skipping amp", name, "here\n")}
}

lapply(seq_along(PSall), function (i) {
    getDiffPlot(PSall[[i]], paste0("heatAmp", i, ".pdf"))
})


### -> WOW we see (without having analysed this statistically) rather
### -> North/South differentiation between the Brandenburg/Czeach
### -> tansects, NOT for HI!


### Looking at the markers all together:

getNASVs <- function(ps, name, n=5){
    foo <- as.data.frame(unclass(otu_table(ps)))
    if (all(dim(foo)>2)){
        foo <- foo[rowSums(foo) > median(rowSums(foo)), ]
        foo <- foo[, colSums(foo) > sort(colSums(foo), decreasing=TRUE)[n]]
        colnames(foo) <- paste0("Amp", 1:ncol(foo), name)
        foo$Run <- rownames(foo)
        foo
    } else {message("skipping amp", name, "here\n")}
}

foo <- lapply(seq_along(PSall), function (i) {
    getNASVs(PSall[[i]], names(PSall)[i], n=3)
})

foo <- foo[!unlist(lapply(foo, is.null))]
bar <- Reduce(function (x, y) merge(x, y, all=TRUE, by="Run"), foo)
rownames(bar) <- bar$Run
bar$Run <- NULL
bar[is.na(bar)] <- 0

annDat <- as.data.frame(unclass(getSampleData(MA9)))

rownames(annDat) <- annDat$Run
annDat <- annDat[, c("HI", "geo_loc_name", "Latitude")]

pdf("combindedGentoypeHeat.pdf", height=12, width=18)
pheatmap(log(bar+1), annotation_row=annDat)
dev.off()
