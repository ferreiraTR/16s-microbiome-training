######### 
###NIAID 16S Microbiome Training 
#########


## This script lists all of the packages that need to be installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("dada2")
BiocManager::install("phangorn")
BiocManager::install("phyloseq")
BiocManager::install("DECIPHER")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")

install.packages("ggplot2")

##Find and Organize our Data
data_location <- "16S-data-processing/raw_files"

library(dada2)

read_indicator <- "_1" #can be different deppending on sequencing (e.g. _R1)

all_fastqs <- list.files(data_location, full.names = T)

r1_fastqs <- all_fastqs[grepl(read_indicator, all_fastqs)]
r2_fastqs <- all_fastqs[!grepl(read_indicator, all_fastqs)]
r1_fastqs[1:10]

##Identifying primers and adapters

library(ShortRead)
library(Biostrings)
FWD <- "GTGCCAGCMGCCGCGGTAA"
REV <- "GGACTACHVGGGTWTCTAAT"

primerHits <- function(primer, fn) {
  start_seq <- sread(readFastq(fn))
  start_seq <- subseq(start_seq, start=1, end=length(primer)+1)
  nhits <- vcountPattern(primer, start_seq, fixed = FALSE)
  return(sum(nhits > 0))
}
sapply(FWD, primerHits, r1_fastqs[1:10])
sapply(REV, primerHits, r2_fastqs[1:10])


##Quality FIltering and Trimming

plotQualityProfile(r1_fastqs[1:10])
plotQualityProfile(r2_fastqs[1])

r1_filt <- gsub("raw_files", "filtered", r1_fastqs)
r2_filt <- gsub("raw_files", "filtered", r2_fastqs)
## In my case, the filtered directory could exist before I make the document -- this `unlink` function just makes sure it's not there when I run filterAndTrim
unlink("filtered/", recursive = T)
out <- filterAndTrim(r1_fastqs, r1_filt, r2_fastqs, r2_filt, truncLen=0, maxN=0, maxEE=2, multithread = F, verbose=T, rm.lowcomplex = 4, rm.phix = T) #can use trimleft to remove bases on the left, can use truncLen to truncate to a specific length, maxLen = max length, truncQ = quality truncation

##Denois Reads

r1_err <- learnErrors(r1_filt, multithread=F, verbose=T)
r2_err <- learnErrors(r2_filt, multithread=F, verbose=T)
plotErrors(r1_err, nominalQ = T)

r1_dada <- dada(r1_filt, r1_err)
r2_dada <- dada(r2_filt, r2_err)

##Merge Reads
#Here, we need to keep in mind our sequencing structure – for the data that we want to work with today, we need to keep in mind that our region of interest is 253 bases, which we can cover with 153 bases on the forward read and 153 bases on the reverse read. Therefore, 306-253=53, suggesting we can’t expect an overlap greater than ~53 bases and we need to keep that in mind.
merged_reads <- mergePairs(r1_dada, r1_filt, r2_dada, r2_filt, verbose=TRUE, minOverlap = 20)

#Merge reads option 2
merged_reads2 <- mergePairs(r1_dada, r1_filt, r2_dada, r2_filt, verbose=TRUE, minOverlap = 20, maxMismatch = 2)

##Make Sequence Table
seq_table <- makeSequenceTable(merged_reads)
dim(seq_table)
table(nchar(getSequences(seq_table)))


##Remove chimeras - 

chimera_check <- removeBimeraDenovo(seq_table, verbose=T)
dim(chimera_check)
sum(chimera_check)/sum(seq_table) #how many reads were retained


##Workflow checkpoint
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(r1_dada, getN), sapply(r2_dada, getN), sapply(merged_reads, getN), rowSums(chimera_check))
colnames(track) <- c("input","filtered","denoised_r1","denoised_r2","merged","nonchimeric")
track



##Assign Taxonomy
taxa <- assignTaxonomy(chimera_check, "16S-data-processing/databases/silva_nr99_v138.1_train_set.fa.gz", multithread=F, verbose=T, minBoot = 80)

taxa_print <- taxa
rownames(taxa_print) <- NULL
head(taxa_print,10)


##Make Phylogenetic Tree

library(DECIPHER)

library(phangorn)
seqs <- as.character(getSequences(chimera_check))
names(seqs) <- seqs
alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- as.phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit <- pml(treeNJ, data=phang.align)

library(ape)
ape::write.tree(fit$tree, "saved_tree.txt")

##Generate Phyloseq Object

library(phyloseq)
dummy_sample_data <- data.frame(samples=rownames(chimera_check), treatment_group=sample(c("A","B"), size=20, replace=TRUE))
rownames(dummy_sample_data) <- dummy_sample_data$samples

phy_obj <- phyloseq(otu_table(chimera_check, taxa_are_rows=F), sample_data(dummy_sample_data), phy_tree(fit$tree), tax_table(taxa))
phy_obj

##VIsualizations
library(ggplot2)
plot_bar(phy_obj, fill="Genus") + guides(fill="none")

ord <- ordinate(phy_obj, method="PCoA", distance = "bray")
plot_ordination(phy_obj, ord)

ord <- ordinate(phy_obj, method="PCoA", distance = "wunifrac")

plot_ordination(phy_obj, ord)

#tree visualization
library(ggtree) #also adds geom_tree to ggplot2
ggtree(fit$tree, layout = "daylight", branch.length = 'none')
ggtree(fit$tree, layout = "equal_angle", branch.length = 'none')


