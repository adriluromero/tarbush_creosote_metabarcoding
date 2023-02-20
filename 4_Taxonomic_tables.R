#Jornada_Creosote_Tarbush_Fungal_Metabarcoding

#Protocol for processing Metabarcoding of soil
#Section 4
#File name: 4_Taxonomic_tables
#Step: Use RStudio to process paired-end sequences and generate a taxonomic table.

#In this script, you will: 
# 1. Read-in output ITS2 sequences and perform quality filtering with dada2.
# 2. Assign taxonomy to the sequences using the most recent UNITE release.
# 3. Export the OTU and taxonomic tables. 


#Clear the environment
rm(list = ls(all.names = TRUE))

#Load the packages
require(dada2); packageVersion("dada2")
require(phyloseq)
require(ggplot2)
require(Biostrings)
require(vegan)
require(dplyr)
require(ape)
require(textshape)

#1. Read-in output ITS2 sequences and perform quality filtering with dada2.

#Set the path to where the 'fastq.gz' files are located.
pathF<- "/project/egcc/metabarcoding_tarbush_creosote/dada2/Read_1_ITSx"
pathR<- "/project/egcc/metabarcoding_tarbush_creosote/dada2/Read_2_ITSx"

#Read in the names of the fastq files.
pathF.names<- sort(list.files(pathF, pattern= ".fastq.gz", full.names = TRUE))
pathR.names<- sort(list.files(pathR, pattern= ".fastq.gz", full.names = TRUE))

#Extract sample names where everything before ".fastq.gz" is the sequence name.
pathF.sample.names <- sapply(strsplit(basename(pathF.names), ".fastq.gz"), `[`, 1)
pathR.sample.names <- sapply(strsplit(basename(pathR.names), ".fastq.gz"), `[`, 1)

#Make folders for filtered files.
filtFs1 <- file.path(pathF, "filtered", paste0(pathF.sample.names, "_filt.fastq.gz"))
filtRs1 <- file.path(pathR, "filtered", paste0(pathR.sample.names, "_filt.fastq.gz"))

#Fix names.
names(filtFs1) <- pathF.sample.names
names(filtRs1) <- pathR.sample.names

#Filter sequences.
out1<- filterAndTrim(pathF.names, filtFs1, pathR.names, filtRs1,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, minLen=100, verbose = TRUE, multithread = FALSE)
saveRDS(out1, "/project/egcc/metabarcoding_tarbush_creosote/dada2/Filtered/Filtered.seqs.rds")

#All 60 samples passed the filter.

table(file.exists(filtFs1))
#TRUE 
#60

table(file.exists(filtRs1))
#TRUE 
#60

exists<- file.exists(filtFs1) & file.exists(filtRs1)

#Only keep the files that do exist. In this case, 60 sequences.
filtFs1<-filtFs1[exists]
filtRs1<-filtRs1[exists]

#Learn the error rates for both Forward (F) and Reverse (R) reads.
errF.1 <- learnErrors(filtFs1, multithread=FALSE)
#113407145 total bases in 460952 reads from 6 samples will be used for learning the error rates.

errR.1 <- learnErrors(filtRs1, multithread=FALSE)
#100701232 total bases in 385827 reads from 5 samples will be used for learning the error rates.

saveRDS(errF.1, "/project/egcc/metabarcoding_tarbush_creosote/dada2/Filtered/Forward.errors.rds")
saveRDS(errR.1, "/project/egcc/metabarcoding_tarbush_creosote/dada2/Filtered/Reverse.errors.rds")

#Examine sequence inference.
dada16S_F <- dada(filtFs1, err = errF.1, multithread=FALSE)
dada16S_R <- dada(filtRs1, err = errR.1, multithread=FALSE)
saveRDS(dada16S_F, "/project/egcc/metabarcoding_tarbush_creosote/dada2/Filtered/Forward.denoise.rds")
saveRDS(dada16S_R, "/project/egcc/metabarcoding_tarbush_creosote/dada2/Filtered/Reverse.denoise.rds")

#Merge forward and reverse reads.
mergers <- mergePairs(dada16S_F, filtFs1, dada16S_R, filtRs1, verbose=TRUE)
saveRDS(mergers, "/project/egcc/metabarcoding_tarbush_creosote/dada2/Filtered/Merged.Seqs.rds")

#Make a sequence table.
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/project/egcc/metabarcoding_tarbush_creosote/dada2/Filtered/Seq.tab.rds")

#Remove chimeras.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE) 

#Identified 1885 bimeras out of 10548 input sequences.

saveRDS(seqtab.nochim,"/project/egcc/metabarcoding_tarbush_creosote/dada2/Filtered/No.chimera.seq.tab.rds")

seqtab.nochim.length.filter <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% 164:600]
saveRDS(seqtab.nochim.length.filter,"/project/egcc/metabarcoding_tarbush_creosote/dada2/Filtered/No.chimera.short.removed.seq.tab.rds")

#How many sequences were retained at each QC step? 
getN <- function(x) sum(getUniques(x))
track <- cbind(out1, sapply(dada16S_F, getN), 
               sapply(dada16S_R, getN), 
               sapply(mergers,getN), 
               rowSums(seqtab.nochim),
               rowSums(seqtab.nochim.length.filter))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim", "length.filter")

View(track)
#The QC parameters look good with no major loss at any step.

#2. Assign taxonomy to the sequences using the most recent UNITE release.

#Download the UNITE general FASTA release for Fungi from https://unite.ut.ee/repository.php 
#At the time of running this pipeline, the most recent UNITE release was created in 2022-10-16.
#This particular release was downloaded from https://doi.org/10.15156/BIO/2483911

#Upload and unpack the UNITE file (*.tgz) into your directory. 
#There will be a 'sh_general_release_dynamic_29.11.2022.fasta' and 'sh_general_release_dynamic_29.11.2022_dev.fasta' file. 
#Use the 'sh_general_release_dynamic_29.11.2022.fasta' file. 

unite.ref <- "/project/egcc/metabarcoding_tarbush_creosote/UNITE/sh_general_release_dynamic_29.11.2022.fasta"
taxa <- assignTaxonomy(seqtab.nochim.length.filter, unite.ref, multithread = FALSE, tryRC = TRUE)
#This step of assigning taxonomy can take hours.

#Check the results. 
taxa.print <- taxa 
rownames(taxa.print) <- NULL 
head(taxa.print)

#which phyla are present?
unique(taxa.print[,2])

# [1] "p__Basidiomycota"            "p__Ascomycota"              
# [3] "p__Mortierellomycota"        "p__Chytridiomycota"         
# [5] "p__Fungi_phy_Incertae_sedis" "p__Calcarisporiellomycota"  
# [7] NA                            "p__Basidiobolomycota"       
# [9] "p__Glomeromycota"            "p__Monoblepharomycota"      
#[11] "p__Mucoromycota"            

#3. Export the OTU and taxonomic tables. 
ps <- phyloseq(otu_table(seqtab.nochim.length.filter, taxa_are_rows=FALSE), 
               tax_table(taxa))

#Rename ASV sequences as ASV# 
#save as a string in case needed later-on as DNA
dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 8663 taxa and 60 samples ]
#tax_table()   Taxonomy Table:    [ 8663 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 8663 reference sequences ]

#Taxonomic filtering
rank_names(ps)
#[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

#Create a table with the number of features for each phyla

table(tax_table(ps)[, "Phylum"], exclude = NULL)

#              p__Ascomycota        p__Basidiobolomycota 
#                       7462                           4 
#           p__Basidiomycota   p__Calcarisporiellomycota 
#                        690                           3 
#         p__Chytridiomycota p__Fungi_phy_Incertae_sedis 
#                        313                          62 
#           p__Glomeromycota       p__Monoblepharomycota 
#                         15                           4 
#       p__Mortierellomycota             p__Mucoromycota 
#                         33                           4 
#                       <NA> 
#                         73 

#Filter out potential artifacts
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "p__Monoblepharomycota", "p__Glomeromycota", "p__Calcarisporiellomycota", "p__Basidiobolomycota"))

#Now, make these into data frames and save the files.
otu.table.not.rare<- data.frame(t(data.frame(otu_table(ps))))
#Look at the sequencing depth of each sample and sort column sums. 
OTU_Col_Sums <- sort(colSums(otu.table.not.rare)) 
tax.table<- data.frame(tax_table(ps))
identical(row.names(otu.table.not.rare), row.names(tax.table)) #row names match up perfect
otu.table.tax.table<- cbind(otu.table.not.rare, tax.table)
write.csv(otu.table.tax.table, "/project/egcc/metabarcoding_tarbush_creosote/OTU-tables/2023-OTUs/otu.tax.table.csv")
write.csv(OTU_Col_Sums, "/project/egcc/metabarcoding_tarbush_creosote/OTU-tables/2023-OTUs/otu.col.sums.csv")