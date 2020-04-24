#' ---
#' title: MiSeq data processing
#' author: Devin R Leopold
#' date: April 10, 2020
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: false
#'      self_contained: true
#'      highlight: zenburn
#' ---

#' # Prepare environment
library(tidyverse)
library(magrittr)
library(dada2)
library(Biostrings)
library(phyloseq)

#' read in and prepare sample meta data
mapping <- read.delim("data/mapping.tab",header=T,row.names=1,as.is=T)
mapping %<>% 
  rownames_to_column("SampID") %>%
  mutate(
    logAge=log(Age*1000),
    Site=factor(Site,levels=c("Thurston","Olaa","Laupahoehoe","Kohala","Kokee")),
    N=grepl("N",Frt_treat),
    P=grepl("P",Frt_treat)
  ) %>%
  column_to_rownames("SampID")

#' identify gene primers for trimming
ITS1F <- "CTTGGTCATTTAGAGGAAGTAA"
ITS1F.rc <- dada2::rc(ITS1F)
ITS2 <- "GCTGCGTTCTTCATCGATGC"
ITS2.rc <- dada2::rc(ITS2)

#' identify paths to raw data
raw.fwd <- list.files("data/MiSeq_raw/",pattern="R1_001.fastq.gz") %>% sort 
raw.rev <- list.files("data/MiSeq_raw/",pattern="R2_001.fastq.gz") %>% sort 

#' identify sample IDs
samps <- raw.fwd %>% strsplit("_") %>% sapply(`[`,1)

#' # Trim reads
if(!dir.exists("output/trim/")){
  dir.create("output/trim/", F, T)
  for(samp in samps){
    in.raw <- paste0(samp,"_") %>% grep(raw.fwd, value = T) %>% list.files("data/MiSeq_raw",pattern=., full.names = T)
    out.trim <- paste0(samp,".R1.fastq.gz") %>% file.path("output/trim",.)
    
    flags <- paste("-o", out.trim,
                   "-a", ITS2.rc,
                   "-e 0.2",
                   "--minimum-length 100","-l 250",
                   "--quality-cutoff 15",
                   "--quiet",
                   in.raw)
    system2("/home/harpua/miniconda3/bin/cutadapt",args=flags)
  }
}

#' identify paths to trimmed reads
trimmed.fwd <- list.files("output/trim",pattern="R1.fastq", full.names = T) %>% sort 

#' preview read quality of 6 randomly selected samples
qual.samps <- sample(1:length(trimmed.fwd),6)
#+ fig.align="center", dpi=49, out.width="600px"
plotQualityProfile(trimmed.fwd[qual.samps])

#' # Quality filter
dir.create("output/dada/filt",F,T)
#' idenitfy file paths for filtering output
filt.fwd <- trimmed.fwd %>% gsub("trim","dada/filt",.)
#+ filt
filt <- filterAndTrim(trimmed.fwd,filt.fwd,
                      truncQ = 2, minLen = 75, maxN=0, maxEE=2, trimLeft=10, 
                      rm.phix=TRUE, compress=TRUE, multithread=TRUE)
filt
#' check quality profiles after trimming
#+ fig.align="center", dpi=49, out.width="600px"
plotQualityProfile(filt.fwd[qual.samps])

#' # Dereplicate
#+ Derep
derep.fwd <- derepFastq(filt.fwd, verbose=F)

#' Trim names of derep objects to just sample names 
names(derep.fwd) %<>% gsub(".R1.fastq.gz","",.)

#' # Learn errors 
#+ Learn, fig.align="center", dpi=49, out.width="600px"
set.seed(808080)
if(!file.exists("output/rds/errMod.rds")){
  err.fwd <- learnErrors(filt.fwd, multithread=TRUE, nbases = 5e08, randomize=T) 
  plotErrors(err.fwd, nominalQ=TRUE)
  saveRDS(err.fwd,"output/rds/errMod.rds")
}else{err.fwd <- readRDS("output/rds/errMod.rds")}

#' # Denoise
#+ Denoise
if(!file.exists("output/rds/dada.rds")){
  dada.fwd <- dada(derep.fwd, err=err.fwd, multithread=TRUE, pool=T,
                   GREEDY=F, OMEGA_A=1e-45, MIN_ABUNDANCE=4)
  saveRDS(dada.fwd,"output/rds/dada.rds")
}else(dada.fwd <- readRDS("output/rds/dada.rds"))

#' Make sequence table
seqtab <- dada.fwd  %>% makeSequenceTable %>% collapseNoMismatch

#' # Remove chimeras
seqtab %<>% removeBimeraDenovo(method="consensus", multithread=TRUE,verbose=T)

#' # Make phyloseq object

#' Extract sequences
seqs <- getSequences(seqtab) %>% DNAStringSet()
names(seqs) <- getSequences(seqtab)

#' Assign taxonomy with UNITE fungal database v8.2
taxa <- assignTaxonomy(seqs,"data/taxDB/sh_general_release_dynamic_s_04.02.2020.fasta", multithread = T)

#' combine
(phy <- phyloseq(otu_table(seqtab,taxa_are_rows = F),
                refseq(seqs),
                tax_table(taxa),
                sample_data(mapping)))

#' # cluster OTUs 
cluster <- function(phy.in,method="single",dissimilarity=0.01){
  require(DECIPHER)
  clusters <- DistanceMatrix(refseq(phy.in), includeTerminalGaps = T, processors=NULL) %>%
    IdClusters(method=method, cutoff=dissimilarity, processors=NULL) 
  clusters[,1] %<>% as.character()
  tax.holder <- tax_table(phy.in)
  tax_table(phy.in) <- clusters %>% as.matrix %>% tax_table
  phy.in %<>% speedyseq::tax_glom("cluster")
  tax_table(phy.in) <- tax.holder
  return(phy.in)
}  
phy %<>% cluster  

#' # Taxonomy filter
#' Use UNITE "all Eukaryotes" database to filter non-fungal amplicons
taxa.filter <- assignTaxonomy(refseq(phy),"data/taxDB/sh_general_release_dynamic_all_04.02.2020.fasta", multithread = T)
(phy %<>% prune_taxa(grepl("k__Fungi",taxa.filter[,1]),.))

#' # Negative controls
#' max seq count in neg controls
maxNegAbund <- phy %>% subset_samples(Control=="neg") %>% otu_table %>% data.frame %>% apply(2,max)

#' find mean seq abund in env. samps
meanAbund <- phy %>% subset_samples(is.na(Control)) %>% otu_table %>% apply(2,mean) 

#' filter likely contaminants (max abund in control samps greater than mean abund in env. samps.)
phy %<>% prune_taxa(maxNegAbund < meanAbund,.) %>%
  subset_samples(is.na(Control)) #Also remove control samps. from phy object

#' set ascending OTU names
tmp.names <- order(taxa_sums(phy),decreasing = T)
names(tmp.names) <- paste0("OTU.",1:ntaxa(phy)) 
taxa_names(phy) <- names(sort(tmp.names))

#' # Subset ERM
#' identify taxa that are in likely ERM orders)
ERM_bytax <- tax_table(phy)[,'Order'] %in% c("Chaetothyriales","Helotiales","Sebacinales","Trechisporales") |
  tax_table(phy)[,'Genus'] %in% c("Oidiodendron","Leohumicola") 
#' subset phyloseq object
(phy.erm <- phy %>% prune_taxa(ERM_bytax,.))

#' # Summary
#' **Total # of sequences =** `r phy %>% sample_sums() %>% sum`  
#' **Total # of ERM sequences =** `r phy.erm %>% sample_sums() %>% sum`   
#' **Proportion of OTUs in ERM subset =** `r ntaxa(phy.erm)/ntaxa(phy)`  
#' **Proportion of sequences in ERM subset =** `r (phy.erm %>% sample_sums() %>% sum)/(phy %>% sample_sums() %>% sum)`  

#' # Save output
dir.create("output/rds",F)
saveRDS(phy,"output/rds/phy.rds")
saveRDS(phy.erm,"output/rds/phy.erm.rds")

