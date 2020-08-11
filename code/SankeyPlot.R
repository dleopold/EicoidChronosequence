#' ---
#' title: Sankey plot
#' author: Devin R Leopold
#' date: April 13, 2020
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: false
#'      self_contained: true
#'      highlight: zenburn
#' ---

#' This script manually constructs a Sankey plot linking the most abundant OTUs to the LSAG chronosequence sites to visualize difference in relative abundance. 

#' # Prepare environment
library(tidyverse)
library(magrittr)
library(phyloseq)
library(foreach)
library(ggforce)
library(labdsv)
library(ggtext)
source("code/Rfunctions.R")

#' Colors for figure
SiteColors <- rev(viridis::viridis_pal(direction=-1,option="inferno")(10)[c(1,2,3,5,8)])

#' Import data and convert to proportional abundance
phy <- readRDS("output/rds/phy.rds") %>% 
  subset_samples(LSAG) %>%
  sweepOTUs %>%
  transform_sample_counts(function(x){x/sum(x)}) 

#' identify top 50 OTUs for plotting
topOTUs = names(sort(taxa_sums(phy), TRUE)[1:50])
phy %<>% prune_taxa(topOTUs,.) %>% transform_sample_counts(function(x){x/sum(x)}) 

#' prepare data in table form
dat_wide <- phy %>% otu_table %>% data.frame %>% 
  rownames_to_column("Site") %>%
  mutate(Site=factor(Site,levels=Site))

#' # indicator species analysis
ind.res <- indval(dat_wide[,-1],sample_data(phy)$Site,numitr=10000)
ind.dat <- data.frame(Taxon=names(ind.res$maxcls),
                      maxcls=levels(sample_data(phy)$Site)[ind.res$maxcls],
                      pvals=ind.res$pval %>% p.adjust("fdr"),
                      stringsAsFactors = F) %>%
  transmute(Taxon=Taxon,
            maxcls=maxcls,
            indicator=ifelse(pvals<.05,maxcls,NA))

#' add indicator data to a tidy summarized data set for plotting
dat_long <- phy %>% merge_samples("Site") %>%
  otu_table %>% data.frame %>% 
  rownames_to_column("Site") %>% 
  mutate(Site=factor(Site,levels=Site)) %>%
  pivot_longer(2:51,names_to="Taxon") %>%
  left_join(ind.dat) 


#' # Parallel set data
#' Manually construct data for plotting as parallel_sets using ggforce.  
#'   
#' Set the order of OTU plotting (sort by site with highest ind val and them by relative abundance)
OTUorder <- foreach(site=levels(dat_long$Site), .combine=c) %do% {
  dat_long %>% filter(maxcls==site) %>% 
    group_by(Taxon) %>% summarize(value=sum(value)) %>%
    arrange(desc(value)) %$% Taxon
}

#' Set parameters and counters for defining ranges for each taxon
siteHeight <- dat_long %>% filter(Site=="Kokee") %$% value %>% sum %>% multiply_by(9)
sep1 <- 50 # Site separation
sep2 <- 14 # Max OTU separation
sep3 <- 5 # Min OTU separation
#' function to scale variation in OTU width for label spacing
OTUrange <- dat_wide[,-1] %>% colSums() %>% range
spacingScaler <- function(x){
  ((sep3-sep2)/diff(OTUrange)*(x-max(OTUrange))+sep3)
}
#' Initiate counters for loop
siteCounter <- list(Kokee=siteHeight,
                    Kohala=2*siteHeight+sep1,
                    Laupahoehoe=3*siteHeight+2*sep1,
                    Olaa=4*siteHeight+3*sep1,
                    Thurston=5*siteHeight+4*sep1)
OTUcounter <- sum(dat_long$value) + sep2*(length(unique(dat_long$Taxon))-1)

#' Loop over each OTU (top to bottom in plotting order) and define coodinates in each set
parallel.dat <- foreach(otu=OTUorder, .combine = bind_rows) %do% {
  otu.dat <- filter(dat_long,Taxon==otu)
  y <- foreach(site=levels(otu.dat$Site), .combine=c) %do% {
    width <- 9*otu.dat$value[otu.dat$Site==site]
    y=c(siteCounter[[site]],siteCounter[[site]]-width,
        OTUcounter,OTUcounter-width)
    siteCounter[[site]] %<>% subtract(width)
    OTUcounter %<>% subtract(width)
    return(y)
  }
  OTUcounter %<>% subtract(spacingScaler(sum(otu.dat$value)))
  data.frame(x=rep(c(1.5,1.5,2,2),5),
             y=y,
             group=rep(paste(otu.dat$Taxon,otu.dat$Site,sep="."),each=4),
             Site=rep(otu.dat$Site,each=4),
             Taxon=rep(otu.dat$Taxon,each=4),
             indicator=otu.dat$indicator,
             stringsAsFactors = F)
}
#' Center plotting ranges for each set
parallel.dat %<>% mutate_by(x,y=y-((max(y)+min(y))/2))
#' Make subset of only indicator species for overplotting
parallel.dat.ind <- parallel.dat %>% filter(Site==indicator)

#' Define coordinates for site and OTU nodes in plot
siteShapes <- parallel.dat %>% filter(x==1.5) %>%
  group_by(Site) %>%
  summarise(ymin=min(y),ymax=max(y),center=mean(c(ymin,ymax)))
otuShapes <- parallel.dat %>% filter(x==2) %>%
  group_by(Taxon) %>%
  summarize(ymin=min(y),ymax=max(y),indicator=indicator[1])
#' Add site labels
siteShapes$age <- c("3.0&times;10^2 yr",
                    "2.1&times;10^3 yr",
                    "2.0&times;10^4 yr",
                    "1.5&times;10^5 yr",
                    "4.1&times;10^6 yr")
siteShapes %<>% mutate(lab=paste0(Site,"<br>",age)) 

#' Format IDs and taxonomy strings fro potting (this code needs to be streamlined)
otuLabs <- tax_table(phy) %>% data.frame %>%
  rownames_to_column("Taxon") 
otuLabs$label <- apply(otuLabs,1,function(x){
  x[!is.na(x)] %>% .[length(.)]
}) %>% gsub("_sp","",.) %>% gsub("_fam_.*","",.)
otuLabs$labTest <- apply(otuLabs,1,function(x){
  which(x[-9]==x[9])
})
otuLabs %<>% transmute(Taxon=Taxon,
                       label=ifelse(labTest>6,paste0("*",label,"*"),label) %>%
                         gsub("_"," ",.))
otuLabs <- parallel.dat %>% filter(x==2) %>%
  group_by(Taxon) %>% summarize(y=mean(c(min(y),max(y)))) %>%
  full_join(otuLabs)
otuLabs$label %<>% paste("|",.)
otuLabs$label %<>% gsub("Trechisporales","Trechisporales*",.)

#erm <- readRDS("output/rds/phy.erm.rds") %>% taxa_names()
#otuLabs$Taxon[otuLabs$Taxon %in% erm] %<>% paste0("**",.,"**")

#' Make plot
sankey <- ggplot(parallel.dat) +
  geom_diagonal_wide(aes(x, y, group = group), alpha=0.5,
                     strength=0.5, fill="grey60")+
  geom_diagonal_wide(data=parallel.dat.ind, color="grey40", strength=0.5,size=0.25,
                     aes(x, y, group = group, fill=indicator))+
  geom_rect(data=siteShapes, color="grey40",size=0.25,
            aes(xmin=1.4,xmax=1.5,ymin=ymin,ymax=ymax,fill=Site))+
  geom_rect(data=otuShapes, fill="grey60",
            aes(xmin=2,xmax=2.05,ymin=ymin,ymax=ymax))+
  geom_rect(data=drop_na(otuShapes), color="grey40",size=0.25,
            aes(xmin=2,xmax=2.05,ymin=ymin,ymax=ymax,fill=indicator))+
  scale_fill_manual("Site age (yr)",
                    values=SiteColors)+
  geom_richtext(data=otuLabs,aes(x=2.07,label=Taxon,y=y), size=4,
                hjust=0, vjust=0.75, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0, 4), "pt"))+
  geom_richtext(data=siteShapes, aes(x=1.39,y=center,label=lab),
                hjust=1,fill = NA, label.color = NA, size=4.8)+
  geom_richtext(data=data.frame(),aes(x=2.315,y=610, label="OTU.id | predicted taxonomy"),
                fill = NA, label.color = NA, size=5) +
  geom_richtext(data=otuLabs,aes(x=2.19,label=label,y=y), size=4,
                hjust=0, vjust=0.75, fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt"))+
  lims(x=c(1.15,2.58)) +
  cowplot::theme_map() +
  theme(legend.position = 'none',
        legend.title = element_markdown(size=16),
        legend.text = element_markdown(size=14),
        legend.key.size = unit(1.2,"line"),
        plot.margin = margin(0,0,0,0)) 
ggsave("output/figs/Sankey.pdf",sankey,width=18,height=28,units = "cm")


