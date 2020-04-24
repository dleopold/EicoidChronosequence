#' ---
#' title: Chronosequence diversity
#' author: Devin R Leopold
#' date: April 20, 2020
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
library(phyloseq)
library(ggtext)
library(gt)
library(cowplot)
source("code/Rfunctions.R")
library(iNEXT)
library(breakaway)

#' Read in and subset data
phy.all <- readRDS("output/rds/phy.rds") %>% 
  subset_samples(LSAG) %>%
  sweepOTUs 
phy.erm <- readRDS("output/rds/phy.erm.rds") %>% 
  subset_samples(LSAG) %>%
  sweepOTUs 

#' # Estimate diversity
div.dat <- bind_rows(iNEXT.phy(phy.all) %>% mutate(Taxa="All fungi"),
                     iNEXT.phy(phy.erm) %>% mutate(Taxa="Putative ErM")) %>% 
  filter(Diversity!="Simpson diversity") %>% # we will focus on Shannon diversity and species richness
  mutate_by(Diversity,Site,Taxa,
            seScale=scale(s.e.))

#' # Global test
#' Test whether richness and diversity vary among chronosequence sites
(div.stat <- div.dat %>%
  select(Diversity,Site,Taxa, Estimator, s.e.) %>%
  group_by(Diversity,Taxa) %>%
  group_modify(~betta_LRT(.x))) %>%
  gt()

#' # Post-hoc test
#' For metrics that vary among sites, test which site are different
(pairwise <- div.dat %>%
  filter(Diversity=="Shannon diversity") %>% 
  select(Diversity,Site,Taxa, Estimator, s.e.) %>%
  group_by(Diversity,Taxa) %>%
  group_modify(~betta_pairwise(.x)))
#' add plotting locations
pairwise %<>%
  left_join(
    div.dat %>% 
      group_by(Diversity,Site,Taxa) %>%
      summarise(Estimator=1+quantile(Estimator,0.85,names=F,type=5))
  )

#' # Plot
#' Get the predicted means, pooling information across samples and accounting for uncrtainty in point estimates.
(div.means <- div.dat %>%
    select(Diversity,Site,Taxa, Estimator, s.e.) %>%
    group_by(Diversity,Taxa) %>%
    group_modify(~betta_means(.x, "~0+Site")) %>%
    mutate(Predictor=gsub("Site","",Predictor)) %>%
    rename(Site=Predictor) %>%
    mutate(Site=factor(Site,levels=levels(div.dat$Site)))) %>%
  gt()

#' Make plot

#' Colors for figure
SiteColors <- viridis::viridis_pal(direction=-1,option="inferno")(10)[c(1,2,3,5,8)]

plotter <- function(taxa){
  dat <- filter(div.dat,Taxa==taxa)
  dat.means <- filter(div.means,Taxa==taxa)
  lets <- filter(pairwise, Taxa==taxa)
  title <- taxa
  ggplot(dat,aes(x=Site,Estimator)) +
    geom_boxplot(color="grey65")+
    scale_color_gradient(low="grey10",high="grey80")+
    geom_pointrange(data=dat.means, aes(y=mean, ymin=LCI, ymax=UCI),
                    shape=21, size=1, fatten=2.5, stroke=0.8, fill="white")+
    scale_fill_manual(values=SiteColors)+
    scale_x_discrete(labels=c("Thurston",
                              "Ola'a",
                              "Laupahoehoe",
                              "Kohala",
                              "Kokee"),
                     name="Site age (yr)")+
    geom_text(data=lets, aes(label=Let), nudge_x = -0.25)+
    geom_hline(aes(yintercept=0))+
    geom_text(aes(x=1,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=2,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=3,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=4,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=5,y=0,label="|"),vjust=1,size=1.75)+
    ggtitle(taxa)+
    facet_wrap(~Diversity,scales="free_y",ncol=1,strip.position ="left") + 
    expand_limits(y=0) +
    theme_cowplot() +
    theme(axis.title = element_blank(),
          axis.text.x = element_markdown(angle=35,hjust=1,vjust=1),
          axis.text.y = element_text(size=10),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(size=12),
          legend.position = "top",
          legend.title = element_blank(),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5))
}
all <- plotter("All fungi")
erm <- plotter("Putative ErM")
plot_grid(all,erm+theme(strip.text.y = element_text(color="white")),nrow=1,labels=c("(a)","(b)"))
ggsave("output/figs/Fig.4.pdf",width=16,height=11,units="cm")







