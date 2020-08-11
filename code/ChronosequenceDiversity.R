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

#' Read in and subset data
phy.all <- readRDS("output/rds/phy.rds") %>% 
  subset_samples(LSAG) %>%
  sweepOTUs 
phy.erm <- readRDS("output/rds/phy.erm.rds") %>% 
  subset_samples(LSAG) %>% 
  prune_samples(sample_names(.)!="B10",.) %>% # Remove one poorly sequenced sample
  sweepOTUs 

#' # Estimate diversity
div.dat <- bind_rows(iNEXT.phy(phy.all, level=5000) %>% mutate(Taxa="All fungi"),
                     iNEXT.phy(phy.erm, level=5000) %>% mutate(Taxa="Putative ErM")) 

#' # Global tests
#' Test whether richness and diversity vary among chronosequence sites
div.stat <- div.dat %>%
  select(Metric, Site, Taxa, Estimate) %>%
  group_by(Metric,Taxa) %>%
  group_map(~ aov_chrono(.x,.y)) %>%
  bind_rows
gt(div.stat)

#' # Post-hoc tests
pairwise <- div.dat %>%
  select(Metric,Site,Taxa, Estimate) %>%
  group_by(Metric,Taxa) %>%
  group_map(~tukeyWrap(.x,.y)) %>%
  bind_rows
#' add plotting locations
pairwise %<>%
  left_join(
    div.dat %>% 
      group_by(Metric,Site,Taxa) %>%
      summarise(Estimate=1+quantile(Estimate,0.85,names=F,type=5))
  )

#' # Plot
#' Bootstrap distributions for plotting uncertainty 
boot.dat <- div.dat %>%
  select(Metric,Site,Taxa, Estimate) %>%
  group_by(Metric,Taxa,Site) %>%
  group_map(~bootFun(.x,.y)) %>%
  bind_rows %>%
  group_by(Metric,Taxa,Site) %>%
  summarize(mean=mean(Estimate),
            LCI=coxed::bca(Estimate)[1],
            UCI=coxed::bca(Estimate)[2])

#' Make plot
plotter <- function(taxa){
  dat <- filter(div.dat,Taxa==taxa) %>%
    mutate(Metric=factor(Metric, levels=(unique(Metric))))
  dat.means <- filter(boot.dat,Taxa==taxa) %>%
    mutate(Metric=factor(Metric,levels=levels(dat$Metric)))
  lets <- filter(pairwise, Taxa==taxa) %>%
    mutate(Metric=factor(Metric,levels=levels(dat$Metric)))
  title <- taxa
  ggplot(dat,aes(x=Site,Estimate)) +
    geom_boxplot(color="grey65")+
    scale_color_gradient(low="grey10",high="grey80")+
    geom_pointrange(data=dat.means, aes(y=mean, ymin=LCI, ymax=UCI),
                    shape=21, size=1, fatten=2, stroke=0.8, fill="white")+
    scale_x_discrete(labels=c("Thurston",
                              "Ola'a",
                              "Laupahoehoe",
                              "Kohala",
                              "Kokee"),
                     name="Site age (yr)")+
    geom_text(data=lets, aes(label=Letters), nudge_x = -0.25)+
    geom_hline(aes(yintercept=0))+
    geom_text(aes(x=1,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=2,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=3,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=4,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=5,y=0,label="|"),vjust=1,size=1.75)+
    ggtitle(taxa)+
    facet_wrap(~Metric,scales="free_y",ncol=1,strip.position ="left") + 
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
ggsave("output/figs/divPlot.pdf",width=16,height=11,units="cm")







