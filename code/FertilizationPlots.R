#' ---
#' title: Fertilization plots
#' author: Devin R Leopold
#' date: April 21, 2020
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
library(gt)
library(cowplot)
library(ggtext)
source("code/Rfunctions.R")
library(iNEXT)
library(breakaway)

#' Read in and subset data
phy.all <- readRDS("output/rds/phy.rds") %>% 
  subset_samples(!is.na(Frt_plot)) %>%
  sweepOTUs 
phy.erm <- readRDS("output/rds/phy.erm.rds") %>% 
  subset_samples(!is.na(Frt_plot)) %>%
  sweepOTUs 

#' # Estimate diversity
div.dat <- bind_rows(iNEXT.phy(phy.all) %>% mutate(Taxa="All fungi"),
                     iNEXT.phy(phy.erm) %>% mutate(Taxa="Putative ErM")) %>% 
  filter(Diversity!="Simpson diversity") %>%# we will focus on Shannon diversity and species richness
  mutate(Frt_treat=factor(Frt_treat,levels=c("C","N","P","NP"))) %>%
  mutate_by(Diversity,Site,Taxa,
            seScale=scale(s.e.))
  
#' # Test fertilization effects 
(div.stat <- div.dat %>%
    select(Diversity,Site,Taxa, Estimator, s.e., N, P) %>%
    group_by(Site,Diversity,Taxa) %>%
    group_modify(~betta_wrap(.x)) %>%
    filter(Predictor!="(Intercept)")) %>%
  gt()

#' # Plot
#' Get mean predicted values for plotting
div.means <- div.dat %>%
  select(Diversity,Site,Taxa, Estimator, s.e., Frt_treat) %>%
  group_by(Site,Diversity,Taxa) %>%
  group_modify(~betta_means(.x, "~0+Frt_treat")) %>%
  mutate(Predictor=gsub("Frt_treat","",Predictor),
         Predictor=factor(Predictor,levels=levels(div.dat$Frt_treat)))
#' We will make separate plots for the 2 fertilized experiments and then combine them
plotter <- function(site){
  dat <- filter(div.dat, Site==site)
  means <- filter(div.means, Site==site)
  title <- ifelse(site=="Thurston","Young (300 yr)","Old (4.1 myr)")
  stats <- filter(div.stat, Site==site) %>%
    mutate(stars=gtools::stars.pval(p.values) %>% 
             gsub("*","&#42;",.,fixed=T) %>%
             gsub(".","+",.,fixed=T),
           Predictor=gsub(":","*&times;*",Predictor,fixed=T),
           label=paste0("*",Predictor,"*", stars)) %>%
    group_by(Diversity,Taxa) %>%
    summarize(label=paste(label,collapse = "<br>"))
  ggplot(dat,aes(x=Frt_treat,y=Estimator))+
    geom_boxplot(color="grey65")+
    #ggbeeswarm::geom_quasirandom(size=1.5,position=position_jitter(width = 0.2),
    #                             shape=19,aes(color=seScale), show.legend = F, alpha=0.5)+
    scale_color_gradient(low="grey20",high="grey80")+
    geom_pointrange(data=means, shape=21, fill="white", size=1, fatten=2,
                    aes(ymin=LCI,ymax=UCI,y=mean,x=Predictor)) +
    geom_richtext(data=stats,aes(x=4.6,y=0,label=label),
                  vjust=0,hjust=-0,fill = NA, label.color = NA)+
    coord_cartesian(clip = 'off') +
    facet_grid(Diversity~Taxa,scales="free",
               switch="y")+
    scale_x_discrete("Fertilization treatment")+
    geom_hline(aes(yintercept=0))+
    geom_text(aes(x=1,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=2,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=3,y=0,label="|"),vjust=1,size=1.75)+
    geom_text(aes(x=4,y=0,label="|"),vjust=1,size=1.75)+
    expand_limits(y=0)+
    ggtitle(title)+
    theme_cowplot() +
    theme(panel.spacing = unit(2.75, "lines"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust=0.5),
          axis.line.x = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          plot.margin = unit(c(7,50,7,7), "pt"))
}
(young <- plotter("Thurston"))
(old <- plotter("Kokee")) + theme(strip.text.y = element_blank())
plot_grid(young, old + theme(strip.text.y = element_blank()), nrow=1,labels=c("(a)","(b)")) %>%
  ggsave("output/figs/Fig.5.pdf",.,width=24,height=12,units="cm")        

