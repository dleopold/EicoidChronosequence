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
div.dat <- bind_rows(iNEXT.phy(phy.all, level=2000) %>% mutate(Taxa="All fungi"),
                     iNEXT.phy(phy.erm, level=2000) %>% mutate(Taxa="Putative ErM")) %>%
  mutate(Frt_treat=factor(Frt_treat,levels=c("C","N","P","NP")),
         Metric=factor(Metric,levels=c("Species richness","Shannon diversity"))) 

#' # Test fertilization effects 
(div.stat <- div.dat %>%
  group_by(Site,Metric,Taxa) %>%
  group_map(~aov_fert(.x,.y)) %>%
  bind_rows) %>%
  gt() %>%
  fmt_number(vars(sumsq, meansq, Fval, Pval),
             decimals = 3)
  
#' Bootstrap distributions for plotting uncertainty 
boot.dat <- div.dat %>%
  select(Metric,Site,Taxa, Frt_treat, Estimate) %>%
  group_by(Metric,Taxa, Frt_treat, Site) %>%
  group_map(~bootFun_fert(.x,.y)) %>%
  bind_rows %>%
  group_by(Metric,Taxa, Predictor, Site) %>%
  summarize(mean=mean(Estimate),
            LCI=coxed::bca(Estimate)[1],
            UCI=coxed::bca(Estimate)[2]) %>%
  mutate(Predictor=factor(Predictor,levels=levels(div.dat$Frt_treat)))

#' # Plot
#' We will make separate plots for the 2 fertilized experiments and then combine them
plotter <- function(site){
  dat <- filter(div.dat, Site==site)
  means <- filter(boot.dat, Site==site)
  title <- ifelse(site=="Thurston","Young (300 yr)","Old (4.1 myr)")
  stats <- filter(div.stat, Site==site) %>%
    mutate(stars=gtools::stars.pval(Pval) %>%
             gsub("*","&#42;",.,fixed=T) %>%
             gsub(".","^+ ",.,fixed=T),
           Predictor=gsub(":","*&times;*",Predictor,fixed=T),
           label=paste0("*",Predictor,"* ", stars)) %>%
    group_by(Metric,Taxa) %>%
    summarize(label=paste(label,collapse = "<br>"))
  ggplot(dat,aes(x=Frt_treat,y=Estimate))+
    geom_boxplot(color="grey65")+
    #ggbeeswarm::geom_quasirandom(size=1.5,position=position_jitter(width = 0.2),
    #                             shape=19,aes(color=seScale), show.legend = F, alpha=0.5)+
    #scale_color_gradient(low="grey20",high="grey80")+
    geom_pointrange(data=means, shape=21, fill="white", size=1, fatten=2,
                    aes(ymin=LCI,ymax=UCI,y=mean,x=Predictor)) +
    geom_richtext(data=stats,aes(x=4.6,y=0,label=label),
                  vjust=0,hjust=-0,fill = NA, label.color = NA)+
    coord_cartesian(clip = 'off') +
    facet_grid(Metric~Taxa,scales="free",
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
  ggsave("output/figs/fertDivPlot.pdf",.,width=24,height=12,units="cm")        

#' # Composition
library(vegan)
adonis_fert <- function(phy.in){
  X <- phy.in %>% sample_data %>% data.frame 
  Y <- phy.in %>% 
    filter_taxa(function(x) {sum(x>0) > 2 }, T) %>%
    transform_sample_counts(function(x){min(sample_sums(.))*x/sum(x)}) %>%
    otu_table %>% data.frame %>% sqrt
  print(adonis(Y~N*P, data=X, permutations = 9999, parallel=parallel::detectCores()))
  adonis2(Y~N*P, data=X, permutations = 9999, parallel=parallel::detectCores(), by=NULL)
}
phy.all %>% subset_samples(Site=="Thurston") %>%
  adonis_fert()
phy.all %>% subset_samples(Site=="Kokee") %>%
  adonis_fert()
phy.erm %>% subset_samples(Site=="Thurston") %>%
  adonis_fert()
phy.erm %>% subset_samples(Site=="Kokee") %>%
  adonis_fert()

#' Ordination
NMDS_fert <- function(phy.in){
  otuTab <- phy.in %>% 
    filter_taxa(function(x) {sum(x>0) > 2 }, T) %>% 
    transform_sample_counts(function(x){min(sample_sums(.))*x/sum(x)}) %>%
    otu_table %>% data.frame
  nmds <-metaMDS(otuTab,trymax=1000,parallel=parallel::detectCores())
  print(nmds)
  nmds  %>%
    scores("sites") %>% data.frame %>%
    mutate(Treatment=sample_data(phy.in)$Frt_treat,
           NMDS1=scale(NMDS1),
           NMDS2=scale(NMDS2)) %>%
    group_by(Treatment) %>%
    summarise(n=n(),
              NMDS1_se=sd(NMDS1)/sqrt(n),
              NMDS1=mean(NMDS1),
              NMDS2_se=sd(NMDS2)/sqrt(n),
              NMDS2=mean(NMDS2))
}
NMDS_all_thurston <- subset_samples(phy.all, Site=="Thurston") %>% NMDS_fert()
NMDS_all_kokee <- subset_samples(phy.all, Site=="Kokee") %>% NMDS_fert()
NMDS_erm_thurston <- subset_samples(phy.erm, Site=="Thurston") %>% NMDS_fert()
NMDS_erm_kokee <- subset_samples(phy.erm, Site=="Kokee") %>% NMDS_fert()

bind_rows(NMDS_all_thurston %>% mutate(Site="Young (300 yr)", Taxa="All Fungi"),
          NMDS_all_kokee %>% mutate(Site="Old (4.1 myr)", Taxa="All Fungi"),
          NMDS_erm_thurston %>% mutate(Site="Young (300 yr)", Taxa="Putative ErMF"),
          NMDS_all_thurston %>% mutate(Site="Old (4.1 myr)", Taxa="Putative ErMF")) %>%
  mutate(Site=factor(Site, levels=c("Young (300 yr)","Old (4.1 myr)"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=Treatment)) +
  geom_errorbarh(aes(xmin=NMDS1-NMDS1_se, xmax=NMDS1+NMDS1_se), height=0)+
  geom_errorbar(aes(ymin=NMDS2-NMDS2_se, ymax=NMDS2+NMDS2_se), width=0)+
  geom_point(size=3) +
  scale_color_manual(name = "Fertilizer treatment",
                    labels=c("Control","Nitrogen","Phosphorus","Both"),
                    values=c("grey80","#1b9e77","#d95f02","#756bb1")) +
  facet_grid(Taxa~Site,scales="free") +
  ggthemes::theme_few() +
  theme(axis.text=element_blank(),
        axis.ticks = element_blank())
ggsave("output/figs/fertCompPlot.pdf",width=15.5,height=10,units="cm")
