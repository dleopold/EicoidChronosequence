#' ---
#' title: Chronosequence composition
#' author: Devin R Leopold
#' date: April 11, 2020
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: false
#'      self_contained: true
#'      highlight: zenburn
#' ---

#' # Prepare environment
#+ warning=F, message=F
library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)
library(ggvegan)
library(cowplot)
library(patchwork)
library(ggtext)
source("code/Rfunctions.R")

#' Colors for figure
SiteColors <- viridis::viridis_pal(direction=-1,option="inferno")(10)[c(1,2,3,5,8)]

#' Read in, subset data, and make proportional
phy.all <- readRDS("output/rds/phy.rds") %>% 
  subset_samples(LSAG) %>%
  sweepOTUs %>%
  filter_taxa(function(x) {sum(x>0) > 2 }, T) %>%
  transform_sample_counts(function(x){min(sample_sums(.))*x/sum(x)})
phy.erm <- readRDS("output/rds/phy.erm.rds") %>% 
  subset_samples(LSAG) %>%
  sweepOTUs %>%
  filter_taxa(function(x) {sum(x>0) > 2 }, T) %>%
  transform_sample_counts(function(x){min(sample_sums(.))*x/sum(x)})

#' # NMDS
#' Fit nmds and output as a data.frame with sample data
NMDS <- function(phy.in,label){
  dat <- phy.in %>% 
    otu_table %>% data.frame
  nmds <-metaMDS(dat,trymax=500)
  print(nmds)
  nmds  %>%
    scores("sites") %>% data.frame %>%
    mutate(Site=sample_data(phy.in)$Site,taxa=label)
}
#+ echo=T, results='hide'
ord.all <- NMDS(phy.all,"All fungi")
ord.erm <- NMDS(phy.erm,"Putative ErM")

#' # Ordisurf
#' Fit non-linear surface representing site age to ordinations and output as a data.frame for plotting
contours <- function(ord.in,X,label){
  surf <- ordisurf(ord.in[,1:2] ~ X, plot = FALSE)
  print(summary(surf))
  surf.out <- expand.grid(x = surf$grid$x, y = surf$grid$y) %>% 
    mutate(z=as.vector(surf$grid$z)) %>%
    na.omit %>% data.frame %>% mutate(taxa=label)
  return(surf.out)
}
surf.all <- contours(ord.all,sample_data(phy.all)$logAge,"All fungi")
surf.erm <- contours(ord.erm,sample_data(phy.erm)$logAge,"Putative ErM")

#' merge data
ord.dat <- bind_rows(ord.all,ord.erm)
surf.dat <- bind_rows(surf.all,surf.erm)

#' # Plot ordination
#+ fig.align = "center", fig.asp=21/12, dpi=64, out.width="600px"
(ord.plot <- ggplot(ord.dat, aes(x = NMDS1, y = NMDS2)) +
    stat_contour(data=surf.dat, 
                 aes(x=x, y=y, z=z, colour=(..level..)),
                 show.legend=F,size=0.4) +
    geom_point(size = 3.5, shape=21, aes(fill=Site)) + 
    scale_fill_manual(name = "Site age (yr)",values=SiteColors,
                      labels=c("3.0&times;10^2",
                               "2.1&times;10^3",
                               "2.0&times;10^4",
                               "1.5&times;10^5",
                               "4.1&times;10^6")) +
    scale_colour_gradient(high = "grey50", low = "grey90",guide=F) +
    geom_text(aes(label=let,x=-Inf,y=Inf),
              data=data.frame(taxa=c("All fungi","Putative ErM"),let=c("(a)","(b)")),
              family="serif",size=4,inherit.aes = F,fontface="bold", vjust=1.8, hjust=-0.45) +
    facet_wrap(~taxa,scales="free",ncol=1)+
    labs(x="NMDS axis 1",y="NMDS axis 2") +
    ggthemes::theme_few() +
    theme(legend.text = element_markdown(),
          #legend.position = 'bottom',
          strip.text = element_text(size=12),
          axis.title = element_text(size=10),
          axis.text = element_text(size=8)))
ggsave("output/figs/ordinations.pdf", width=11.5, height=16, unit="cm") 

#' # perMANOVA
adonis_phy <- function(phy.in,x){
  X <- phy.in %>% sample_data %>% data.frame %>% select(all_of=x) %>% .[,1]
  Y <- phy.in %>% otu_table %>% data.frame %>% sqrt
  adonis(Y~X, permutations=9999, parallel=parallel::detectCores())
}
#' All fungi
adonis_phy(phy.all,"Site")
adonis_phy(phy.all,"logAge")
#' Putative ErM
adonis_phy(phy.erm,"Site")
adonis_phy(phy.erm,"logAge")

#' # Betadispersion
betadisper_phy <- function(phy.in){
  X <- phy.in %>% sample_data %>% data.frame %$% Site
  dist <- phy.in %>% otu_table %>% data.frame %>% sqrt %>% vegdist 
  betadisper(dist,X)
}
#' All fungi
beta.all <- betadisper_phy(phy.all)
boxplot(beta.all, main="Betadispersion all fungi")
permutest(beta.all, pairwise = TRUE, permutations = 999)  
#' Putative Erm
beta.erm <- betadisper_phy(phy.erm)
boxplot(beta.erm, main="Betadispersion putative ErM")
permutest(beta.erm, pairwise = TRUE, permutations = 999)

