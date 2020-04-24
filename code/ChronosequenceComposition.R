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

#' Read in and subset data
phy.all <- readRDS("output/rds/phy.rds") %>% 
  subset_samples(LSAG) %>%
  sweepOTUs 
phy.erm <- readRDS("output/rds/phy.erm.rds") %>% 
  subset_samples(LSAG) %>%
  sweepOTUs

#' Community dissimilatiry (Jensen-Shannon distance)
jsd.all <- phy.all %>% phyloseq::distance("jsd") %>% sqrt
jsd.erm <- phy.erm %>% phyloseq::distance("jsd") %>% sqrt
  
#' # NMDS
#' Fit nmds and output as a data.frame with sample data
NMDS <- function(dist.in,phy.in,label){
  invisible(capture.output(nmds <-metaMDS(dist.in,trymax=500)))
  print(nmds)
  nmds  %>%
    scores("sites") %>% data.frame %>%
    mutate(Site=sample_data(phy.in)$Site,taxa=label)
}
#+ echo=T, results='hide'
ord.all <- NMDS(jsd.all,phy.all,"All fungi")
ord.erm <- NMDS(jsd.erm,phy.erm,"Putative ErM")

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
    scale_colour_gradient(high = "grey40", low = "grey90",guide=F) +
    geom_text(aes(label=let,x=-Inf,y=Inf),
              data=data.frame(taxa=c("All fungi","Putative ErM"),let=c("(a)","(b)")),
              family="serif",size=4,inherit.aes = F,fontface="bold", vjust=1.8, hjust=-0.45) +
    facet_wrap(~taxa,scales="free",ncol=1)+
    labs(x="NMDS axis 1",y="NMDS axis 2") +
    ggthemes::theme_few() +
    theme(legend.text = element_markdown(),
          legend.position = "bottom"))
ggsave("output/figs/Fig.3.pdf", width=12, height=21, unit="cm") 

#' # perMANOVA
#' All fungi
adonis(jsd.all~Site,data=data.frame(sample_data(phy.all)))
adonis(jsd.all~logAge,data=data.frame(sample_data(phy.all)))
#' Putative ErM
adonis(jsd.erm~Site,data=data.frame(sample_data(phy.erm)))
adonis(jsd.erm~logAge,data=data.frame(sample_data(phy.erm)))

#' # Betadispersion
#' All fungi
beta.all <- betadisper(jsd.all,sample_data(phy.all)$Site)
boxplot(beta.all, main="Betadispersion all fungi")
permutest(beta.all, pairwise = TRUE, permutations = 999)  
#' Putative Erm
beta.erm <- betadisper(jsd.erm,sample_data(phy.erm)$Site)
boxplot(beta.erm, main="Betadispersion putative ErM")
permutest(beta.erm, pairwise = TRUE, permutations = 999)



