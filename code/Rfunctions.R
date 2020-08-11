# remove empty OTUs from phyloseq object
sweepOTUs <- function(phy){
  prune_taxa(taxa_sums(phy) > 0,phy)
}

# mutate over ragged array
mutate_by <- function(.data, group_vars, ...) {
  gvs <- rlang::enquos(group_vars)
  .data %>%
    group_by_at(vars(!!!gvs)) %>%
    mutate(...) %>%
    ungroup
}

# Get estimated diversity and richness estimates from phyloseq object with iNEXT
iNEXT.phy <- function(phy.in,level, use.cores=parallel::detectCores()){
  require(foreach)
  require(doMC)
  registerDoMC(cores=use.cores)
  otuTab <- phy.in %>% otu_table %>% t %>% data.frame 
  out <- foreach(i=1:ncol(otuTab), .combine=bind_rows) %dopar% {
    estimateD(otuTab[,i],level=level) %>% 
      tibble %>%
      mutate(SampID=colnames(otuTab)[i])
  }
  out %>% 
    filter(order!=2) %>%
    transmute(SampID=SampID,
              Metric=case_when(order==0 ~ 'Species richness',
                               order==1 ~ 'Shannon diversity'),
              Estimate=qD,
              LCL=qD.LCL,
              UCL=qD.UCL) %>%
    left_join(sample_data(phy.in) %>% data.frame %>%
                rownames_to_column("SampID"))
}

### Stats helper functions ###

# anova testing differences among sites
aov_chrono <- function(x,y){
  mod <- aov(Estimate~Site, data=x) %>% broom::tidy()
  tibble(Metric=y$Metric,
         Taxa=y$Taxa,
         df=mod$df[1],
         df.resid=mod$df[2],
         Fval=mod$statistic[1],
         Pval=mod$p.value[1])
}

# Tukey HSD wrapper to order letters according to factor levels
tukeyWrap <- function(x,y){
  aovMod <- aov(Estimate~Site, data=x)
  factorLevels <- aovMod$xlevels$Site
  hsd <- TukeyHSD(aovMod)
  Tukey.levels <- hsd[["Site"]][,4]
  Tukey.labels <- data.frame(multcompView::multcompLetters(Tukey.levels)['Letters'],stringsAsFactors = F)
  Tukey.labels %<>% rownames_to_column("Site")
  Tukey.labels %<>% .[match(factorLevels,.$Site) ,]
  Tukey.labels$Site %<>% factor(levels=factorLevels)
  uniqLets <- Tukey.labels$Letters %>% as.character() %>% strsplit("") %>% unlist %>% unique
  for(i in 1:nrow(Tukey.labels)){
    foo <- Tukey.labels$Letters[i] %>% as.character() %>% strsplit("") %>% unlist
    Tukey.labels$Letters[i] <- letters[which(uniqLets %in% foo)] %>% paste0(collapse="")
  }
  tibble(Tukey.labels) %>%
    mutate(Taxa=y$Taxa,
           Metric=y$Metric)
}

# anova of fertilizer effects
aov_fert <- function(x,y){
  mod <- aov(Estimate~N*P, data=x) %>% broom::tidy() 
  tibble(Metric=y$Metric,
         Taxa=y$Taxa,
         Site=y$Site,
         Predictor=mod$term,
         df=mod$df,
         sumsq=mod$sumsq,
         meansq=mod$meansq,
         Fval=mod$statistic,
         Pval=mod$p.value)
}

# bootstrap CIs for plotting
bootFun <- function(x,y,n=10000){
  boots <- replicate(n, mean(sample(x$Estimate, replace=T)))
  tibble(Metric=y$Metric,
         Site=y$Site,
         Taxa=y$Taxa,
         Estimate=boots)
}
bootFun_fert <- function(x,y,n=10000){
  boots <- replicate(n, mean(sample(x$Estimate, replace=T)))
  tibble(Metric=y$Metric,
         Site=y$Site,
         Taxa=y$Taxa,
         Predictor=y$Frt_treat,
         Estimate=boots)
}

