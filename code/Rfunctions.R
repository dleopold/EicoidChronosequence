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

# Get asympotic diversity and richness estimated from phyloseq object with iNEXT
iNEXT.phy <- function(phy.in){
  phy.in %>% otu_table %>% data.frame %>% t %>%
    iNEXT(conf=0.95) %$% AsyEst %>% rename(SampID=Site) %>%
    mutate(SampID=rep(sample_names(phy.in),each=3)) %>%
    left_join(sample_data(phy.in) %>% data.frame %>%
                rownames_to_column("SampID"))
}

#############################################################################
### Wrapper functions for the betta() function form the breakaway package ###
#############################################################################

# Basic wrapper to run betta() and extract the results table
betta_wrap <- function(x){
  betta(chats=x$Estimator,
        ses=x$s.e.,
        X=model.matrix(~x$N*x$P)) %$% 
    table %>% data.frame %>% 
    rownames_to_column("Predictor") %>%
    mutate(Predictor=gsub("x$","",Predictor,fixed = T),
           Predictor=gsub("TRUE","",Predictor,fixed = T))
}

# Perform a global test of a multi-lvel factor using a likelihood ratio test and AICc
betta_LRT <- function(x){
  m0 <- betta(chats=x$Estimator,ses=x$s.e.)
  m1 <- betta(chats=x$Estimator,ses=x$s.e.,
              X=model.matrix(~x$Site))
  LR <- -2*(m0$loglikelihood-m1$loglikelihood)
  df <- nrow(m1$table)-nrow(m0$table)
  p <- pchisq(LR,df=df,lower.tail=FALSE)
  deltaAICc <- m0$aicc - m1$aicc
  data.frame(deltaAICc=deltaAICc,
             LR=LR,
             df=df,
             pval=p)
}

#' Test for pairwise differences in diversity accoutning for uncertainty in estimates
betta_pairwise <- function(x){
  # Pair-wise post-hoc test (https://github.com/adw96/breakaway/issues/84#issuecomment-613153006)
  div.test <- betta(chats=x$Estimator,
                    ses=x$s.e.,
                    X=model.matrix(~x$Site-1))
  A <- matrix(c(1, -1, 0, 0, 0, 1, 0, -1, 0, 0, 1, 0, 0, -1, 0, 1, 0, 0, 0, -1,
                0, 1, -1, 0, 0, 0, 1, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, 1, -1, 0,
                0, 0, 1, 0, -1, 0, 0, 0, 1, -1), ncol = 5, byrow = TRUE)
  stats <- MASS::ginv(expm::sqrtm(A %*% div.test$cov %*% t(A))) %*% 
    A %*% div.test$table[,'Estimates']
  pvals <- 2*(1 - stats %>% abs %>% pnorm) %>% p.adjust("fdr")
  names(pvals) <- apply(A,1,function(y){paste(levels(x$Site)[which(y==1)],
                                              levels(x$Site)[which(y==-1)],sep="-")})
  multcompView::multcompLetters(pvals) %$% Letters %>% 
    data.frame(stringsAsFactors = F) %>% rownames_to_column("Site") %>% rename(Let='.')
}

# Extract predicted means and CIs for plotting
betta_means <- function(x,X){
  mod.mx <- model.matrix(as.formula(X),data=x)
  betta(x$Estimator, x$s.e., X=mod.mx) %$% table %>%
    data.frame %>%
    rownames_to_column("Predictor") %>%
    transmute(Predictor=Predictor,
              mean=Estimates,
              LCI=mean-(1.95*Standard.Errors),
              UCI=mean+(1.95*Standard.Errors))
}
