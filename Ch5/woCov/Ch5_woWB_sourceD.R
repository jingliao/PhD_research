###******** sampling distribution of genotypes using Dirichlet model
sampDist <- post_p_ij(data=tableSource, alpha=1,
                      simNumber=simulation)

sampDistST <- sampDist %>% select(Poultry, Other, Ruminants, ST)
source_name <- colnames(sampDistST)

#' Water data, pW_j=P^w(source j) in logit formula
nS_sourcej <- dim(tableSource)[2]
pw_source <- rnorm(n=nS_sourcej-1, mean=0, sd=sigma)

#' delogit the initial parameter for source P, O, R
pW_j <- delogit(par_vector=pw_source)

#' delogit the initial parameter for source P, O, W, R
nH_sourcej <- dim(freqtable)[2]-1
ph_source <- rnorm(n=nH_sourcej-1, mean=0, sd=sigma)
pH_j <- delogit(par_vector=ph_source)
