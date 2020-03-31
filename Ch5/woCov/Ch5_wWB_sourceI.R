###******** sampling distribution of genotypes using island model
islandtypes <- alltypes
islandtypes$Source <- fct_collapse(islandtypes$Source,
                                   HumanWater = c("Human", "Water"))

st2=st_fit(formula=Source~ST, non_primary="HumanWater",
           method="island", data=islandtypes,
           sequences =~ASP+GLN+GLT+GLY+PGM+TKT+UNC,
           samples=simulation)

sampDistlist <- list()

for(s in 1:simulation){

    sampDistSingle <- as.data.frame(st2$sampling_distribution[,,s])
    sampDistSingle$ST <- rownames(sampDistSingle)
    sampDistlist[[s]] <- sampDistSingle[order(as.numeric(sampDistSingle$ST)),]

}

sampDist <- do.call(rbind, sampDistlist)
#' Poultry, Other, WBird, Ruminants(b/l)
sampDistST <- sampDist %>% select(Poultry, Other, WildWaterBird, Ruminants, ST)

#' Water data, pW_j=P^w(source j) in logit formula
#' sources: Poultry, Other, Wbirds, Ruminants
nS_sourcej <- dim(tableSource)[2]
pw_source <- rnorm(n=nS_sourcej-1, mean=0, sd=sigma)

#' delogit the initial parameter for source
#' Poultry, Other, Wbirds
pW_j <- delogit(par_vector=pw_source)

#' delogit the initial parameter for source Poultry, Other,
#' Wbirds, Water, Ruminants
nH_sourcej <- dim(freqtable)[2]-1
ph_source <- rnorm(n=nH_sourcej-1, mean=0, sd=sigma)

pH_j <- delogit(par_vector=ph_source)

