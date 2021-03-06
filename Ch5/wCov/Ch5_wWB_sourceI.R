###******** load the data
alltypes <- read.csv("human_types_wbirds.csv")

###******** divide the data into human, source(P,R,O) and water
#' frequency table of ST vs. sources
freqtable <- as.data.frame.matrix(table(alltypes$ST,
                                          alltypes$Source))

SeqTypes <- rownames(freqtable)

#' human data, h_{i}
tableHum <- freqtable %>% select(Human)

#' source data, x_{ij}
tableSource <- freqtable %>% select(Other:Ruminants, WildWaterBird)

#' water data
tableWater <- freqtable %>% select(Water)

###******** starting points of parameters
#' Water data
thetaW <- numeric()

nS_sourcej <- dim(tableSource)[2]


for (j in 1:(nS_sourcej-1)){

    thetaW[j] <- rnorm(n=1, mean=0, sd=sigma)

}

#' delogit thetaW
pW_j <- delogit(par_vector=thetaW)

#' for human data
nH_sourcej <- dim(freqtable)[2]-1

humantypes <- alltypes %>% filter(Source=="Human")

humanCov <- humantypes %>% select(ST, UR2006_num) %>%
    rename("Rurality"=UR2006_num) %>% na.omit

humanST <- as.numeric(humanCov$ST)

#' design matrix for covariate Rurality (x1)
designR_mtrx <- model.matrix(ST~Rurality, data=humanCov)

thetaH_R <- sample_theta(para_rowno=ncol(designR_mtrx))

#' compute F
ph_capF_R <- computeF(theta=thetaH_R, H=designR_mtrx)

####***** choose which model for simulating sampling distribution
#' island model

islandtypes <- alltypes
islandtypes$Source <- fct_collapse(islandtypes$Source,
                                   HumanWater = c("Human", "Water"))

st2=st_fit(formula=Source~ST, non_primary="HumanWater",
           method="island", data=islandtypes,
           sequences = ~ ASP+GLN+GLT+GLY+PGM+TKT+UNC,
           samples=simulation)

sampDistlist <- list()

for(s in 1:simulation){

    sampDistSingle <- as.data.frame(st2$sampling_distribution[,,s])
    sampDistSingle$ST <- rownames(sampDistSingle)
    sampDistlist[[s]] <- sampDistSingle[order(as.numeric(sampDistSingle$ST)),]

}

sampDist <- do.call(rbind, sampDistlist)
sourceST <- as.numeric(unique(sampDist$ST))

sampDistST <- sampDist %>%
    select(Other, Poultry, WildWaterBird, Ruminants) %>%
    rename("WBirds"=WildWaterBird)
