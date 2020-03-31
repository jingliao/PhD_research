##**** required library ----
library(dplyr)
library(purrr) # pull out theta or elements from lists

##**** load data, specify y and x ----
base_path <- "data"
alltypes <- read.csv(file.path(base_path, "human_types_age.csv"))

## baseline: ruminants
freqtable <- as.data.frame.matrix(table(alltypes$ST, alltypes$Source)) %>% select(Human, Other, Poultry, Water, Ruminants)

## x, number of ST typed from source
x <- as.matrix(freqtable %>% select(-Human))

source_name <- colnames(x)
nsources <- dim(x)[2]
ntypes <- dim(x)[1]

## human data with covariates and ST; y[h] represents the row in x/pi_s/pi_h that represents the type that human [h] has.
## i.e. rownames(x)[y[h]] = humandata$ST[h].
humandata <- alltypes %>% filter(Source=="Human") %>% select(ST, UR_num) %>% na.omit

## human data
y <- match(humandata$ST, rownames(x))

# check we got this right
all.equal(rownames(x)[y], as.character(humandata$ST))

#### w/o covariate
H  <- model.matrix(ST ~ 1, data=humandata)

#### w/ covariate, Rurality
H  <- model.matrix(ST ~ UR_num, data=humandata)

################################ different setting before running
################ choose setting of c
IQR=0.01
c_scale=IQR/2

#### c=0
## to test if the model is valid as it would be same as the full bayesian model
source("update_eta.R")
c = matrix(0, nrow=ntypes, ncol=nsources)
source("sa_full_bayes_eta.R")


#### random c
c_sdprop=0.05
source("update_etaNc.R")

#' option A. c is random without seed
c = matrix(rcauchy(ntypes*nsources, 0, scale=c_scale), nrow=ntypes, ncol=nsources)

#' option B. with seed
set.seed(123)
c = matrix(rcauchy(ntypes*nsources, 0, scale=c_scale), nrow=ntypes, ncol=nsources)

set.seed(321)
c = matrix(rcauchy(ntypes*nsources, 0, scale=c_scale), nrow=ntypes, ncol=nsources)

source("sa_full_bayes_etaNc.R")

## **** MCMC algorithm ----
n_iter <- 5000
n_thin <- 10
n_burn <- 2000 

posterior <- postThin <- list()

## Main MCMC 

for (i in -n_burn:n_iter){

    for (j in 1:n_thin){

# update our current state

################################ choose c setting                 
        state <- update_eta(data, state, prior)   # for c=0

        state <- update_etaNc(data, state, prior) # for random c
################################
        
        state <- update_theta(data, state, prior)
        
        postThin[[j]] <- state

    }
    
    if(i>0){
        
        posterior[[i]] <- postThin[[n_thin]] # Now sample our posterior

    }
}


