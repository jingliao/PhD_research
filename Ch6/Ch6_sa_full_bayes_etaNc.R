## load sources
source("Ch6_likelihoods.R")
source("Ch6_update_theta.R")

##**** initial parameters ----
## hyper priors ----
## prior on source prob, which is same as prior on eta ~ gamma(alpha, beta=1) 
pi_alpha = 1

## sample parameters ----
## eta, for source, sample it then transform to pi_s with function compute_pi_s
eta=matrix(rgamma(ntypes*nsources, pi_alpha), nrow=ntypes, ncol=nsources)
# name row/colunms to ensure they are the same as source data
rownames(eta) <- rownames(x)
colnames(eta) <- colnames(x)
# compute source probability
pi_s = compute_pi(eta) 

# compute pi_h for a matrix
pi_h <- compute_pi(eta=eta*exp(c))

## theta, for human, parameters in modeling F on logit scale; para_rowno indicates number of parameters
sample_theta <- function(para_rowno){ 
  
  matrix(rep(rnorm(nsources-1,0,1), para_rowno), nrow=para_rowno)
  
}

# if only consider one covariate, rurality
theta <- sample_theta(para_rowno=ncol(H))
# name row and column for theta, it should be same as 
rownames(theta) <- colnames(H)
# NOTE: assumes the 4th source is the default (this is consistent with computeF)
colnames(theta) <- colnames(x)[1:(nsources-1)] 

## reject/accept in the loop
accept <- reject <- 0
accept_theta <- c(0,0)
accept_etaNc <- c(0,0)

################ name rows and columns to ensure they are same as source data
rownames(c) <- rownames(x)
colnames(c) <- colnames(x)

## State object for the MCMC. The idea is we store everything we need to know in a list, and then just update that accordingly
state <- list(eta=eta, 
              pi_s=pi_s, # result from eta
              c=c, 
              pi_h=pi_h, # result from pis and c
              theta=theta, 
              lh_h=logLH_h(y, H, theta, pi_h), # Stored likelihood as this is reasonably slow to compute
              accept_etaNc=accept_etaNc,
              accept_theta=accept_theta)

## Date object for the MCMC. This makes it easier to pass lots of stuff around
data <- list(x=x, # Animal counts,
             y=y, # humandata, Human ST lookup table
             H=H  # design matrix
             )

## Prior object for the MCMC. Makes it easier to pass parameters around
prior <- list(pi_alpha=pi_alpha,        # Dirichlet prior on pi_s, which is also Gamma prior on eta, given beta_eta=1
              mu_theta=0,        # Normal prior mu on theta
              sigma_theta=1,  # Normal prior sigma on theta
              c_scale=c_scale           # Cauchy prior on c
              )

