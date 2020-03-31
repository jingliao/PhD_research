## Functions for updating eta

# pi = eta/sum(eta)
compute_pi <- function(eta) {
  sum_eta = colSums(eta)
  sweep(eta, 2, FUN='/', sum_eta)
}
# same as above, just a single column
compute_pi_j <- function(eta) {
  eta/sum(eta)
}

# log-likelihood l(p; x) for sources
logLH_s <- function(datas, probs){
  sum(datas*log(probs))
}

## log prior ratio for eta ~ gamma(shape=pi_alpha, scale=1)
loggamma_ratio <- function(x1, x2, shape=1, rate=1) {
  sum((shape-1)*log(x1/x2) + (x2-x1)*rate)
}

## log proposal ratio for eta ~ lognormal(mu, sigma)
loglnorm_proposal_ratio <- function(eta_star, eta) {
  sum(log(eta_star/eta))
}

## Function to update eta
update_etaNc <- function(data, state, prior) {

  # We update all the etas one at a time in a random order:
  source_order <- sample(ncol(state$eta))
  type_order <- sample(nrow(state$eta))

  # local accept/reject for this loop
  accept <- reject <- 0

  for (source in source_order) {

    for (type in type_order) {

      # optimise_eta
      # propose eta_star[type] for this source using log normal random walk
      eta_star <- state$eta[,source]
      eta_star[type] <- rlnorm(1, log(state$eta[type, source]), 1)

      # optimise_c
      # propose c_star[type] for this source using normal random walk
      c_star <- state$c[,source]
      c_star[type] <- rnorm(1, mean=state$c[type, source], sd=c_sdprop)  

        
      # update variables that depend on eta
      pi_s_star <- compute_pi_j(eta=eta_star)

      # compute proposal ratio, prior ratio
      ratio_propE <- loglnorm_proposal_ratio(eta_star=eta_star[type], eta=state$eta[type,source])
      ratio_priorE <- loggamma_ratio(x1=eta_star[type], x2=state$eta[type, source], shape=prior$pi_alpha)
      ratio_priorC <- logcauchy_ratio(c_star=c_star[type],
                                      c=state$c[type, source],
                                      c0=0,
                                      scale=c_scale)
        
      # compute source likelihood ratio
      ratio_lh <- logLH_s(data$x[,source], pi_s_star) - logLH_s(data$x[,source], state$pi_s[,source])

      # compute human likelihood ratio
      if (nrow(data$H) > 0) {
        # have humans
        pi_h_star <- state$pi_h
        pi_h_star[,source] <- compute_pi_j(eta=eta_star*exp(c_star))
        lh_h_star <- logLH_h(data$y, data$H, state$theta, pi_h_star)
        ratio_lh_human <- lh_h_star - state$lh_h
        ratio_lh <- ratio_lh + ratio_lh_human
      }
      # compute acceptance, update state if accepted
      logalpha <- ratio_lh + ratio_propE + ratio_priorE + ratio_priorC

      if (log(runif(1, 0, 1)) < min(0, logalpha)) {
        # accept - update eta
        # optimise_eta
        state$eta[,source] <- eta_star
        state$pi_s[,source] <- pi_s_star
        state$c[,source] <- c_star  

        if (nrow(data$H) > 0) {
          state$pi_h <- pi_h_star
          state$lh_h  <- lh_h_star
        }

        accept <- accept+1

      } else {
        # reject - leave state as is
        reject <- reject+1
      }
    }
  }
    state$accept_etaNc <- state$accept_etaNc + c(accept, reject)
  return(state)
}
