## Functions for updating eta

update_theta <- function(data, state, prior) {

  # We update all the thetas one at a time in a random order:
  theta_order <- sample(length(state$theta))

  # local accept/reject for this loop
  accept <- reject <- 0

  for (i in theta_order) {

    # update theta[i] using a normal random walk
    theta_star = state$theta
    theta_star[i] <- rnorm(1, mean=state$theta[i], sd=1)

    ## Compute LogLik for theta
    ratio_priorT <- lognormal_ratio(theta_star=theta_star, theta=state$theta, mu=prior$mu_theta, sigma=prior$sigma_theta)

    # compute human likelihood ratio, only update theta
    lh_h_star <- logLH_h(data$y, data$H, theta_star, state$pi_h)
    ratio_lh_human <- lh_h_star - state$lh_h

    # compute acceptance, update state if accepted
    logalpha <- ratio_lh_human + ratio_priorT

    if (log(runif(1, 0, 1)) < min(0, logalpha)) {
      # accept - update theta
      state$theta <- theta_star
      state$lh_h  <- lh_h_star

      accept <- accept+1
      
    } else {
      
      reject <- reject + 1
      
    }
  }

  state$accept_theta <- state$accept_theta+c(accept, reject)
  
  return(state)
  
}
