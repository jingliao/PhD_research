## Likelihood functions/ratio helpers

## log-likelihood functions ----

## F, for human, source attribution probablty, which is inversed by logit(f(theta))
computeF <- function(theta, H){

  f = H %*% theta
  ef <- exp(cbind(f, 0))
  ef/rowSums(ef)

}

## log-likelihood l(theta, pi_h; y, H) for humans
logLH_h <- function(y, H, theta, pi_h) {
    
    F_star = computeF(theta, H)
    
    LH <- rowSums(F_star * pi_h[y,])
    
    return(sum(log(LH)))    
}

## log prior ratio for c ~ Cauchy(c0=0, scale) 
logcauchy_ratio <- function(c_star, c, c0=0, scale){

  sum(log(1+((c-c0)/scale)^2)-log(1+((c_star-c0)/scale)^2))

}

## log prior ratio for theta, theta^(t) ~ N(0,1)
lognormal_ratio <- function(theta_star, theta, mu=0, sigma=1){

  sum(1/2*1/sigma^2*((theta-mu)^2-(theta_star-mu)^2))

}
