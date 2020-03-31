#####################################################################
##########  Date: 2018.02.22                               ##########
##########  Note: functions for New theory (2nd paper)     ##########
#####################################################################

post_p_ij <- function(data, alpha, simNumber){

    sim_p_ij_list <- list()

    ## posterior dirichlet parameter, alpha_star
    ## no need to extract 1 behind the addition
    ## it is default in Dirichlet commond
    alpha_star <- data+alpha

    ## simulate posterior Dirichlet p(ST_i|source_j)

    for (i in 1:dim(data)[2]) {

        sim_p_ij_list[[i]] <- as.numeric(t(
            rdirichlet(n=simNumber, alpha=alpha_star[,i])))
    }

    ## combine all list(simulations)
    sim_p_ij <- as.data.frame(do.call(cbind, sim_p_ij_list))

    ## rename column names and match ST no. to each p_ij
    colnames(sim_p_ij) <- colnames(data)
    sim_p_ij$ST <- rep(as.integer(rownames(data)), simNumber)

    return(sim_p_ij)
}


delogit <- function(par_vector){

    ## the last one parameter is zero
    full_par <- c(par_vector, 0)  
    
    delogit_full_par <- exp(full_par)/sum(exp(full_par))
    
    return(delogit_full_par)

}

delogit_mtrx <- function(par_mtrx){
   
    full_par <- cbind(par_mtrx, 0) 
    
    delogit_full_par <- exp(full_par)/rowSums(exp(full_par))
    
    return(delogit_full_par)
}


mcmcalgm <- function(qsigma, sigma, thetaW, thetaH, log_LH,
                     post_pij){

    output <- matrix(NA, nrow=Iteration, ncol=(length(thetaW)+
                                               length(thetaH)))

    number <- 1

    for (iMH in 1:totalMH){

        #' update thetaW and thetaH
        thetaW_prime <- thetaW
        thetaH_prime <- thetaH

        #' sample one of thetaW and one of thetaH to update
        id_w <- sample(x=seq(from=1, to=nS_sourcej-1, by=1), size=1)
        id_h <- sample(x=seq(from=1, to=nH_sourcej-1, by=1), size=1)

        thetaW_prime[id_w] <- rnorm(n=1, mean=thetaW[id_w],
                                       sd=qsigma)
        thetaH_prime[id_h] <- rnorm(n=1, mean=thetaH[id_h],
                                       sd=qsigma)

        #' inverse logit
        pW_j_prime <- delogit(par_vector=thetaW_prime)
        pH_j_prime <- delogit(par_vector=thetaH_prime)

        #' update par_pW and par_pH
        par_pW_prime <- as.numeric(as.matrix(post_pij)%*%pW_j_prime)
        
        par_pH1_prime <- par_pW_prime*pH_j_prime[3] #' water
        par_pH2_prime <- as.numeric(as.matrix(post_pij)%*%pH_j_prime[-3])
        
        par_pH_prime <- par_pH1_prime + par_pH2_prime

        #' joint loglikelihood
        logLH_pW_prime <- sum(tableWater$Water*log(par_pW_prime))
        logLH_pH_prime <- sum(tableHum$Human*log(par_pH_prime))
        jointlogLH_prime <- logLH_pW_prime + logLH_pH_prime

        #' prior ratio on MH algorithm
        log_prior_ratio_w <- (thetaW[id_w]*thetaW[id_w]-
                       thetaW_prime[id_w]*thetaW_prime[id_w])/
                       (2*sigma*sigma)
        
        log_prior_ratio_h <- (thetaH[id_h]*thetaH[id_h]-
                       thetaH_prime[id_h]*thetaH_prime[id_h])/
                       (2*sigma*sigma)

        #' log_LH ratio
        log_LH_ratio <- jointlogLH_prime - jointlogLH

        #' alpha, acceptance rate
        log_alpha <- log_LH_ratio +
                     log_prior_ratio_w + log_prior_ratio_h

        log_u <-  log(runif(1))

        if(log_u < min(0, log_alpha)){

            thetaW <- thetaW_prime
            thetaH <- thetaH_prime
            jointlogLH <- jointlogLH_prime

            samplelist <- c(thetaW_prime, thetaH_prime)
            
            accept <- accept + 1

        }

        else {

            reject <- reject + 1
        }


        if( (iMH>Burnin) && (iMH%%Thinning==0) ){

            output[number,] <- samplelist

            number <- number + 1            

        }

    }

    #' print the acceptance rate
    if(iMH%%Thinning==0){

        print(sprintf("Iteration: %d Acceptance rate: %.3f", iMH,
                      accept/iMH))
        print(sprintf("Iteration: %d Rejection rate: %.3f", iMH,
                      reject/iMH))
    }

    return(output)
}

## mcmcalgmI <- function(qsigma, sigma, thetaW, thetaH, log_LH,
##                      post_pij){

##     output <- matrix(NA, nrow=Iteration, ncol=(length(thetaW)+
##                                                length(thetaH)))

##     number <- 1

##     for (iMH in 1:totalMH){

##         #' update thetaW and thetaH
##         thetaW_prime <- thetaW
##         thetaH_prime <- thetaH

##         #' sample one of thetaW and one of thetaH to update
##         id_w <- sample(x=seq(from=1, to=nS_sourcej-1, by=1), size=1)
##         id_h <- sample(x=seq(from=1, to=nH_sourcej-1, by=1), size=1)

##         thetaW_prime[id_w] <- rnorm(n=1, mean=thetaW[id_w],
##                                        sd=qsigma)
##         thetaH_prime[id_h] <- rnorm(n=1, mean=thetaH[id_h],
##                                        sd=qsigma)

##         #' inverse logit
##         pW_j_prime <- delogit(par_vector=thetaW_prime)
##         pH_j_prime <- delogit(par_vector=thetaH_prime)

##         #' update par_pW and par_pH
##         par_pW_prime <- data.frame(ST=SeqTypesI,
##                                    p_i=rowSums(t(apply(X=post_pij,
##                                                        MARGIN=1,
##                                                        FUN=function(x){x*pW_j_prime}))))

## #' human
##         par_pH1_prime <- par_pW_prime$p_i*pH_j_prime[4]
##         par_pH2_prime <- rowSums(t(apply(X=post_pij,
##                                          MARGIN=1,
##                                          FUN=function(x){x*pH_j_prime[-4]})))

##         par_pH_prime <- data.frame(ST=SeqTypesI,
##                                    q_i=par_pH1_prime + par_pH2_prime)


##         #' joint loglikelihood
##         tabsW <- na.omit(join(tableWater, par_pW_prime, type="full"))
##         tabsW$prod <- with(tabsW, Water*log(p_i))
##         logLH_pW_prime <- sum(tabsW$prod)

##         tabsH <- na.omit(join(tableHum, par_pH_prime, type="full"))
##         tabsH$prod <- with(tabsH, Human*log(q_i))
##         logLH_pH_prime <- sum(tabsH$prod)

##         jointlogLH_prime <- logLH_pW_prime + logLH_pH_prime

##         #' prior ratio on MH algorithm
##         log_prior_ratio_w <- (thetaW[id_w]*thetaW[id_w]-
##                        thetaW_prime[id_w]*thetaW_prime[id_w])/
##                        (2*sigma*sigma)
        
##         log_prior_ratio_h <- (thetaH[id_h]*thetaH[id_h]-
##                        thetaH_prime[id_h]*thetaH_prime[id_h])/
##                        (2*sigma*sigma)

##         #' log_LH ratio
##         log_LH_ratio <- jointlogLH_prime - jointlogLH

##         #' alpha, acceptance rate
##         log_alpha <- log_LH_ratio +
##                      log_prior_ratio_w + log_prior_ratio_h

##         log_u <-  log(runif(1))

##         if(log_u < min(0, log_alpha)){

##             thetaW_source <- thetaW_prime
##             thetaH_source <- thetaH_prime
##             jointlogLH <- jointlogLH_prime

##             samplelist <- c(thetaW_prime, thetaH_prime)
            
##             accept <- accept + 1

##         }

##         else {

##             reject <- reject + 1
##         }


##         if( (iMH>Burnin) && (iMH%%Thinning==0) ){

##             output[number,] <- samplelist

##             number <- number + 1            

##         }

##     }

##     #' print the acceptance rate
##     if(iMH%%Thinning==0){

##         print(sprintf("Iteration: %d Acceptance rate: %.3f", iMH,
##                       accept/iMH))
##         print(sprintf("Iteration: %d Rejection rate: %.3f", iMH,
##                       reject/iMH))
##     }

##     return(output)
## }

mcmcalgmWB <- function(qsigma, sigma, thetaW, thetaH, log_LH,
                     post_pij){

    output <- matrix(NA, nrow=Iteration, ncol=(length(thetaW)+
                                               length(thetaH)))

    number <- 1

    for (iMH in 1:totalMH){

        #' update thetaW and thetaH
        thetaW_prime <- thetaW
        thetaH_prime <- thetaH

        #' sample one of thetaW and one of thetaH to update
        id_w <- sample(x=seq(from=1, to=nS_sourcej-1, by=1), size=1)
        id_h <- sample(x=seq(from=1, to=nH_sourcej-1, by=1), size=1)

        thetaW_prime[id_w] <- rnorm(n=1, mean=thetaW[id_w],
                                       sd=qsigma)
        thetaH_prime[id_h] <- rnorm(n=1, mean=thetaH[id_h],
                                       sd=qsigma)

        #' inverse logit
        pW_j_prime <- delogit(par_vector=thetaW_prime)
        pH_j_prime <- delogit(par_vector=thetaH_prime)

        #' update par_pW and par_pH
        par_pW_prime <- as.numeric(as.matrix(post_pij)%*%pW_j_prime)

        par_pH1_prime <- par_pW_prime*pH_j_prime[4]
        par_pH2_prime <- as.numeric(as.matrix(post_pij)%*%pH_j_prime[-4])
        
        par_pH_prime <- par_pH1_prime + par_pH2_prime

        #' joint loglikelihood
        logLH_pW_prime <- sum(tableWater$Water*log(par_pW_prime))
        logLH_pH_prime <- sum(tableHum$Human*log(par_pH_prime))
        jointlogLH_prime <- logLH_pW_prime + logLH_pH_prime

        #' prior ratio on MH algorithm
        log_prior_ratio_w <- (thetaW[id_w]*thetaW[id_w]-
                       thetaW_prime[id_w]*thetaW_prime[id_w])/
                       (2*sigma*sigma)
        
        log_prior_ratio_h <- (thetaH[id_h]*thetaH[id_h]-
                       thetaH_prime[id_h]*thetaH_prime[id_h])/
                       (2*sigma*sigma)

        #' log_LH ratio
        log_LH_ratio <- jointlogLH_prime - jointlogLH

        #' alpha, acceptance rate
        log_alpha <- log_LH_ratio +
                     log_prior_ratio_w + log_prior_ratio_h

        log_u <-  log(runif(1))

        if(log_u < min(0, log_alpha)){

            thetaW <- thetaW_prime
            thetaH <- thetaH_prime
            jointlogLH <- jointlogLH_prime

            samplelist <- c(thetaW_prime, thetaH_prime)
            
            accept <- accept + 1

        }

        else {

            reject <- reject + 1
        }


        if( (iMH>Burnin) && (iMH%%Thinning==0) ){

            output[number,] <- samplelist

            number <- number + 1            

        }

    }

    #' print the acceptance rate
    if(iMH%%Thinning==0){

        print(sprintf("Iteration: %d Acceptance rate: %.3f", iMH,
                      accept/iMH))
        print(sprintf("Iteration: %d Rejection rate: %.3f", iMH,
                      reject/iMH))
    }

    return(output)
}
