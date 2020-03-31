#####################################################################
##########  Date: 2020.02.06                               ##########
##########  Note: functions for New theory                 ##########
#####################################################################

computeF <- function(theta, H){

    f = H %*% theta
    ef <- exp(cbind(f, 0))

    return(ef/rowSums(ef))

}

## just use matrix multiplier =_=

## computeP_i <- function(pi, G){

##     prod <- pi %>% apply(., MARGIN=1, function(x) x*G)
##     sum <- as.numeric(colSums(prod))
##     prod <- pi %*% G
##     return(sum)    
## }

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

delogit_mtrx <- function(par_mtrx){
   
    full_par <- cbind(par_mtrx, 0) 
    
    delogit_full_par <- exp(full_par)/rowSums(exp(full_par))
    
    return(delogit_full_par)
}

delogit <- function(par_vector){

    ## the last one parameter is zero
    full_par <- c(par_vector, 0)  
    
    delogit_full_par <- exp(full_par)/sum(exp(full_par))
    
    return(delogit_full_par)

}


sample_theta <- function(para_rowno){ 

    matrix(rbind(rnorm(nH_sourcej-1,0,1),rnorm(nH_sourcej-1,0,1)),
           nrow=para_rowno)
}


mcmcalgmCov <- function(qsigma, sigma, Wtheta, Htheta, log_LH,
                        post_pij, design_mtrx, y){

    number <- 1

    for (iMH in 1:totalMH){

        #' update thetaW and thetaH
        thetaW_prime <- Wtheta #' P, O (B/L:Rum)
        thetaH_prime <- Htheta #' P, O, W (B/L:Rum)

        #' sample one of thetaW and one of thetaH to update
        id_w <- sample(x=seq(from=1, to=length(Wtheta), by=1),
                       size=1)
        id_h <- sample(x=seq(from=1, to=length(Htheta), by=1),
                       size=1)

        thetaW_prime[id_w] <- rnorm(n=1, mean=Wtheta[id_w],
                                       sd=qsigma)
        thetaH_prime[id_h] <- rnorm(n=1, mean=Htheta[id_h],
                                       sd=qsigma)

        #' inverse logit
        ## water, P, O, R
        pW_j_prime <- delogit(par_vector=thetaW_prime)

        ## human, Rurality model, P, O, W, R(B/L)
        ph_capF_prime <- computeF(theta=thetaH_prime, H=design_mtrx)
        
        #' update par_pW according to pW_j_prime, P, O, R
        par_pW_prime <- post_pij %*% pW_j_prime

        #' update par_pH according to par_pW_prime and ph_capF
        #' s3: F_Water
        comb_pH_prime <- cbind(par_pW_prime, post_pij)

        par_pH_prime <- as.numeric(rowSums(comb_pH_prime[y,]*ph_capF_prime[,c(3,1,2,4)]))
                 
        #' update joint loglikelihood
        logLH_pW_prime <- sum(tableWater$Water*log(par_pW_prime))
        logLH_pH_prime <- sum(log(par_pH_prime))
        jointlogLH_prime <- logLH_pW_prime + logLH_pH_prime

        #' prior ratio on MH algorithm
        log_prior_ratio_w <- (Wtheta[id_w]*Wtheta[id_w]-
                       thetaW_prime[id_w]*thetaW_prime[id_w])/
                       (2*sigma*sigma)
        
        log_prior_ratio_h <- (Htheta[id_h]*Htheta[id_h]-
                       thetaH_prime[id_h]*thetaH_prime[id_h])/
                       (2*sigma*sigma)

        #' log_LH ratio
        log_LH_ratio <- jointlogLH_prime - jointlogLH

        #' alpha, acceptance rate
        log_alpha <- log_LH_ratio +
                     log_prior_ratio_w + log_prior_ratio_h

        log_u <-  log(runif(1))

        if(log_u < min(0, log_alpha)){

            Wtheta <- thetaW_prime
            Htheta <- thetaH_prime
            jointlogLH <- jointlogLH_prime

            samplelist <- c(jointlogLH_prime,  logLH_pW_prime,
                            logLH_pH_prime, thetaW_prime,
                            t(thetaH_prime))
            
            accept <- accept + 1

        }

        else {

            reject <- reject + 1
        }


        if( (iMH>Burnin) && (iMH%%Thinning==0) ){

            output[[number]] <- list(samplelist,
                                      par_pW_prime,
                                      par_pH_prime)

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

cidataCov <- function(data, colSource, scaleRurality, CIpercent){
   
    s <- data[,c(1,colSource)]
    means <- mean(s[s$V1==scaleRurality,2])
    Cis <- HPDinterval(mcmc(s[s$V1==scaleRurality,2]), CIpercent/100)
    listci <- list(mean=means,
                   lower=Cis[1],
                   upper=Cis[2])

    return(listci)
}



cidataCovMed <- function(data, colSource, scaleRurality, CIpercent){
   
    s <- data[,c(1,colSource)]
    means <- median(s[s$V1==scaleRurality,2])
    Cis <- HPDinterval(mcmc(s[s$V1==scaleRurality,2]), CIpercent/100)
    listci <- list(mean=means,
                   lower=Cis[1],
                   upper=Cis[2])

    return(listci)
}


col_alpha <- function(col, alpha){
    rgb(t(col2rgb(col)/255), alpha=alpha)}



mcmcalgmCovWB <- function(qsigma, sigma, Wtheta, Htheta, log_LH,
                          post_pij, design_mtrx, y){

    output <- matrix(NA, nrow=Iteration, ncol=(3+length(Wtheta)+
                                               length(Htheta)))

    number <- 1

    for (iMH in 1:totalMH){

        #' update thetaW and thetaH
        thetaW_prime <- Wtheta #' P, O, WB, R(B/L)
        thetaH_prime <- Htheta #' P, O, WB, W, R(B/L)

        #' sample one of thetaW and one of thetaH to update
        id_w <- sample(x=seq(from=1, to=length(Wtheta), by=1),
                       size=1)
        id_h <- sample(x=seq(from=1, to=length(Htheta), by=1),
                       size=1)

        thetaW_prime[id_w] <- rnorm(n=1, mean=Wtheta[id_w],
                                       sd=qsigma)
        thetaH_prime[id_h] <- rnorm(n=1, mean=Htheta[id_h],
                                       sd=qsigma)

        #' inverse logit
        ## water, P, O, R
        pW_j_prime <- delogit(par_vector=thetaW_prime)

        ## human, Rurality model, P, O, W, R(B/L)
        ph_capF_prime <- computeF(theta=thetaH_prime, H=design_mtrx)
        
        #' update par_pW according to pW_j_prime, P, O, WB, R
        par_pW_prime <- post_pij %*% pW_j_prime
                  
        #' update par_pH according to par_pW_prime and ph_capF
        #' s4: F_Water
        comb_pH_prime <- cbind(par_pW_prime, post_pij)
        
        par_pH_prime <- as.numeric(rowSums(
            comb_pH_prime[y,]*ph_capF_prime[,c(4,1,2,3,5)]))

        #' update joint loglikelihood
        logLH_pW_prime <- sum(tableWater$Water*log(par_pW_prime))
        logLH_pH_prime <- sum(log(par_pH_prime))
        jointlogLH_prime <- logLH_pW_prime + logLH_pH_prime

        #' prior ratio on MH algorithm
        log_prior_ratio_w <- (Wtheta[id_w]*Wtheta[id_w]-
                       thetaW_prime[id_w]*thetaW_prime[id_w])/
                       (2*sigma*sigma)
        
        log_prior_ratio_h <- (Htheta[id_h]*Htheta[id_h]-
                       thetaH_prime[id_h]*thetaH_prime[id_h])/
                       (2*sigma*sigma)

        #' log_LH ratio
        log_LH_ratio <- jointlogLH_prime - jointlogLH

        #' alpha, acceptance rate
        log_alpha <- log_LH_ratio +
                     log_prior_ratio_w + log_prior_ratio_h

        log_u <-  log(runif(1))

        if(log_u < min(0, log_alpha)){

            Wtheta <- thetaW_prime
            Htheta <- thetaH_prime
            jointlogLH <- jointlogLH_prime

            samplelist <- c(jointlogLH_prime,  logLH_pW_prime,
                            logLH_pH_prime, thetaW_prime,
                            t(thetaH_prime))
            
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


fun_Fj <- function(theta, design_mtrx){

    f <- design_mtrx %*% theta
    F <- as.data.frame(na.omit(delogit_mtrx(par_mtrx=f)))
   return(F)
}


myhpd <- function(data, CIpercent){
    
    means <- mean(data)
    
    orderdata <- data[order(data)]
    
    ncis <- round(length(orderdata)*(1-CIpercent/100))
    size <- numeric(ncis)
    
    for(i in 1:ncis){
        
        lower <- orderdata[i]
        upper <- orderdata[length(orderdata)-ncis+i]
        size[i] <- upper-lower
        
    }
    
    sizemin <- which.min(size)
    
    dataframe <- data.frame(mean=means, 
                   lower=orderdata[sizemin], 
                   upper=orderdata[length(orderdata)-ncis+sizemin])

    return(dataframe)
}


tidydata <- function(dataname, sourcename, CIpercent){

      
        data1 <- dataname %>% group_by(V1) %>%
            do(myhpd(data=.$V2, CIpercent)) %>%
            mutate(Source=sourcename[1])

        data2 <- dataname %>% group_by(V1) %>%
            do(myhpd(data=.$V3, CIpercent)) %>%
            mutate(Source=sourcename[2])

        data3 <- dataname %>% group_by(V1) %>%
            do(myhpd(data=.$V4, CIpercent)) %>%
            mutate(Source=sourcename[3])

        finaldata <- rbind(data1, data2, data3)

        if(length(sourcename)>3){

            data4 <- dataname %>% group_by(V1) %>%
                do(myhpd(data=.$V5, CIpercent)) %>%
                mutate(Source=sourcename[4])

            finaldata <- rbind(data1, data2, data3, data4)

        }
  
    return(finaldata)

}


tidydataRA <- function(dataname, sourcename, CIpercent){

    data1 <- dataname %>% group_by(V1, Age) %>%
        do(myhpd(data=.$V2, CIpercent)) %>%
        mutate(Source=sourcename[1])

    data2 <- dataname %>% group_by(V1, Age) %>%
        do(myhpd(data=.$V3, CIpercent)) %>%
        mutate(Source=sourcename[2])

    data3 <- dataname %>% group_by(V1, Age) %>%
        do(myhpd(data=.$V4, CIpercent)) %>%
        mutate(Source=sourcename[3])

    finaldata <- rbind(data1, data2, data3)

    if(length(sourcename)>3){

        data4 <- dataname %>% group_by(V1, Age) %>%
            do(myhpd(data=.$V5, CIpercent)) %>%
            mutate(Source=sourcename[4])

        finaldata <- rbind(data1, data2, data3, data4)

    }

    return(finaldata)

}



tidydataWB <- function(dataname, sourcename, CIpercent){

    data1 <- dataname %>% group_by(V1) %>%
        do(myhpd(data=.$V2, CIpercent)) %>%
        mutate(Source=sourcename[1])

    data2 <- dataname %>% group_by(V1) %>%
        do(myhpd(data=.$V3, CIpercent)) %>%
        mutate(Source=sourcename[2])

    data3 <- dataname %>% group_by(V1) %>%
        do(myhpd(data=.$V4, CIpercent)) %>%
        mutate(Source=sourcename[3])

    data4 <- dataname %>% group_by(V1) %>%
        do(myhpd(data=.$V5, CIpercent)) %>%
        mutate(Source=sourcename[4])

    finaldata <- rbind(data1, data2, data3, data4)

    if(length(sourcename)>4){

        data5 <- dataname %>% group_by(V1) %>%
            do(myhpd(data=.$V6, CIpercent)) %>%
            mutate(Source=sourcename[5])

        finaldata <- rbind(data1, data2, data3, data4, data5)

    }

    return(finaldata)

}

