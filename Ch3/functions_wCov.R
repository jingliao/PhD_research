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


mcmcalgmCov <- function(qsigma, sigma, Wtheta, Htheta, log_LH,
                        post_pij, design_mtrx){

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
        ph_fmodel_prime <- design_mtrx %*% thetaH_prime

        ph_capF_prime <- as.data.frame(cbind(humanCov$ST,
                                             delogit_mtrx(par_mtrx=ph_fmodel_prime)))

        names(ph_capF_prime) <- c("ST","s1","s2","s3","s4")

        #' update par_pW according to pW_j_prime, P, O, R
        par_pW_prime <- data.frame(ST=SeqTypes,
                                   p_i=rowSums(t(apply(X=post_pij[,-4], MARGIN=1, FUN=function(x){x*pW_j_prime}))))

        #' update par_pH according to par_pW_prime and ph_capF
        par_pH1_M_prime <- merge(x=ph_capF_prime,
                               y=par_pW_prime, all.x=TRUE)

        #' s3: Water
        par_pH1_prime <- with(par_pH1_M_prime,
                              data.frame(ST=ST, pH1=s3*p_i))

        par_pH2_M_prime <- merge(x=ph_capF_prime,
                                 y=post_pij, all.x=TRUE)

        par_pH2_prime <- with(par_pH2_M_prime,
                              data.frame(ST=ST,
                                         prodS1=s1*sampS1,
                                         prodS2=s2*sampS2,
                                         prodS3=s4*sampS3))

        par_pH_CM_prime <- cbind(par_pH1_prime , par_pH2_prime[,-1])

        par_pH_prime <- with(par_pH_CM_prime,
                                 data.frame(
                                     ST=ST,
                                     q_i=pH1+prodS1+prodS2+prodS3))

        #' update joint loglikelihood
        logLH_pW_prime <- sum(tableWater$Water*log(par_pW_prime$p_i))
        logLH_pH_prime <- sum(log(na.omit(par_pH_prime$q_i)))
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
                          post_pij, design_mtrx){

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
        ph_fmodel_prime <- design_mtrx %*% thetaH_prime

        ph_capF_prime <- as.data.frame(cbind(humanCov$ST,
                                             delogit_mtrx(par_mtrx=ph_fmodel_prime)))

        names(ph_capF_prime) <- c("ST","HOther","HPoultry","HWBirds","HWater", "HRuminants")
        
        #' update par_pW according to pW_j_prime, P, O, WB, R
        par_pW_prime <- data.frame(ST=post_pij$ST,
                                   p_i=as.numeric(as.matrix(post_pij[,-5])%*%pW_j_prime))

        #' update par_pH according to par_pW_prime and ph_capF
        par_pH1_M_prime <- merge(x=ph_capF_prime,
                               y=par_pW_prime, all.x=TRUE)

        #' s4: Water
        par_pH1_prime <- with(par_pH1_M_prime,
                              data.frame(ST=ST, prodW=HWater*p_i))

        par_pH2_M_prime <- merge(x=ph_capF_prime,
                                 y=post_pij, all.x=TRUE)

        par_pH2_prime <- with(par_pH2_M_prime,
                              data.frame(ST=ST,
                                         prodO=HOther*Other,
                                         prodP=HPoultry*Poultry,
                                         prodWB=HWBirds*WBirds,
                                         prodR=HRuminants*Ruminants))

        par_pH_CM_prime <- cbind(par_pH1_prime , par_pH2_prime[,-1])

        par_pH_prime <- data.frame(ST=par_pH_CM_prime$ST,
                                   q_i=rowSums(par_pH_CM_prime[,-1]))

        #' update joint loglikelihood
        logLH_pW_prime <- sum(tableWater$Water*log(par_pW_prime$p_i))
        logLH_pH_prime <- sum(log(na.omit(par_pH_prime$q_i)))
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


tidydata <- function(dataname, sourcename, CIpercent, H){

    if(ncol(H)<3){
        
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
    }

    else {
        
        data1 <- dataname %>% group_by(V1, V2) %>%
            do(myhpd(data=.$V3, CIpercent)) %>%
            mutate(Source=sourcename[1])

        data2 <- dataname %>% group_by(V1, V2) %>%
            do(myhpd(data=.$V4, CIpercent)) %>%
            mutate(Source=sourcename[2])

        data3 <- dataname %>% group_by(V1, V2) %>%
            do(myhpd(data=.$V5, CIpercent)) %>%
            mutate(Source=sourcename[3])

        finaldata <- rbind(data1, data2, data3)

        if(length(sourcename)>3){

            data4 <- dataname %>% group_by(V1, V2) %>%
                do(myhpd(data=.$V6, CIpercent)) %>%
                mutate(Source=sourcename[4])

            finaldata <- rbind(data1, data2, data3, data4)

        }
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

