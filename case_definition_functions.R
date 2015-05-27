## PROGRAM: functions for malaria case definition 
## PURPOSE: The functions are used in the script 'case_definition'

## WRITTEN BY: Amos Thairu
## DATE: June 2014

#function to determine the log likelihoods from a logistic regression model at different values of tau in the power fxn
tauloglik <- function(data, taus=seq(0.1,1,0.1)){ 
  result <- NULL
  for (i in 1:length(taus)){
    tau <- taus[i]
    pfmcltau <- data$pfmcl^tau
    glmfit <- glm2(data$febrile ~ pfmcltau, family=binomial(link="logit"))
    result[[i]] <- data.frame(tau=tau,loglik=logLik(glmfit))
  }
  result <- adply(result,1)[,2:3]
  result <- list(result=result,maxloglik=result[result$loglik==max(result$loglik),])
}

#function to find the tau (power function) that yields the largest log likelihoods
bestTau.fit <- function(data){ 
  taumax <- tauloglik(data)$maxloglik$tau
  if(taumax==0.1){taumax <- taumax+0.01}
  tauloglik(data,taus=seq(taumax-0.1,taumax+0.1,0.01))$maxloglik$tau
}

#function to calculate the Sensitivity and specificity of case definitions using risk ratios
specSens <- function(data){
  #try values of tau between 0.1-0.9, then refine at differences of 0.01
  data$pfmcltau <- suppressWarnings(data$pfmcl^bestTau.fit(data)) 
  
  glmfit <- glm2(febrile ~ pfmcltau, family=binomial(link="logit"), data=data)
  #RR
  data$rr <- exp(glmfit$coeff[2] * data$pfmcltau)
  datafeb <- data[data$febrile==1,]      #subset of fever cases
  N <- nrow(datafeb)         #total no. of fever cases
  lambda <- 1/N*sum((datafeb$rr-1)/datafeb$rr) #proportion of fever cases attributable to malaria
  
  getAF <- function(c){   #c is parasite density cutoff    
    subdatafeb <- datafeb[datafeb$pfmcl>=c,] #subset of fever cases diagnosed as malaria cases by cutoff
    lambda_c <- 1/nrow(subdatafeb)*sum((subdatafeb$rr-1)/subdatafeb$rr)  #proportion of diagnosed cases 
    #that are attributable to malaria  
    N.lambda <- lambda*N             # total no of malaria attributable cases from the model
    n_c <- nrow(subdatafeb)  #No. fever cases diagnosed as malaria cases by cutoff
    n_c.lambda_c <- lambda_c*n_c           #no. of cases identified correctly by the case definition
    
    sens <- 100*(n_c.lambda_c/N.lambda)          #sensitivity
    
    num.spec <- n_c*(1-lambda_c)
    den.spec <- N*(1-lambda)
    spec <- 100*(1 - (num.spec/den.spec))              #specificity
    
    return(data.frame(threshold=c, AF=lambda_c*100, sensitivity=sens, specificity=spec))
  }
  return(data.frame(t(sapply(thresh,getAF))))
}

#function to bootstrap the data and generate confidence intervals for the MAFs at various thresholds
boot.AF <- function(data){
  data <- subset(data, select=c(febrile,pfmcl,agecat))
  
  boot.fn <- function(bootsample){
    bootsample$pfmcltau <- suppressWarnings(bootsample$pfmcl^bestTau.fit(bootsample))     #get tau
    glm.fit <- glm2(febrile ~ pfmcltau, family=binomial(link="logit"), data=bootsample)
    bootsample$rr <- exp(glm.fit$coeff[2] * bootsample$pfmcltau)      #odds ratio
    datafeb <- bootsample[bootsample$febrile==1,]       #subset with only cases
    
    lambda.c <- function(c){
      subdatafeb <- datafeb[datafeb$pfmcl>=c,] #subset of fever cases diagnosed as malaria cases by c
      lambda_c <- 1/nrow(subdatafeb)*sum((subdatafeb$rr-1)/subdatafeb$rr)  #proportion of cases diagnosed as malaria by c
    }
    
    AF <- sapply(thresh,lambda.c) 
    return(data.frame(thresh=thresh, AF=AF))
  }
  
  n <- nrow(data)
  
  #get AFs from the original data before running the bootstrap
  result <- try(ddply(data,.(agecat),boot.fn))
  R <- 100  
  #subresult <- matrix(nrow=length(thresh)*length(levels(data$agecat)) ,ncol=R)
  subresult <- matrix(nrow=length(thresh)*4 ,ncol=R)
  for(i in 1:R){
    bootsample <- data[sample(n, n, replace = TRUE),]    
    subresult[,i] <- ddply(bootsample,.(agecat),boot.fn)$AF
  }
  
  confint.mafs <- function(mafs){
    mafs.sorted <- sort(mafs)      #sort the MAFs estimates to obtain bootstrap CI
    #bootstrap percentile confidence intervals
    CI.mafs <- c(mafs.sorted[0.025*R], mafs.sorted[0.975*R+1])  
  }
  
  CI.mafs <- apply(subresult,1,confint.mafs) 
  
  return(cbind(result, lowerAF=CI.mafs[1,], upperAF=CI.mafs[2,]))
}