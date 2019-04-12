library(rjags)
library(randomForest)

load("./BayesianRelativeAbundanceModelingBACTERIA.RData")

dat2 <- read.csv("./data/ML_fungi_noALFU_withShannonSimpson.csv", stringsAsFactors = F)

#Bring in trait data
dat5 <- read.csv("./data/alldata2.csv", stringsAsFactors = F)

dat <- merge(dat5, 
             dat2,
             by.x = "plant",
             by.y = "samps")
dim(dat)

#Imputation of missing values for predictors of interest
#will be performed for each response variable
#but getting this in order, like so, helps streamline code later

dat4 <- data.frame(
  dat$treament,
  dat$plant,
  dat$shannon,
  dat$simpson,
  dat$height,
  dat$width,
  dat$length,
  dat$sla_avg,
  dat$percN,
  dat$percC,
  dat$ndfa,
  dat$X..Swainsonine
)
names(dat4) <-  c(
  "treatment",
  "samps",
  "shannon",
  "simpsons",
  "height",
  "width",
  "length",
  "sla_avg",
  "percN",
  "percC",
  "ndfa",
  "swain")

dat4 <- dat4[order(dat4$treatment),]  
####################
# Define Model Strings #
####################

linearModel <- "model{for(i in 1:N){
 for(j in start[i]:end[i]){ 
  response[j]~dnorm(alpha[j], tau.y[i])
  
  alpha[j] <- mu[i] + 
    beta1[i]*term1[j]+
    beta2[i]*term2[j]+
    beta3[i]*term3[j]+
    beta4[i]*term4[j]+
    beta5[i]*term5[j]
    #beta6[i]*term6[j]
    #beta7[i]*term7[j]+
    #beta8[i]*term8[j]
  }    
  
  tau.y[i] ~ dgamma(.01,.01)
  
  #priors for beta coefficients
  beta1[i] ~ dnorm(beta1mu, beta1tau)
  beta2[i] ~ dnorm(beta2mu, beta2tau)
  beta3[i] ~ dnorm(beta3mu, beta3tau)
  beta4[i] ~ dnorm(beta4mu, beta4tau)
  beta5[i] ~ dnorm(beta5mu, beta5tau)
  #beta6[i] ~ dnorm(beta6mu, beta6tau)
  #beta7[i] ~ dnorm(beta7mu, beta7tau)
  #beta8[i] ~ dnorm(beta8mu, beta8tau)
  mu[i] ~ dnorm(mumu, mutau)
  }
  #hyperpriors
  beta1tau ~ dgamma(0.1,0.1) 
  beta2tau ~ dgamma(0.1,0.1)
  beta3tau ~ dgamma(0.1,0.1)
  beta4tau ~ dgamma(0.1,0.1)
  beta5tau ~ dgamma(0.1,0.1)
 # beta6tau ~ dgamma(0.1,0.1)
  #beta7tau ~ dgamma(0.1,0.1)
  #beta8tau ~ dgamma(0.1,0.1)
  
  mutau ~ dgamma(0.1,0.1)
  
  beta1mu ~ dnorm(0,0.001)	
  beta2mu ~ dnorm(0,0.001)
  beta3mu ~ dnorm(0,0.001)
  beta4mu ~ dnorm(0,0.001)
  beta5mu ~ dnorm(0,0.001)
  #beta6mu ~ dnorm(0,0.001)
  #beta7mu ~ dnorm(0,0.001)
  #beta8mu ~ dnorm(0,0.001)
  
  mumu ~ dnorm(0,0.001)
}"

####################
# Define functions #
####################

results <- function(x, m){
  out <- matrix(ncol = 4,
                nrow = 9)
  
  
  #calculate avg and > 0 stat for all samples across all treatments
  out[1, 1] <- round(mean(x[, , 1:2]),2)
  out[1, 2] <-
    round(100 * length(which(x[, , 1:2] > 0)) / length(x[, , 1:2]), 
          digits = 1)
  out[1, 3] <- paste(paste(out[1, 1], "(", sep = " "),
                     out[1, 2], "%)",
                     sep = "")
  out[1,4] <- "acrossTreatments"
  
  k <- 2
  for(i in 1:m){
    # use this for realtime feedback
    # print(i)
    # print("ML estimate of beta")
    # print(mean(x[i,,1:2]))
    # print("perc > 0")
    # print(length(which(x[i,,1:2]>0))/length(x[i,,1:2]))
    # print(quantile(x[i,,1:2], c(0.025,0.975)))
    out[k,1] <- round(mean(x[i,,1:2]),2)
    out[k,2] <- round(100*length(which(x[i,,1:2] > 0))/length(x[i,,1:2]),digits = 1)
    out[k,3] <- paste(paste(out[k,1],"(", sep = " "),
                      out[k,2], "%)",
                      sep = "")
    out[k,4] <- unique(treatment)[i]
    k <- k + 1
  }
  return(out)
}

#define start and end points for model
bounder <- function(x) {
  #x is a vector of treatments, it should be ordered!
  start <- vector()
  end <- vector()
  for (i in unique(x)) {
    start <- c(start, min(which(x == i)))
    end <- c(end, max(which(x == i)))
  }
  return(list(start,
              end))
}

modelRun <- function (x, datz, starts, ends){
  sim.mod.jags <- jags.model(
  textConnection(linearModel),
  data = list(
    response = x,
    term1 = c(scale(
      datz$height * datz$width * datz$length,
      center = T,
      scale = T
      )),
    term2 = c(scale(datz$sla_avg,
                  center = T,
                  scale = T)),
    term3 = c(scale(datz$percN,
                  center = T,
                  scale = T)),
    term4 = c(scale(datz$percC,
                  center = T,
                  scale = T)),
    term5 = c(scale(datz$ndfa,
                  center = T,
                  scale = T)),
    start = unlist(starts),
    end = unlist(ends),
    N = 8
  ),
  n.chains = 2,
  n.adapt = 0
)

#adapt
iter_needed <- 0
y=FALSE
while(y==FALSE){
  y <-  adapt(sim.mod.jags, 
              n.iter=1000,
              end.adaptation=FALSE)
  iter_needed <- 1000 + iter_needed
  print(iter_needed)
  if(iter_needed > 8000){break}
}

#burn
update(sim.mod.jags, 
       n.iter = 5000)

#sample
#if one were interested in modeling variance here, then one could do that. 
#I dont see any reason to compare between and within variance as one does in an Anova
#to determine if means differ because we are directly modeling the means themselves. 

sim.mod.sam <- jags.samples(model=sim.mod.jags, 
                            variable.names=c(
                              "mu",
                              "beta1",
                              "beta2",
                              "beta3",
                              "beta4",
                              "beta5",
                              #"beta6",
                              #"beta7",
                              #"beta8",
                              "beta1mu",
                              "beta2mu",
                              "beta3mu",
                              "beta4mu",
                              "beta5mu"
#                              "beta6mu"
                              #"beta7mu",
                              #"beta8mu"
                            ), 
                            n.iter=4000,
                            thin=4)

#Diagnostic statistics. 

#Calculate the Gelman-Rubin statistic, spot check
gr <- gelman.diag(as.mcmc.list(sim.mod.sam$beta4), 
                  multivariate=F)
print(gr)

gk <- geweke.diag(as.mcmc.list(sim.mod.sam$beta4), 
                  frac1=0.1, frac2=0.5) 
print(gk)
#see which parameters didn't have a long enough burn in. 
print(names(which(2*pnorm(-abs(gk[[1]]$z))<0.05)))

print( results(sim.mod.sam$beta1,8))
print(  results(sim.mod.sam$beta2,8))
print(  results(sim.mod.sam$beta3,8))
print(  results(sim.mod.sam$beta4,8))
print(  results(sim.mod.sam$beta5,8))
print(  results(sim.mod.sam$mu,8))

}

###############
# Do analysis #
###############
#Use this to avoid imputation
# dat3 <- na.omit(dat4)
# rmv <- unique(which(is.na(dat4), arr.ind = T)[,1])
# 
# starts <-  bounder(dat$treament[-rmv])[1]
# ends <- bounder(dat$treament[-rmv])[2]
#modelRun(dat2[-rmv,26], unlist(starts), unlist(ends))

#Impute data that is missing. 
dat5 <- dat4
test <- rfImpute(x = dat5[,5:11], y = dat5$shannon)
dat5[,5:11] <- test[,2:8]

starts <-  bounder(dat5$treatment)[1]
ends <- bounder(dat5$treatment)[2]

modelRun(x = dat5$shannon,
         datz = dat5,
         starts = unlist(starts), 
         ends = unlist(ends))

dat5 <- dat4
test <- rfImpute(x = dat5[,5:11], y = dat5$simpsons)
dat5[,5:11] <- test[,2:8]

starts <-  bounder(dat5$treatment)[1]
ends <- bounder(dat5$treatment)[2]

modelRun(x = dat5$simpsons,
         datz = dat5,
         starts = unlist(starts), 
         ends = unlist(ends))

############################################################
#Repeat with swainsonine, but for only Afulva + treatments#
############################################################

linearModel2 <- "model{for(i in 1:N){
 for(j in start[i]:end[i]){ 
  response[j]~dnorm(alpha[j], tau.y[i])
  
  alpha[j] <- mu[i] + 
    beta1[i]*term1[j]+
    beta2[i]*term2[j]+
    beta3[i]*term3[j]+
    beta4[i]*term4[j]+
    beta5[i]*term5[j]+
    beta6[i]*term6[j]
    #beta7[i]*term7[j]+
    #beta8[i]*term8[j]
  }    
  
  tau.y[i] ~ dgamma(.01,.01)
  
  #priors for beta coefficients
  beta1[i] ~ dnorm(beta1mu, beta1tau)
  beta2[i] ~ dnorm(beta2mu, beta2tau)
  beta3[i] ~ dnorm(beta3mu, beta3tau)
  beta4[i] ~ dnorm(beta4mu, beta4tau)
  beta5[i] ~ dnorm(beta5mu, beta5tau)
  beta6[i] ~ dnorm(beta6mu, beta6tau)
  #beta7[i] ~ dnorm(beta7mu, beta7tau)
  #beta8[i] ~ dnorm(beta8mu, beta8tau)
  mu[i] ~ dnorm(mumu, mutau)
  }
  #hyperpriors
  beta1tau ~ dgamma(0.1,0.1) 
  beta2tau ~ dgamma(0.1,0.1)
  beta3tau ~ dgamma(0.1,0.1)
  beta4tau ~ dgamma(0.1,0.1)
  beta5tau ~ dgamma(0.1,0.1)
  beta6tau ~ dgamma(0.1,0.1)
  #beta7tau ~ dgamma(0.1,0.1)
  #beta8tau ~ dgamma(0.1,0.1)
  
  mutau ~ dgamma(0.1,0.1)
  
  beta1mu ~ dnorm(0,0.001)	
  beta2mu ~ dnorm(0,0.001)
  beta3mu ~ dnorm(0,0.001)
  beta4mu ~ dnorm(0,0.001)
  beta5mu ~ dnorm(0,0.001)
  beta6mu ~ dnorm(0,0.001)
  #beta7mu ~ dnorm(0,0.001)
  #beta8mu ~ dnorm(0,0.001)
  
  mumu ~ dnorm(0,0.001)
}"

modelRun_swain <- function(x, datz, starts, ends){
  sim.mod.jags <- jags.model(
  textConnection(linearModel2),
  data = list(
    response = x,
    term1 = c(scale(
      datz$height * datz$width * datz$length,
      center = T,
      scale = T
    )),
    term2 = c(scale(
      datz$sla_avg,
      center = T,
      scale = T
    )),
    
    term3 = c(scale(
      datz$percN,
      center = T,
      scale = T
    )),
    term4 = c(scale(
      datz$percC,
      center = T,
      scale = T
    )),
    term5 = c(scale(
      datz$ndfa,
      center = T,
      scale = T
    )),
     term6 = c(scale(
      as.numeric(datz$swain),
      center = T,
      scale = T
    )),
    N = 4,
  start = starts,
  end = ends),
  n.chains = 2,
  n.adapt = 0
)

#adapt
iter_needed <- 0
y=FALSE
while(y==FALSE){
  y <-  adapt(sim.mod.jags, 
              n.iter=1000,
              end.adaptation=FALSE)
  iter_needed <- 1000 + iter_needed
  print(iter_needed)
  if(iter_needed > 3000){break}
}

#burn

update(sim.mod.jags, 
       n.iter=5000)

#sample
#if one were interested in modeling variance here, then one could do that. 
#I dont see any reason to compare between and within variance as one does in an Anova
#to determine if means differ because we are directly modeling the means themselves. 

sim.mod.sam <- jags.samples(model=sim.mod.jags, 
                            variable.names=c(
                              "mu",
                              "beta1",
                              "beta2",
                              "beta3",
                              "beta4",
                              "beta5",
                              "beta6",
                              #"beta7",
                              #"beta8",
                              "beta1mu",
                              "beta2mu",
                              "beta3mu",
                              "beta4mu",
                              "beta5mu",
                             "beta6mu"
                              #"beta7mu",
                              #"beta8mu"
                            ), 
                            n.iter=4000,
                            thin=4)

#Diagnostic statistics. 

#Calculate the Gelman-Rubin statistic, spot check
gr <- gelman.diag(as.mcmc.list(sim.mod.sam$beta4), 
                  multivariate=F)
print(gr)

gk <- geweke.diag(as.mcmc.list(sim.mod.sam$beta4), 
                  frac1=0.1, frac2=0.5) 
gk
#see which parameters didn't have a long enough burn in. 
names(which(2*pnorm(-abs(gk[[1]]$z))<0.05))

print(results(sim.mod.sam$beta1,4))
print(results(sim.mod.sam$beta2,4))
print(results(sim.mod.sam$beta3,4))
print(results(sim.mod.sam$beta4,4))
print(results(sim.mod.sam$beta5,4))
print(results(sim.mod.sam$beta6,4))
print(results(sim.mod.sam$mu,4))
}

#wrangle, impute, etc.
#Impute data that is missing. 
#Do not impute swainsonine, bc I 
#Do not think imputation will be very accurate.

rmv <- which(is.na(dat4$swain))
#TO avoid imputation:
# rmv <- unique(which(is.na(dat4), arr.ind = T)[,1])
# dat3 <- dat4[-rmv,]

dat4 <- dat4[-rmv,]

dat5 <- dat4
test <- rfImpute(x = dat5[,5:11], 
                 y = dat5$shannon)
dat5[,5:11] <- test[,2:8]

starts <-  bounder(dat5$treatment)[1]
ends <- bounder(dat5$treatment)[2]

modelRun_swain(x = dat5$shannon,
         datz = dat5,
         starts = unlist(starts), 
         ends = unlist(ends))

#REMEMBER the first row in the output is the mean
dat5 <- dat4
test <- rfImpute(x = dat5[,5:11], 
                 y = dat5$simpsons)
dat5[,5:11] <- test[,2:8]

starts <-  bounder(dat5$treatment)[1]
ends <- bounder(dat5$treatment)[2]

modelRun_swain(x = dat5$simpsons,
         datz = dat5,
         starts = unlist(starts), 
         ends = unlist(ends))

