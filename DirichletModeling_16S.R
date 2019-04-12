#J. Harrison
set.seed(666)

library(rjags)

#read data
dat <- read.csv("./data/alldata2.csv", stringsAsFactors = F)

#Define start and end points for model indexing.
#These correspond to different treatment groups. 
#Note I am ignoring control treated plants for now. 
#Control plants were planted as seeds next to some agar with ALFU. 
#will model these at the bottom of the script.

dat <- dat[order(dat$treament),]
#remove NAs
dat <- dat[!is.na(dat[,84]),]

treatment <- dat$treament

orderModeled <- dat$plant
dat <- dat[,grep("^bactOtu\\d+$", names(dat))]
dim(dat)

dat <- dat[,-which(colSums(dat) == 0)]
dim(dat)
sum(colSums(dat))
#####################
#Model specification#
#####################

community.model.level <- "model{
  for(i in 1:N){
    for(j in start[i]:end[i]){
      datamatrix[j,] ~ dmulti(p[j,], nreads[j])
      p[j,1:notus] ~ ddirch(pi[i,]*theta[i])   
    }

    theta[i] ~ dunif(1.0E-3, max(nreads)) 
    pi[i,1:notus] ~ ddirch(alpha*hypertheta)
  }

  alpha ~ ddirch(hyperalpha)
  hypertheta ~ dunif(1.0E-3, max(nreads)) 

  for(k in 1:notus){
    hyperalpha[k] <-1/notus
  }
}"

####################
# Define functions #
####################

#This function computes start and end values for rjags
indexer <- function(x){
  start <- NA
  end <- NA
  k <- 1
  for(i in unique(x)){
    start[k] <- min(which(x == i))
    end[k] <- max(which(x == i ))
    k <- k +1
  }
  return(list(start, 
              end))
}
#Calculate MCMC diagnostics
mcmcdiag <- function(x) {
  #x is an mcmc object
  #nparams is number of params in the object
  Gr <- vector(length = dim(x[[1]])[2])
  GK <- vector(length = dim(x[[1]])[2])
  k <- 1
  a <- character(0)
  
  while (k <= dim(x[[1]])[2]) {
    m <- x[1:length(x)][, k]
    
    gr <- gelman.diag(m) #have to have at least two chains
    print(paste("Feature", k, sep = " "))
    print("Gelman-Rubin")
    print(gr)
    if (gr[[1]][1] <= 2) {
      Gr[k] <- "passed"
    } else{
      Gr[k] <- "failed"
    }
    
    gk <- geweke.diag(m,
                      frac1 = 0.1,
                      frac2 = 0.5)
    suspectGK <-  names(which(2 * pnorm(-abs(gk[[1]]$z)) < 0.08))
    if (identical(a, suspectGK)) {
      GK[k] <- "passed"
    } else if (suspectGK == "var1") {
      GK[k] <- "failed"
    }
    
    k <- k + 1
  }
  return(list(Gr, 
              GK))
}

modelRun <- function(dat, treatment){
  sim.mod.jags <- jags.model(textConnection(community.model.level),
                           data = list(
                             #Addition of a constant avoids infinite density slicer errors
                             datamatrix = 1 + dat,
                             notus = length(dat),               
                             nreads = rowSums(dat),  
                             N = 8,
                             start = indexer(treatment)[[1]],
                             end = indexer(treatment)[[2]]
                           ), 
                           #inits = initcalculator(dat, 
                           #                       indexer(treatment)[[1]],
                           #                       indexer(treatment)[[2]]),
                           n.adapt = 0,
                           n.chains = 2)

#Adapt model
iter_needed <- 0
y=FALSE
while(y==FALSE){
  y <-  adapt(sim.mod.jags,
              n.iter = 1000,
              end.adaptation=FALSE)
  iter_needed <- 1000 + iter_needed
  if(iter_needed > 30000){break}
}
iter_needed

update(sim.mod.jags,
       n.iter = 100000)

#Extract MCMC samples
sim.mod.sam <- jags.samples(model = sim.mod.jags,
                            variable.names = c(
                              "pi",
                              "alpha",
                              "p",
                              "theta"
                            ),
                            n.iter = 4000,
                            thin = 4)

#Determine MCMC performance
zed <- as.mcmc.list(sim.mod.sam$pi)
out <- mcmcdiag(x = zed)

#these are the ones that failed
fails <- dimnames(zed[[1]])[[2]][which(out[[1]] == "failed")]

return(list(sim.mod.sam, 
       out, 
       fails)
      )
}

####################
# Perform analysis #
####################

results<- modelRun(dat, treatment)
save.image("./BayesianRelativeAbundanceModelingBACTERIA.RData")