
newdat <- read.csv("data/alldata2.csv")

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

pdf(width = 6, height = 8, file = "./visuals/culturingPlot.pdf")
par(oma = c(9,4,2,0),
    mar = c(1,5,1,0),
    mfrow = c(3,1))

#plot morphospecies
stripchart(newdat$morphospecies_not_counting_Undifilum ~ newdat$treament,
           vertical = TRUE,
           data = dat,
           method = "jitter",
           pch = 20,
           bty = "n",
           cex.lab = 1.5,
           col = add.alpha("gray", alpha = 0.8),
           cex=2,
           xaxt="n",
           yaxt="n",
           ylab="Morphospecies",
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)

boxplot(newdat$morphospecies_not_counting_Undifilum ~ newdat$treament, 
        las = 2, 
        xaxt="n",
        frame.plot=F,
        col= c(NA, NA, 
               add.alpha("light gray", alpha=0.6),
               add.alpha("light gray", alpha=0.6),
               add.alpha("light gray", alpha=0.6),
               add.alpha("light gray", alpha=0.6),
               NA,NA),
        add = T,
        outline = F)
mtext("a)", side = 2,
      line = 2.5, 
      padj = -4,
      cex = 2,
      las = 2,
      xpd = NA)

#Plot percentage colonization
stripchart(newdat$notAlfuInfected/newdat$num_leaf_segments ~ newdat$treament,
           vertical = TRUE,
           data = dat,
           method = "jitter",
           pch = 20,
           col = add.alpha("gray", alpha = 0.8),
           cex=2,
           xaxt="n",
           cex.lab = 1.5,
           yaxt="n",
           ylab="% infected",
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)

boxplot(newdat$notAlfuInfected/newdat$num_leaf_segments ~ newdat$treament, 
        las = 2, 
        xaxt="n",
        col= c(NA, NA, 
               add.alpha("light gray", alpha=0.6),
               add.alpha("light gray", alpha=0.6),
               add.alpha("light gray", alpha=0.6),
               add.alpha("light gray", alpha=0.6),
               NA,NA),
        frame.plot=F,
        add = T,
        outline = F)
mtext("b)", side = 2,
      line = 2.5, 
      padj = -4,

      cex = 2,
      las = 2)

#Plot A. fulva colonization
stripchart(newdat$Undifilum/newdat$num_leaf_segments ~ newdat$treament,
           vertical = TRUE,
           data = dat,
           method = "jitter",
           pch = 20,
           col = add.alpha("gray", alpha = 0.8),
           cex=2,
           xaxt="n",
           cex.lab = 1.5,
           yaxt="n",
           bty = "n",
           ylab= substitute(paste("% ", italic("A. fulva"), " infected", sep = "")),
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)

boxplot(newdat$Undifilum/newdat$num_leaf_segments ~ newdat$treament, 
        las = 2, 
        col= c(NA, NA, 
               add.alpha("light gray", alpha=0.6),
               add.alpha("light gray", alpha=0.6),
               add.alpha("light gray", alpha=0.6),
               add.alpha("light gray", alpha=0.6),
               NA,NA),
        add = T,
        frame.plot=F,
        outline = F, 
        names = c("Control -", 
                   "Control +", 
                   "Treated, Control - ",
         "Treated, Control + ",
        "Treated - ",
         "Treated +",
        "Untreated -",
        "Untreated +"
        ),
        cex = 1.5)

mtext("c)", side = 2,
      line = 2.5, 
      padj = -4,
      cex = 2,
      las = 2)

dev.off()


##################
#Perform analysis#
##################


meanModeler <- "model{
#model likelihood by sampling group
for(j in 1:groups){
for(i in start[j]:end[j]){
response[i] ~ dnorm(mu[j], tau[j])
}
#priors
mu[j] ~ dnorm(0, 0.0001)
tau[j] <- pow(sigma[j], -2)
sigma[j] ~ dunif(0, 100)
}
#calculate differences between sampling groups
delta12 <- mu[1]-mu[2]
delta13 <- mu[1]-mu[3]
delta14 <- mu[1]-mu[4]
delta23 <- mu[2]-mu[3]
delta24 <- mu[2]-mu[4]
delta34 <- mu[3]-mu[4]
delta56 <- mu[5]-mu[6]
delta78 <- mu[7]-mu[8]
}"

modeler <- function(x){
  sim.mod.jags <- jags.model(textConnection(meanModeler),
                             data = list(
                               response = x,
                               groups = 8,
                               start = c(1,9,14,17,21,57,88,119),
                               end = c(8,13,16,20,56,87,118,153)
                             ), 
                             n.chains=2, 
                             n.adapt=0)
  
  
  #adapt
  iter_needed <- 0
  y=FALSE
  while(y==FALSE){
    y <-  adapt(sim.mod.jags, 
                n.iter=1000,
                end.adaptation=FALSE)
    iter_needed <- 1000 + iter_needed
    if(iter_needed > 5000){break}
  }
  
  #burn
  burn<-5000
  update(sim.mod.jags, 
         n.iter=burn)
  
  #sample
  #if one were interested in modeling variance here, then one could do that. 
  #I dont see any reason to compare between and within variance as one does in an Anova
  #to determine if means differ because we are directly modeling the means themselves. 
  
  sim.mod.sam <- jags.samples(model=sim.mod.jags, 
                              variable.names=c(
                                "mu",
                                "delta12",
                                "delta13",
                                "delta14",
                                "delta23",
                                "delta24",
                                "delta34",
                                "delta56",
                                "delta78"
                              ), 
                              n.iter=4000,
                              thin=4)
  
  #Diagnostic statistics. 
  #I do not calculate these for delta, because delta depends on mu. If mus are
  #good, then deltas are good. 
  
  #Calculate the Gelman-Rubin statistic
  
  gr <- gelman.diag(as.mcmc.list(sim.mod.sam$mu), 
                    multivariate=F)
  print(gr)
  
  
  #Calculate temporal autocorrelation in the chains
  #this function determines if the mean between the first and last parts of the chain differ
  #if they don't differ then burn in was long enough. 
  #The values output are z scores (diff in means divided by SE)
  #The "frac" options  determine the fractions of the chains to compare
  
  gk <- geweke.diag(as.mcmc.list(sim.mod.sam$mu), 
                    frac1=0.1, frac2=0.5) 
  print(gk)
  
  print("Parameters that didn't have a long enough burn in.")
  print(names(which(2*pnorm(-abs(gk[[1]]$z))<0.05)))
  return(sim.mod.sam)
}

###################
#Get them Results!#
###################

newdat <- newdat[order(newdat$treament),]
sim.mod.sam <- modeler(newdat$morphospecies_not_counting_Undifilum)
levels(newdat$treament)

# quantile(sim.mod.sam$delta12, probs = c(0.05,0.5,0.95))
# quantile(sim.mod.sam$delta34, probs = c(0.05,0.5,0.95))
quantile(sim.mod.sam$delta78, probs = c(0.05,0.5,0.95))
quantile(sim.mod.sam$delta56, probs = c(0.05,0.5,0.95))

sim.mod.sam <- modeler(newdat$Undifilum/newdat$num_leaf_segments)
quantile(sim.mod.sam$delta78, probs = c(0.05,0.5,0.95))
quantile(sim.mod.sam$delta56, probs = c(0.05,0.5,0.95))

sim.mod.sam <- modeler(newdat$notAlfuInfected/newdat$num_leaf_segments)
quantile(sim.mod.sam$delta78, probs = c(0.05,0.5,0.95))
quantile(sim.mod.sam$delta56, probs = c(0.05,0.5,0.95))