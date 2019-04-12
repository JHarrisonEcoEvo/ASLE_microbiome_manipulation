#Here we estimate the population means for each treatment group for each trait of interest
#note that we model variance separately for each treatment group.
#results are plotted at the bottom of this script. 

library(rjags)
set.seed(666)

#read data
dat <- read.csv("./data/alldata2.csv", stringsAsFactors = F)

#Define start and end points for model indexing.
#These correspond to different treatment groups. 
#Note I am ignoring control treated plants for now. 
#Control plants were planted as seeds next to some agar with ALFU. 
#will model these at the bottom of the script.

dat <- dat[order(dat$treament),]
dat <- dat[which(dat$treament %in% c("treated_neg",
                                     "untreated_neg", 
                                     "treated_plus", 
                                     "untreated_plus")),
           ]

#dat$treament <- droplevels(dat$treament)

#Use this to remove swainsonine positive plants with reduced a. fulva from consideration
#these are likely those plants for which treatment was unsuccessful.
#dat <- dat[!(dat$X..Swainsonine > 0.001 & dat$treament %in% c("untreated_neg","treated_neg")),]

#Use this to keep those plants but recode them
# tochange <- which(dat$X..Swainsonine > 0.001 & dat$treament == "treated_neg")
# dat$treament[tochange] <- "treated_plus"
# tochange <- which(dat$X..Swainsonine > 0.001 & dat$treament == "untreated_neg")
# dat$treament[tochange] <- "untreated_plus"
# 
# #use this to remove swainsonine negative plants that were supposed to have A. fulva
# dat <- dat[-which(dat$X..Swainsonine < 0.001 & dat$treament == "treated_plus"),]
# dat <- dat[-which(dat$X..Swainsonine < 0.001 & dat$treament == "untreated_plus"),]
# dat <- dat[order(dat$treament),]

indexer <- function(x){
  start <- c(min(which(x =="treated_neg")),
             min(which(x =="treated_plus")),
             min(which(x=="untreated_neg")),
             min(which(x =="untreated_plus"))
  )
  end <- c(max(which(x =="treated_neg")),
           max(which(x =="treated_plus")),
           max(which(x =="untreated_neg")),
           max(which(x =="untreated_plus"))
  ) 
  return(list(start, end))
}

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
}"

#run model function
modelRun <- function(x, start, end){
  sim.mod.jags <- jags.model(textConnection(meanModeler),
                             data = list(
                               response = x,
                               groups = 4,
                               start = start,
                               end = end
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
                                "sigma",
                                "tau"
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
  
  ###################
  #Get them Results!#
  ###################
  
  #Get estimates of mus and deltas
  
  calcCertainty <- function(x){
    plot(density(x))
    abline(v = 0)
    
    #calculate estimate of certainty for the threshold lovers
    p <- 1-length(which(x > 0))/length(x)
    return(p)
  }
  print("Comparisons")
  print(calcCertainty(sim.mod.sam$delta12))
  print(calcCertainty(sim.mod.sam$delta13))
  print(calcCertainty(sim.mod.sam$delta14))
  print(calcCertainty(sim.mod.sam$delta23))
  print(calcCertainty(sim.mod.sam$delta24))
  print(calcCertainty(sim.mod.sam$delta34))
  
  return(list(quantiles = list(rowMeans(apply(sim.mod.sam$mu[1,,1:2], 2, quantile, 
                                              probs = c(0.025,0.5, 0.975))),
                               #Treatment group 2
                               rowMeans(apply(sim.mod.sam$mu[2,,1:2], 2, quantile, 
                                              probs = c(0.025,0.5, 0.975))),
                               
                               #Treatment group 3
                               rowMeans(apply(sim.mod.sam$mu[3,,1:2], 2, quantile, 
                                              probs = c(0.025,0.5, 0.975))),
                               
                               #Treatment group 4
                               rowMeans(apply(sim.mod.sam$mu[4,,1:2], 2, quantile, 
                                              probs = c(0.025,0.5, 0.975)))),
              means = list(mean(apply(sim.mod.sam$mu[1,,1:2], 2, FUN=mean)),
                           mean(apply(sim.mod.sam$mu[2,,1:2], 2, FUN=mean)),
                           mean(apply(sim.mod.sam$mu[3,,1:2], 2, FUN=mean)),
                           mean(apply(sim.mod.sam$mu[4,,1:2], 2, FUN=mean))),
              tau = sim.mod.sam$tau
  ))
}

out_mat <- matrix(nrow = 9, ncol = 5)

levels(dat$treament)
dat$vol <- dat$width*dat$height*dat$Length

dat_vol <- dat[-which(is.na(dat$width)),]
start <- indexer(dat_vol$treament)[[1]]
end <- indexer(dat_vol$treament)[[2]]
out <- modelRun(dat_vol$vol, start, end)

out_mat[1,] <- c("Size", 
                 paste(round(unlist(out$means),2),
                       " (",
                       round(sapply(out$quantiles, "[[", 1), 2),
                       ",",
                       round(sapply(out$quantiles, "[[", 3), 2),
                       ")", 
                       sep=""))

dat_leaves <- dat[-which(is.na(dat$X.LeavesCollection)),]
start <- indexer(dat_leaves$treament)[[1]]
end <- indexer(dat_leaves$treament)[[2]]
out <- modelRun(dat$X.LeavesCollection, start, end)
out_mat[2,] <- c("Leaves", 
                 paste(round(unlist(out$means),2),
                       " (",
                       round(sapply(out$quantiles, "[[", 1), 2),
                       ",",
                       round(sapply(out$quantiles, "[[", 3), 2),
                       ")", 
                       sep=""))

dat_leafarea <- dat[-which(is.na(dat$leafarea_avg)),]
start <- indexer(dat_leafarea$treament)[[1]]
end <- indexer(dat_leafarea$treament)[[2]]
out <- modelRun(dat$leafarea_avg, start, end)
out_mat[3,] <- c("Leaf area",  
                 paste(round(unlist(out$means),2),
                       " (",
                       round(sapply(out$quantiles, "[[", 1), 2),
                       ",",
                       round(sapply(out$quantiles, "[[", 3), 2),
                       ")", 
                       sep=""))

dat_sla <- dat[-which(is.na(dat$sla_avg)),]
start <- indexer(dat_sla$treament)[[1]]
end <- indexer(dat_sla$treament)[[2]]
out <- modelRun(dat$sla_avg, start, end)
out_mat[4,] <- c("SLA", 
                 paste(round(unlist(out$means),2),
                       " (",
                       round(sapply(out$quantiles, "[[", 1), 2),
                       ",",
                       round(sapply(out$quantiles, "[[", 3), 2),
                       ")", 
                       sep=""))

dat_15 <- dat[-which(is.na(dat$d15N....vs..air.)),]
start <- indexer(dat_15$treament)[[1]]
end <- indexer(dat_15$treament)[[2]]
out <- modelRun(dat$d15N....vs..air., start, end)
out_mat[5,] <- c("delta N15", 
                 paste(round(unlist(out$means),2),
                       " (",
                       round(sapply(out$quantiles, "[[", 1), 2),
                       ",",
                       round(sapply(out$quantiles, "[[", 3), 2),
                       ")", 
                       sep=""))

#Use this to remove the a. fulva plus plants with no swainsonine
#It is less likely that these represent failures in treatment, but it is plausible
dat_swain <- dat[!is.na(dat$X..Swainsonine),]
dat_swain$X..Swainsonine[dat_swain$X..Swainsonine == "#VALUE!"] <- 0
dat_swain$X..Swainsonine <- as.numeric(dat_swain$X..Swainsonine)

# dat_swain <- dat_swain[-which(dat_swain$X..Swainsonine < 0.001 & dat_swain$treament == "treated_plus"),]
# dat_swain <- dat_swain[-which(dat_swain$X..Swainsonine < 0.001 & dat_swain$treament == "untreated_plus"),]
# dat_swain <- dat_swain[order(dat_swain$treament),]

boxplot(as.numeric(dat_swain$X..Swainsonine) ~ dat_swain$treament)
summary(lm(as.numeric(dat_swain$X..Swainsonine) ~ dat_swain$treament))
start <- indexer(dat_swain$treament)[[1]]
end <- indexer(dat_swain$treament)[[2]]
out <- modelRun(as.numeric(dat_swain$X..Swainsonine), start, end)
out_mat[6,] <- c("% swainsonine",
                 paste(round(unlist(out$means),2),
                       " (",
                       round(sapply(out$quantiles, "[[", 1), 2),
                       ",",
                       round(sapply(out$quantiles, "[[", 3), 2),
                       ")", 
                       sep=""))

dat_percN <- dat[-which(is.na(dat$percN)),]
start <- indexer(dat_percN$treament)[[1]]
end <- indexer(dat_percN$treament)[[2]]
out <- modelRun(dat_percN$percN, start, end)
out_mat[7,] <- c("% N", 
                 paste(round(unlist(out$means),2),
                       " (",
                       round(sapply(out$quantiles, "[[", 1), 2),
                       ",",
                       round(sapply(out$quantiles, "[[", 3), 2),
                       ")", 
                       sep=""))

dat_percC <- dat[-which(is.na(dat$percC)),]
start <- indexer(dat_percC$treament)[[1]]
end <- indexer(dat_percC$treament)[[2]]
out <- modelRun(dat$percC, start, end)
out_mat[8,] <- c("% C",  
                 paste(round(unlist(out$means),2),
                       " (",
                       round(sapply(out$quantiles, "[[", 1), 2),
                       ",",
                       round(sapply(out$quantiles, "[[", 3), 2),
                       ")", 
                       sep=""))
dat_ndfa <- dat[-which(is.na(dat$ndfa)),]
start <- indexer(dat_ndfa$treament)[[1]]
end <- indexer(dat_ndfa$treament)[[2]]
out <- modelRun(dat$ndfa, start, end)
out_mat[9,] <- c("NDFA",
                 paste(round(unlist(out$means),2),
                       " (",
                       round(sapply(out$quantiles, "[[", 1), 2),
                       ",",
                       round(sapply(out$quantiles, "[[", 3), 2),
                       ")", 
                       sep=""))
library(xtable)
print(xtable(out_mat,
             type = "latex"),
      file = "./traitMeans.tex",
      floating.environment='table',
      include.rownames = FALSE)

################################################################
# Time to draw pretty pictures!
################################################################

dat <- read.csv("./data/alldata2.csv", header = T)

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

colz <- c("cyan4", "cyan4", "darkgoldenrod2", "darkgoldenrod2")

dat <- dat[-grep("control",dat$treament),]
dat$treament <- droplevels(dat$treament)
# 
#plot
pdf(width=8,height=6, "./visuals/Fig1treatmentEffectsBoxplots.pdf")

dat$vol <- dat$height*dat$width*dat$length

par(oma=c(2,2,1,0), mar=c(3,4.5,1,1), mfrow=c(3,3))
stripchart(dat$vol ~ droplevels(dat$treament),
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
           ,ylim=c(0,75000)
           ,frame.plot=F
)
boxplot(dat$vol ~ droplevels(dat$treament), 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
text(c("ac","ab","b","c"),
     x = c(1,2,3,4),
     y = c(75000),
     cex = 1.2,
     font = 2,
     xpd = NA)

axis(1,las=2, labels = c("","","",""), at = c(1,2,3,4))
axis(2,las=2, 
     at = c(seq(0,75000,by = 25000)),
     labels = c(seq(0,75000,by = 25000))
)
title(ylab=expression(paste("Plant volume (",cm^3,")",sep="")), cex.lab=1.5, line=3.2,xpd = NA)

stripchart(dat$X.LeavesCollection ~ dat$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
           ,ylim=c(0,180)
           ,frame.plot=F
)
boxplot(dat$X.LeavesCollection ~ dat$treament, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
text(c("a","a","b","a"),
     x = c(1,2,3,4),
     y = 180,
     cex = 1.2,
     font = 2,
     xpd = NA)

axis(1,las=2, labels = c("","","",""), at = c(1,2,3,4))
axis(2,las=2,
     at = c(seq(0,180,by = 60)),
     labels = c(seq(0,180,by = 60))
)
title(ylab="Num. of leaves", cex.lab=1.5, line=2.8)
#with outliers: untreated neg different from untreated plus, treated neg, and almost different from treated plus
#wout outliers gets more complex. 


dat$X..Swainsonine <- as.character(dat$X..Swainsonine)
dat$X..Swainsonine[which(dat$X..Swainsonine=="#VALUE!")] <- 0

stripchart(as.numeric(dat$X..Swainsonine) ~ droplevels(dat$treament),
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
           ,ylim=c(0,0.12)
           ,frame.plot=F
)
boxplot(as.numeric(dat$X..Swainsonine) ~ droplevels(dat$treament), 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
axis(1,las=2, labels = c("","","",""), at = c(1,2,3,4))
axis(2,las=2,
     at = c(seq(0,.12,by = 0.02)),
     labels = c(seq(0,.12,by = 0.02))
)
title(ylab=" % Swainsonine", cex.lab=1.5, line=2.8)

text(c("a","b","a","b"),
     x = c(1,2,3,4),
     y = 0.13,
     cex = 1.2,
     font = 2,
     xpd = NA)

stripchart(dat$leafarea_avg ~ droplevels(dat$treament),
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
           ,ylim=c(0.3,0.8)
           ,frame.plot=F
)
boxplot(dat$leafarea_avg ~ droplevels(dat$treament), 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)

text(c("a","b","b","b"),
     x = c(1,2,3,4),
     y = 0.8,
     cex = 1.2,
     font = 2,
     xpd = NA)

axis(1,las=2, labels = c("","","",""), at = c(1,2,3,4))
axis(2,las=2, 
     at = c(seq(0.3,0.8,by = 0.1)),
     labels = c(seq(0.3,0.8,by = 0.1)))
title(ylab=expression(paste("Leaf area (c",m^2,")",sep="")), cex.lab=1.5, line=2.8)

stripchart(dat$sla_avg ~ droplevels(dat$treament),
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
           ,ylim=c(0.08,0.18)
           ,frame.plot=F
)
boxplot(dat$sla_avg ~ droplevels(dat$treament), 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
text("No certain effect\nof treatment",
     x = 2.5,
     y = 0.15,
     cex = 1.2,
     font = 2,
     xpd = NA)
axis(1,las=2, labels = c("","","",""), at = c(1,2,3,4))
axis(2,las=2,
     at = c(seq(0.08,0.18,by = 0.02)),
     labels = c(seq(0.08,0.18,by = 0.02)))
title(ylab=expression(paste("(c",m^2,"/mg)",sep="")), 
      cex.lab=1.3, line=2.8)
mtext("Specific leaf area", 
      line = 4.2,
      side = 2)


stripchart(dat$percC ~ droplevels(dat$treament),
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
           ,ylim=c(9,16)
           ,frame.plot=F
)
boxplot(dat$percC ~ droplevels(dat$treament), 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
text("No certain effect\nof treatment",
     x = 2.5,
     y = 15,
     cex = 1.2,
     font = 1,
     xpd = NA)

axis(1,las=2, labels = c("","","",""), at = c(1,2,3,4))
axis(2,las=2,
     at = c(seq(9,16,by = 1)),
     labels = c(seq(9,16,by = 1)))
title(ylab="% C", cex.lab=1.5, line=2.8)

stripchart(dat$percN ~ droplevels(dat$treament),
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
           ,ylim=c(0.4,1.8)
           ,frame.plot=F
)
boxplot(dat$percN ~ droplevels(dat$treament), 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)

text(c("a","ab","b","a"),
     x = c(1,2,3,4),
     y = 1.8,
     cex = 1.2,
     font = 2,
     xpd = NA)
axis(1,las=2, labels = c("","","",""), at = c(1,2,3,4))
axis(2,las=2, at = c(seq(0.6,1.8, by = 0.2)),
     labels = c(seq(0.6,1.8, by = 0.2)))
title(ylab="% N", cex.lab=1.5, line=2.8)
#with outliers untreated plus and untreated neg different. Yay!

stripchart(dat$d15N....vs..air. ~ droplevels(dat$treament),
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
           ,ylim=c(-2.2,6)
           ,frame.plot=F
)
boxplot(dat$d15N....vs..air. ~ droplevels(dat$treament), 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
text(c("a","b","b","b"),
     x = c(1,2,3,4),
     y = 4,
     cex = 1.2,
     font = 1,
     xpd = NA)

axis(1,las=2, labels = c("","","",""), at = c(1,2,3,4))
axis(2,las=2, at = c(seq(-2,6, by = 2)))
title(ylab = expression(paste(delta^15, "N")), cex.lab=1.5, line = 2.8)
#with outliers untreated plus and untreated neg different. Yay!

stripchart(dat$ndfa ~ droplevels(dat$treament),
           vertical = TRUE,
           data = dat,
           method = "jitter",
           pch = 20,
           col = add.alpha(colz, alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)
boxplot(dat$ndfa ~ droplevels(dat$treament),
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
text(c("a","b","ab","ab"),
     x = c(1,2,3,4),
     y = 100,
     cex = 1.2,
     font = 1,
     xpd = NA)

axis(1,las=2, labels = c("","","",""), at = c(1,2,3,4))
axis(2,las=2)
title(ylab="NDFA", cex.lab=1.5, line=2.8)

dev.off()

#####################
#Supplemental figure#
#####################

dat <- read.csv("./data/alldata2.csv")
dat$vol <- dat$width*dat$height*dat$Length

controlPloter <- function(ypos = 0.9*par("yaxp")[2],
                          linepos = 0.8*par("yaxp")[2]){
  text("Controls", x = 2.5, 
       y = ypos, 
       cex = 1.5)
  segments(x0 = 0.8, 
           y0 = linepos,
           x1 = 4.2, 
           y1 = linepos)
  segments(x0 = 0.8, 
           y0 = linepos - linepos*(0.02),
           x1 = 0.8, 
           y1 = linepos + linepos*(0.02))
  segments(x0 = 4.2, 
           y0 = linepos - linepos*(0.02),
           x1 = 4.2, 
           y1 = linepos + linepos*(0.02))
  
}

#####################
#SUPPLEMENTAL plot  #
#showing controls   #
#####################

pdf(width=10,height=5, "./visuals/FigS1treatmentEffectsBoxplotsWithCONTROLS.pdf")

par(oma=c(4,4,1,1), mar=c(1,4.5,1,1),mfrow=c(3,3))
stripchart(dat$vol/100^3~dat$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.8),
           cex=1,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)
boxplot(dat$vol/100^3~dat$treament, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","","","","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
axis(1,las=2, labels = c("","","","","","","",""), at = c(1,2,3,4,5,6,7,8))
axis(2,las=2)
title(ylab=expression(paste("Plant volume (",m^3,")",sep="")), cex.lab=1.5, line=3.4,xpd = NA)
controlPloter()

stripchart(dat$X.LeavesCollection~dat$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.8),
           cex=1,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)
boxplot(dat$X.LeavesCollection~dat$treament, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","","","","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
axis(1,las=2, labels = c("","","","","","","",""), at = c(1,2,3,4,5,6,7,8))
axis(2,las=2)
title(ylab="Num. of leaves", cex.lab=1.5, line=2.8)
controlPloter()



stripchart(dat$leafarea_avg~dat$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.8),
           cex=1,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0,8.5)
           ,ylim=c(0.35,.9)
           ,frame.plot=F
)
boxplot(dat$leafarea_avg~dat$treament, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","","","","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
axis(1,las=2, labels = c("","","","","","","",""), at = c(1,2,3,4,5,6,7,8))
axis(2,las=2)
title(ylab=expression(paste("Leaf area (c",m^2,")",sep="")), cex.lab=1.5, line=2.8)
controlPloter(ypos = 0.85, 
              linepos = 0.78)

stripchart(dat$sla_avg~dat$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.8),
           cex=1,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)
boxplot(dat$sla_avg~dat$treament, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","","","","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
axis(1,las=2, labels = c("","","","","","","",""), at = c(1,2,3,4,5,6,7,8))
axis(2,las=2)
title(ylab=expression(paste("SLA",sep="")), cex.lab=1.3, line=2.8)
controlPloter(ypos = 0.19, 
              linepos = 0.17)

dat$X..Swainsonine <- as.character(dat$X..Swainsonine)
dat$X..Swainsonine[which(dat$X..Swainsonine=="#VALUE!")] <- 0

stripchart(as.numeric(dat$X..Swainsonine)~dat$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.8),
           cex=1,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)
boxplot(as.numeric(dat$X..Swainsonine)~dat$treament, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","","","","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
axis(1,las=2, labels = c("","","","","","","",""), at = c(1,2,3,4,5,6,7,8))
axis(2,las=2)
title(ylab=" % Swainsonine", cex.lab=1.5, line=2.8)
#with outliers. untreated plus and untreated neg different. Yay!
controlPloter(ypos = 0.12, 
              linepos = 0.105)

#cor.test(as.numeric(dat$X..Swainsonine),dat$vol)
#i think the lack of swain is because the fungus hadnt grown into the tissue?

stripchart(dat$percC~dat$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.8),
           cex=1,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)
boxplot(dat$percC~dat$treament, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","","","","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
axis(1,las=2, labels = c("","","","","","","",""), at = c(1,2,3,4,5,6,7,8))
axis(2,las=2)
title(ylab="% C", cex.lab=1.5, line=2.8)
controlPloter(ypos = 10, 
              linepos = 9.2)


stripchart(dat$percN~dat$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.8),
           cex=1,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)
boxplot(dat$percN~dat$treament, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","","","","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
axis(1,las=2, labels = c("","","","","","","",""), at = c(1,2,3,4,5,6,7,8))
axis(2,las=2, at = c(seq(0.6,1.8, by = 0.2)))
title(ylab="% N", cex.lab=1.5, line=2.8)

controlPloter()


stripchart(dat$d15N....vs..air.~dat$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.8),
           cex=1,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)
boxplot(dat$d15N....vs..air.~dat$treament, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","","","","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
axis(1,las=2, labels = c("","","","","","","",""), at = c(1,2,3,4,5,6,7,8))
axis(2,las=2, at = c(seq(-2,6, by = 2)))
title(ylab = expression(paste(delta^15, "N")), cex.lab=1.5, line = 2.8)

controlPloter()
stripchart(dat$ndfa~dat$treament,
           vertical = TRUE,
           data = dat,
           method = "jitter",
           pch = 20,
           col = add.alpha(colz, alpha=0.8),
           cex=1,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)
boxplot(dat$ndfa~dat$treament,
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","","","","","",""),
        axes=F,
        col=add.alpha(colz, alph = 0.6),
        ylab=""
)
axis(1,las=2, labels = c("","","","","","","",""), at = c(1,2,3,4,5,6,7,8))
axis(2,las=2)
title(ylab="NDFA", cex.lab=1.5, line=2.8)
controlPloter(ypos = 17, 
              linepos = 8)

dev.off()
