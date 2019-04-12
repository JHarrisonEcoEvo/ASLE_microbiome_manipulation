#J. Harrison

set.seed(666) 
##########################
# Load packages and data #
##########################

library(vegan)

load("./withAfulva_BayesianRelativeAbundanceModelingFUNGInoplant_with_p_theta.RData")

mlEst <- matrix(nrow = 146,
                ncol = 13)
for(i in 1:146){
  mlEst[i,1:11] <- 
    apply(results_with_Afulva[[1]]$p[i,,,], 1, mean)
   
  mlEst[i, 12] <- exp(vegan::diversity(mlEst[i,1:11],
                                  index = "shannon"))
  mlEst[i, 13] <- 1 / vegan::diversity(mlEst[i,1:11],
                                   index = "simpson")
}
mlEst <- data.frame(orderModeled, mlEst) 
names(mlEst) <- c("orderModeled", names(dat), "shannon", "simpson")

write.csv(mlEst, "./data/ML_fungi_ALFU_withShannonSimpson.csv")

load("withOUTAfulva_BayesianRelativeAbundanceModelingFUNGInoplant_with_p_theta.RData")

mlEst <- matrix(nrow = 145,
                ncol = 10)
for(i in 1:145){
  mlEst[i,1:8] <- 
    apply(results_withOut_Afulva[[1]]$p[i,,,], 1, mean)
   
  mlEst[i, 9] <- exp(vegan::diversity(mlEst[i,1:8],
                                  index = "shannon"))
  mlEst[i, 10] <- 1 / vegan::diversity(mlEst[i,1:8],
                                   index = "simpson")
}
mlEst <- data.frame(orderModeled_less, mlEst) 

names(mlEst) <- c("samps", names(dat), "shannon", "simpson")

write.csv(mlEst, "./data/ML_fungi_noALFU_withShannonSimpson.csv")

#Do for bacteria
rm(list=ls())
load("./BayesianRelativeAbundanceModelingBACTERIA.RData")

mlEst <- matrix(nrow = 145,
                ncol = 55)
for(i in 1:145){
  mlEst[i,1:53] <- apply(results[[1]]$p[i,,,], 1, mean)

  mlEst[i, 54] <- mean(vegan::diversity(mlEst[i,1:53] ,
                                   index = "shannon"))
  mlEst[i, 55] <- mean(vegan::diversity(mlEst[i,1:53] ,
                                   index = "simpson"))
}
mlEst <- data.frame(orderModeled, mlEst) 

names(mlEst) <- c("samps", names(dat), "shannon", "simpson")

write.csv(mlEst, "./data/ML_bacteria_withShannonSimpson.csv")
