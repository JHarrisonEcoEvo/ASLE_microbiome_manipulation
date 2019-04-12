
load("./withOUTAfulva_BayesianRelativeAbundanceModelingFUNGInoplant_with_p_theta.RData")
library(compositions)

# #for bacteria
# taxa <- names(dat)
# smallfry <- which(colSums(dat) < 10)
# 
# if(length(smallfry) > 0 ){
#   taxa2 <- taxa[-smallfry]
#   results[[1]]$pi <- results[[1]]$pi[,-smallfry , , 1:2] 
# }

########################
# Perform calculations #
########################

#1"controlNeg"         2 "controlplus"        3 "treated_controlneg"  4"treated_controlplus"
#[5] "treated_neg"       6  "treated_plus"      7  "untreated_neg"      8 "untreated_plus"     

# clr <- function(x){
#   log(x / (prod(x))^(1/length(x)))
# }

diff_out <- matrix(nrow = 8, 
                   ncol = 8)

for(i in 1:8){
  for(j in 1:8){
    out_i <- t(apply(cbind(results_with_Afulva[[1]]$pi[i,,,1],
                           results_with_Afulva[[1]]$pi[i,,,2]),
                     2, 
                     FUN = clr))
    out_j <- t(apply(cbind(results_with_Afulva[[1]]$pi[j,,,1],
                           results_with_Afulva[[1]]$pi[j,,,2]),
                     2, 
                     FUN = clr))
    #calculate credible intervals of differences
    diffs <- apply(out_i - out_j, 2, quantile, 
                   probs = c(0.025,0.975))
    
    print(apply(out_i - out_j, 2, mean))
    
    sigtaxa <- c(which(diffs[1,] < 0 & diffs[2,] < 0),
                 which(diffs[1,]  > 0 & diffs[2,]  > 0))
    print(paste(i, j , sep = "_"))
    
    print(sigtaxa)
    diff_out[i,j] <- length(sigtaxa)
    
  }
}
diff_out

#untreated minus  and plus
diff_out <- matrix(nrow = 8, 
                   ncol = 8)

for(i in 1:8){
  for(j in 1:8){
    diffs <-
      apply(results_withOut_Afulva[[1]]$pi[i, , , 1:2] - results_withOut_Afulva[[1]]$pi[j, , , 1:2],
            1,
            quantile,
            probs = c(0.025, 0.975))
    
    sigtaxa <- c(which(diffs[1,] < 0 & diffs[2,] < 0),
                 which(diffs[1,]  > 0 & diffs[2,]  > 0))
    print(paste(i, j , sep = "_"))
    
    print(sigtaxa)
    diff_out[i,j] <- length(sigtaxa)
  }
}
diff_out
unique(treatment)
names(dat)


