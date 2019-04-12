
load("./BayesianRelativeAbundanceModelingBACTERIA.RData")

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

  
#untreated minus  and plus
diff_out <- matrix(nrow = 8, 
                   ncol = 8)

for(i in 1:8){
  for(j in 1:8){
    diffs <-
      apply(results[[1]]$pi[i, , , 1:2] - results[[1]]$pi[j, , , 1:2],
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

