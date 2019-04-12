load("./withOUTAfulva_BayesianRelativeAbundanceModelingFUNGInoplant_with_p_theta.RData")

#Function to compute diversity indices for each MCMC sample
rich_div_est <- function(y, index) {
  #richest <- NA
  Shannon <- NA
  Simpson <- NA
  #
  k <- 1
  for (i in 1:1000) {
    
    #Indexing is treatment, taxon, mcmc sample, chain
    props1 <- y[index, , i, 1]
    props2 <- y[index, , i, 2]
    
    Shannon[k] <- exp(vegan::diversity(props1,
                                   index = "shannon"))
    Simpson[k] <- 1 / vegan::diversity(props1,
                                   index = "simpson")
    k <- k + 1
    
    Shannon[k] <- exp(vegan::diversity(props2,
                                   index = "shannon"))
    Simpson[k] <- 1 / vegan::diversity(props2,
                                   index = "simpson")
    k <- k + 1
  }
  
  return(list(
         shannon = Shannon, 
         simpson = Simpson))
}

#To determine if Shannon's differs between treatment groups
#treatment group order: 
unique(treatment)

quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi, 1)[[1]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,2)[[1]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,3)[[1]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,4)[[1]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[1]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[1]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[1]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[1]], c(0.025, 0.975))

quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[1]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[1]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[1]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[1]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,1)[[1]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,3)[[1]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,2)[[1]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,4)[[1]], c(0.025, 0.975))

#Simpsons

quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,1)[[2]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,2)[[2]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,3)[[2]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,4)[[2]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[2]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[2]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[2]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[2]], c(0.025, 0.975))

quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[2]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[2]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[2]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[2]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,1)[[2]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,3)[[2]], c(0.025, 0.975))
quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,2)[[2]] - rich_div_est(y = results_withOut_Afulva[[1]]$pi,4)[[2]], c(0.025, 0.975))


###############
#Fig for paper#
###############

#Function to make credible interval lines
limner <- function(coords, index) {
  #add vertical lines
  lines(
    x = c(index, index),
    y = c(coords[[1]], coords[[2]]),
    lty = 1,
    lwd = 2
  )
  #add horizontal end lines
  lines(
    x = c(index - 0.05,
          index + 0.05),
    y = c(coords[[1]],
          coords[[1]]),
    lty = 1,
    lwd = 2
  )
  lines(
    x = c(index - 0.05,
          index + 0.05),
    y = c(coords[[2]],
          coords[[2]]),
    lty = 1,
    lwd = 2
  )
  points(index, coords[[3]], 
         pch = 23,
         col = "black",
         bg = "black",
         cex = 2)
}

pdf(width=5,height=6, file = "./visuals/Supplemental_RichDiv_vs_treatment_withControl_ITS.pdf")
par(mfrow = c(2,1), 
    oma=c(6,0,2,0),
    mar = c(1,5,0,0),
    xpd = NA)
plot(NULL,
     xlim = c(0.8,4.7),
     ylim = c(1.1,1.3),
     ylab = "Shannon's",
     xlab = "", 
     las = 2, 
     xaxt = "n", 
     bty = "n",
     cex.lab = 1.3, 
     mar=c(0,0,0,0))

coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[1]]))
limner(coords, 1)

coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[1]]))
limner(coords, 1.5)
coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[1]]))
limner(coords, 2)
coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[1]]))
limner(coords, 2.5)

#controls
coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,3)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,3)[[1]]))
limner(coords, 3)
coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,4)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,4)[[1]]))
limner(coords, 3.5)
coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,1)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,1)[[1]]))
limner(coords, 4)
coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,2)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,2)[[1]]))
limner(coords, 4.5)

plot(NULL,
     xlim = c(0.8,4.7),
     ylim = c(0,30),
     ylab = "Simpson's",
     xlab = "", 
     las = 2, 
     xaxt = "n", 
     yaxt = "n", 
     bty = "n",
     cex.lab = 1.3, 
     mar=c(0,0,0,0))
axis(side=2, 
     at = seq(0,30, by = 5),
     labels = c(seq(0,30, by = 5)),
     las = 1
)
axis(side=1, 
     at = seq(1,4.5, by = 0.5),
     labels = F
)
text(seq(1,4.5, by = 0.5),
     par("usr")[3]-5.5,
     #pos = 2,
     xpd = NA, 
     srt = 45,
     labels = c( expression(paste("No ", italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste("No ", italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste("No ", italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 
                 expression(paste("No ", italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = ""))
     )
)

text(3.7, 
     labels = "Controls",
     xpd= NA,
     par("usr")[3]-12.5
     )
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
rect(col=add.alpha("grey", alpha = 0.2),
     xleft = 0.8,
     xright = 1.8,
     ybottom = 0,
     ytop = 66)
text("Inoculum\ntreated",
     x = 1.3, 
     y = 70)

rect(col=add.alpha("grey", alpha = 0.2),
     xleft = 2.8,
     xright = 3.8,
     ybottom = 0,
     ytop = 66)
text("Inoculum\ntreated",
     x = 3.3, 
     y = 70)

coords <-  c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[2]], c(0.025, 0.975)), 
             mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[2]]))
limner(coords, 1)
coords <-  c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[2]], c(0.025, 0.975)), 
             mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[2]]))
limner(coords, 1.5)
coords <-  c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[2]], c(0.025, 0.975)), 
             mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[2]]))
limner(coords, 2)
coords <-  c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[2]], c(0.025, 0.975)), 
             mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[2]]))
limner(coords, 2.5)

coords <-  c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,3)[[2]], c(0.025, 0.975)), 
             mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,3)[[2]]))
limner(coords, 3)
coords <-  c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,4)[[2]], c(0.025, 0.975)), 
             mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,4)[[2]]))
limner(coords, 3.5)
coords <-  c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,1)[[2]], c(0.025, 0.975)), 
             mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,1)[[2]]))
limner(coords, 4)
coords <-  c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,2)[[2]], c(0.025, 0.975)), 
             mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,2)[[2]]))
limner(coords, 4.5)

dev.off()

####################
#No control plotting
####################
load("./withOUTAfulva_BayesianRelativeAbundanceModelingFUNGInoplant_with_p_theta.RData")

pdf(width=5,height=8, file = "./visuals/Fig2_Div_vs_treatment.pdf")
par(mfrow = c(4,1), 
    oma=c(8,3,4,0),
    mar = c(1,5,0,0),
    xpd = NA)
plot(NULL,
     xlim = c(0.8,2.8),
     ylim = c(1.1,1.3),
     ylab = "Shannon's",
     xlab = "", 
     las = 2, 
     #line = 3.5,
     xaxt = "n", 
     #yaxt = "n", 
     bty = "n",
     cex.lab = 1.5, 
     mar=c(0,0,0,0))

coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[1]]))
limner(coords, 1)

coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[1]]))
limner(coords, 1.5)
coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[1]]))
limner(coords, 2)
coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[1]]))
limner(coords, 2.5)

mtext(side = 2, 
      "Fungi", 
      cex = 1.5, 
      padj = -3.1,
      adj = -0.3,
      line = 2,
      xpd = NA)


plot(NULL,
     xlim = c(0.8,2.7),
     ylim = c(10,30),
     ylab = "Simpson's",
     xlab = "", 
     las = 2, 
    # line = 3.5,
     xaxt = "n", 
     yaxt = "n", 
     bty = "n",
     cex.lab = 1.5, 
     mar=c(0,0,0,0))

axis(side=2, 
     at = seq(10,30, by = 5),
     labels = c(seq(10,30, by = 5)),
     las = 1
)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[2]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,5)[[2]]))
limner(coords, 1)

coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[2]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,6)[[2]]))
limner(coords, 1.5)
coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[2]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,7)[[2]]))
limner(coords, 2)
coords <- c(quantile(rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[2]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results_withOut_Afulva[[1]]$pi,8)[[2]]))
limner(coords, 2.5)

###############
#Plot bacteria#
###############

load("./BayesianRelativeAbundanceModelingBACTERIA.RData")

# pdf(width=5,height=6, file = "./visuals/Fig2b_RichDiv_vs_treatment16s.pdf")
# par(mfrow = c(2,1), 
#     oma=c(6,0,2.5,0),
#     mar = c(1,5,0,0),
#     xpd = NA)
plot(NULL,
     xlim = c(0.8,2.8),
     ylim = c(45,55),
     ylab = "Shannon's",
     xlab = "", 
     las = 2, 
    # line = 3.5,
     xaxt = "n", 
     bty = "n",
     cex.lab = 1.5, 
     mar=c(0,0,0,0))
# axis(side=1, 
#      at = seq(1,2.5, by = 0.5),
#      labels = c(expression(paste(italic("A. fulva"), "; inoculum", sep = "")),
#                 expression(paste("No ", italic("A. fulva"), "; inoculum", sep = "")),
#                 expression(paste(italic("A. fulva"), "; no inoculum", sep = "")),
#                 expression(paste("No ", italic("A. fulva"), "; no inoculum", sep = ""))
#                 ),
#      las = 2
# )
coords <- c(quantile(rich_div_est(y = results[[1]]$pi,5)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results[[1]]$pi,5)[[1]]))
limner(coords, 1)

coords <- c(quantile(rich_div_est(y = results[[1]]$pi,6)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results[[1]]$pi,6)[[1]]))
limner(coords, 1.5)
coords <- c(quantile(rich_div_est(y = results[[1]]$pi,7)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results[[1]]$pi,7)[[1]]))
limner(coords, 2)
coords <- c(quantile(rich_div_est(y = results[[1]]$pi,8)[[1]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results[[1]]$pi,8)[[1]]))
limner(coords, 2.5)

mtext(side = 2, 
      "Bacteria", 
      cex = 1.5, 
      padj = -3.1,
      line = 2,
      adj = -0.5,
      xpd = NA)


plot(NULL,
     xlim = c(0.8,2.7),
     ylim = c(1.015,1.03),
     ylab = "Simpson's",
     xlab = "", 
     las = 2, 
     xaxt = "n", 
     yaxt = "n", 
     bty = "n",
     cex.lab = 1.5, 
     line = 3.5,
     mar=c(0,0,0,0))
axis(side=2, 
     at = c(seq(1.015,1.03, by = 0.005)),
     labels = c(seq(1.015,1.03, by = 0.005)),
     las = 1
)
axis(side=1, 
     at = seq(1, 2.5, by = 0.5),
     labels = F
)
text(seq(1,2.5, by = 0.5),
     par("usr")[3] - 0.005,
     #pos = 2,
     xpd = NA, 
     srt = 45,
     cex = 2,
     labels = c(expression(paste(italic("A. fulva"),
                                 #"; inoculum", 
                                 sep = "")),
                expression(paste("No ", italic("A. fulva"), 
                                 #"; inoculum", 
                                 sep = "")),
                expression(paste(italic("A. fulva"), 
                                 #"; no inoculum",
                                 sep = "")),
                expression(paste("No ", italic("A. fulva"), 
                                 #"; no inoculum", 
                                 sep = ""))
     )
)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

rect(col=add.alpha("grey", 
                   alpha = 0.2),
     xleft = 0.8,
     xright = 1.7,
     xpd = NA,
     ybottom = 1.015,
     ytop = 1.0822)
text("Inoculum\ntreated",
     x = 1.3,
     y = 1.086,
     cex = 1.8)

coords <- c(quantile(rich_div_est(y = results[[1]]$pi,5)[[2]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results[[1]]$pi,5)[[2]]))
limner(coords, 1)

coords <- c(quantile(rich_div_est(y = results[[1]]$pi,6)[[2]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results[[1]]$pi,6)[[2]]))
limner(coords, 1.5)
coords <- c(quantile(rich_div_est(y = results[[1]]$pi,7)[[2]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results[[1]]$pi,7)[[2]]))
limner(coords, 2)
coords <- c(quantile(rich_div_est(y = results[[1]]$pi,8)[[2]], c(0.025, 0.975)), 
            mean(rich_div_est(y = results[[1]]$pi,8)[[2]]))
limner(coords, 2.5)

dev.off()


