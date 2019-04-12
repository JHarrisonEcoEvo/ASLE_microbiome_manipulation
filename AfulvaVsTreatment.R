load("./withAfulva_BayesianRelativeAbundanceModelingFUNGInoplant_with_p_theta.RData")

#to see OTUs
names(dat)

#plot A. fulva for untreated plants, 
#calc averages
apply(results_with_Afulva[[1]]$pi[1,6,,1:2], MARGIN = 2, FUN=mean)
apply(results_with_Afulva[[1]]$pi[1,9,,1:2], MARGIN = 2, FUN=mean)
apply(results_with_Afulva[[1]]$pi[2,6,,1:2], MARGIN = 2, FUN=mean)
apply(results_with_Afulva[[1]]$pi[2,9,,1:2], MARGIN = 2, FUN=mean)

pdf(width=8, height = 8, file="./visuals/AfulvaDensityPlots.pdf")
#Untreated is 7(-) and 8 (+)
plot(density(results_with_Afulva[[1]]$pi[7,9,,1:2]),
     col="green",
     lwd=2.5, 
     ylim = c(0,1000), 
     xlim= c(0.0,0.15),
     main="",
     xlab="Proportion of total fungal reads",
     ylab = "MCMC sample density")

lines(density(results_with_Afulva[[1]]$pi[7,6,,1:2]),
      col="dark green",
      lwd=2.5)
lines(density(results_with_Afulva[[1]]$pi[8,9,,1:2]),
      col="cadet blue",
      lwd=2.5)
lines(density(results_with_Afulva[[1]]$pi[8,6,,1:2]),
      col="dark blue",
      lwd=2.5)
#Treated is 5 (-) and 6 (+)
lines(density(results_with_Afulva[[1]]$pi[5,9,,1:2]),
      col=" green",
      lwd=2.5,
      lty=3)
lines(density(results_with_Afulva[[1]]$pi[5,6,,1:2]),
      col="dark green",
      lwd=2.5,
      lty=3)
lines(density(results_with_Afulva[[1]]$pi[6,9,,1:2]),
      col="cadet blue",
      lwd=2.5,
            lty=3)
lines(density(results_with_Afulva[[1]]$pi[6,6,,1:2]),
      col="dark blue",
      lwd=2.5,
            lty=3)

legend("topright", legend = c("A. fulva genotype 1 - No inoculum; seed intact",
         "A. fulva genotype 2 - No inoculum; seed intact",
         "A. fulva genotype 1 - No inoculum; embryo excised",
         "A. fulva genotype 2 - No inoculum; embryo excised",
         "A. fulva genotype 1 - Inoculum applied; seed intact",
         "A. fulva genotype 2 - Inoculum applied; seed intact",
         "A. fulva genotype 1 - Inoculum applied; embryo excised",
         "A. fulva genotype 2 - Inoculum applied; embryo excised"),
         lty=c(1,1,1,1,3,3,3,3), 
       lwd=2.5, 
       bty = "n",
       col=c("dark blue",
             "cadet blue",
             "dark green",
             "green",
             "dark blue",
             "cadet blue",
             "dark green",
             "green"))
dev.off()

#sanity check that I am plotting the correct thing
#may not work without object from modeling script
# test <- data.frame(cbind(newotus$Otu9/rowSums(newotus[,2:27]), treat))
# aggregate(test$V1~test$treat, FUN=mean)


pdf(width=8, height = 8, file="./visuals/AfulvaDensityPlotsCONTROLS.pdf")
#Untreated is 7(-) and 8 (+)
plot(density(results_with_Afulva[[1]]$pi[1,9,,1:2]),
     col="green",
     lwd=2.5, 
     ylim = c(0,1000), 
     xlim= c(0.0,0.3),
     main="",
     xlab="Proportion of total fungal reads",
     ylab = "MCMC sample density")

lines(density(results_with_Afulva[[1]]$pi[1,6,,1:2]),
      col="dark green",
      lwd=2.5)
lines(density(results_with_Afulva[[1]]$pi[2,9,,1:2]),
      col="cadet blue",
      lwd=2.5)
lines(density(results_with_Afulva[[1]]$pi[2,6,,1:2]),
      col="dark blue",
      lwd=2.5)
#Treated is 5 (-) and 6 (+)
lines(density(results_with_Afulva[[1]]$pi[4,9,,1:2]),
      col=" green",
      lwd=2.5,
      lty=3)
lines(density(results_with_Afulva[[1]]$pi[4,6,,1:2]),
      col="dark green",
      lwd=2.5,
      lty=3)
lines(density(results_with_Afulva[[1]]$pi[3,9,,1:2]),
      col="cadet blue",
      lwd=2.5,
            lty=3)
lines(density(results_with_Afulva[[1]]$pi[3,6,,1:2]),
      col="dark blue",
      lwd=2.5,
            lty=3)

legend("topright", legend = c("A. fulva genotype 1 - No inoculum; agar +",
         "A. fulva genotype 2 - No inoculum; agar +",
         "A. fulva genotype 1 - No inoculum; agar -",
         "A. fulva genotype 2 - No inoculum; agar -",
         "A. fulva genotype 1 - Inoculum applied; agar +",
         "A. fulva genotype 2 - Inoculum applied; agar +",
         "A. fulva genotype 1 - Inoculum applied; agar -",
         "A. fulva genotype 2 - Inoculum applied; agar -"),
         lty=c(1,1,1,1,3,3,3,3), 
       lwd=2.5, 
       bty = "n",
       col=c("dark blue",
             "cadet blue",
             "dark green",
             "green",
             "dark blue",
             "cadet blue",
             "dark green",
             "green"))
dev.off()
