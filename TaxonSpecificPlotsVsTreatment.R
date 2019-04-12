
load("./withOUTAfulva_BayesianRelativeAbundanceModelingFUNGInoplant_with_p_theta.RData")

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#for fungi
taxa <- names(dat)
smallfry <- which(colSums(dat) < 50)

#Remove less abundant taxa
if(length(smallfry) > 0 ){
  taxa2 <- taxa[-smallfry]
  results_withOut_Afulva[[1]]$p <- results_withOut_Afulva[[1]]$p[,-smallfry , , 1:2] 
}

#Calculate max likelihood of each taxon in each replicate
#Only doing this for the abundant taxa
mls <- matrix(nrow = 145,
              ncol = length(taxa2))
#for fungi
for(i in 1:145){
  mls[i,] <- apply(results_withOut_Afulva[[1]]$p[i,,,1:2], MARGIN = 1, FUN = mean)
}

forplot <- data.frame(mls, treatment)

########
# Plot #
########

pdf(width = 8, height = 8, file = "./visuals/fungalTaxonSpecific.pdf")
par(mfrow = c(2,2),
    mar = c(1,3,0,2),
    oma = c(7, 4,1,0))
for(i in 1:4){
  
  #if else is to put different y axis range on first plot
  if(i == 1){
   boxplot(forplot[,i] ~ forplot$treatment,
        las = 2,
        xaxt = "n",
        yaxt = "n",
        outline = F,
        ylim = c(0.9, round(max(forplot[,i]),3)),
        col=c(add.alpha("light gray", alpha=0.6), 
                add.alpha("light gray", alpha=0.6), 
                         NA, 
                         NA))
    axis(side = 2,
         at =  c(0.9, 0.95, 0.993),
         labels = c(0.9, 0.95,  0.993),
         las = 2
    )
  }else{
  
    boxplot(forplot[,i] ~ forplot$treatment,
        las = 2,
        xaxt = "n",
        yaxt = "n",
        outline = F,
        ylim = c(0, round(max(forplot[,i]),3)),
        col=c(add.alpha("light gray", alpha=0.6), 
                add.alpha("light gray", alpha=0.6), 
                         NA, 
                         NA))
    axis(side = 2,
         at = c(0,
               round(max(forplot[,i]) / 2, 3),
                round(max(forplot[,i]),3)
         ),
         labels = c(0,
                round(max(forplot[,i]) / 2, 3),
               round( max(forplot[,i]),3)
         ),
         las = 2
    )
 
  if(i %in% c(3,4)){
    axis(side = 1, 
         las = 2, 
         at = seq(1,8, by = 1),
         labels = c("Control -", 
                    "Control +", 
                    "Treated, Control - ",
                    "Treated, Control + ",
                    "Treated - ",
                    "Treated +",
                    "Untreated -",
                    "Untreated +"
         ),
         cex = 1.5)
  }
}

# add taxon titles, this reflects what will be correct after rerunning models. Nov 21
if(i == 1){
  text(substitute(italic("L. taurica")), #Otu4
        x = 1.8,
        y = 0.9,
       cex = 1.5,
        xpd = NA)
}
if(i == 2){
  text(substitute(italic("L. taurica")),#Otu58
        x = 1.8,
        y = 0,
       cex = 1.5,
        xpd = NA)
}
if(i == 3){
  text(substitute(italic("L. taurica")),#Otu77
        x = 1.8,
        y = 0,
       cex = 1.5,
        xpd = NA)
}
if(i == 4){
  text(substitute(italic("Phyllactinia")),#Otu70
        x = 2,
        y = 0.086,
       cex = 1.5,
        xpd = NA)
}
}
dev.off()


