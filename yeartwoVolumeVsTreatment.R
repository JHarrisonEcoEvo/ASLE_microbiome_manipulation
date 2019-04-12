
dat <- read.csv("./data/alldata2.csv", header = T)
size2018 <- read.csv("./data/size2ndyear2018.csv")
size2018$vol <- size2018$height * size2018$width * size2018$length
dat2 <- merge(dat, size2018, by.y = "plant", by.x = "plant", all.x = T)
dat <- dat2

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#Function to use Anova to test for treatment effects
#
# anovator <- function(x, y){
#   aovout <- aov(x~y)
#   #print(TukeyHSD(aovout))
#   print(car::leveneTest(x,y))
#   #remove outliers and do over
#   outliers <- boxplot.stats(x)$out
#   x[which(x %in% outliers)] <- NA
#   aovout <- aov(x~y)
#   print(TukeyHSD(aovout))
#   print(car::leveneTest(x,y))
# }

colz <- c("gray", "dark gray")

#dat <- dat[-grep("control",dat$treament),]

#plot
pdf(width=10,height=8, "./visuals/FigSX_yeartwoVolumeVsTreatment.pdf")

par(oma=c(8,4,0,1), mar=c(1,4.5,1,1))
stripchart(dat$vol.y~dat$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colz, alpha=0.8),
           cex=2.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,8.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)
boxplot(dat$vol.y~dat$treament, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        names= c("","","","","","","",""),
        axes=F,
        col=c(add.alpha("light gray", alpha=0.6), 
              add.alpha("light gray", alpha=0.6), 
              NA, 
              NA),
        ylab=""
)
axis(1,
     las=2, 
     cex.axis = 1.2, 
     labels = c("Control -",
                         "Control +",
                          "Treated, Control - ",
                          "Treated, Control + ",
                         "Treated -",
                         "Treated +",
                         "Untreated -",
                         "Untreated + "), 
     at = c(1,2,3,4,5,6,7,8))
axis(2,las=2)
title(ylab=expression(paste("Plant volume (c",m^3,") year two",sep="")), cex.lab=2, line=3.7,xpd = NA)

dev.off()

