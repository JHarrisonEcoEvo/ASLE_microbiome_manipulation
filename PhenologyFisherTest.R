dat <- read.csv("data/alldata2.csv")
flowering <- table(dat$Flowering.,dat$treament)[,5:8]

fisher.test(cbind(c(9,32),
                  c(5,35)))
fisher.test(cbind(c(17,70),
                  c(11,70)))
#Aphids
dat2 <- read.csv("./data/size2ndyear2018.csv")
dat3 <- merge(dat, dat2, 
      by.y = "plant",
      by.x = "plant")

table(dat3$aphids,dat3$treament)

