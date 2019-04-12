#J. Harrison
#Nov 11, 2018
set.seed(666)
options(scipen=99)

#This program calculates variation in technical replicates and examines
#how well we capture the expected biodiversity in the mock community

####################################
# Load data, packages, and wrangle #
####################################

its_otus <- read.csv( "./data/fungiOTUtableZotus.csv")

its_otus$newotus.sampleName <- as.character(its_otus$newotus.sampleName)

techs <- its_otus[grep("\\d+[ab]", its_otus$newotus.sampleName),]

techs <- techs[order(techs$newotus.sampleName),]

#convert to proportions

props <- data.frame(matrix(nrow = 14,
                ncol = 26))

for(i in 1:14){
  x <- techs[i,2:length(techs)]

  props[i,] <- x/sum(x)
}

#great success
data.frame(techs$newotus.sampleName,
           props)

############################
# Check out mock community #
############################

otus <- read.table("./data/16s_data/otuTableZotus_raw (1).txt",
                   header= T)
mock <- otus[which(otus[,98:99] > 0), c(1,98:99)]
mock <- mock[1:17,]; mock


otus <- read.table("./data/16s_data/otuTable97otus_raw (1).txt",
                   header= T)
mock <- otus[which(otus[,98:99] > 0), c(1,98:99)]
mock <- mock[1:17,]

#for fungi

otus <- read.table("./data/ITSdata/otuTableZotus_raw (1).txt",
                   header= T)
mock <- otus[which(otus[,194] > 0), c(1,194)]
mock


otus <- read.table("./data/ITSdata/otuTable97otus_rawITS (1).txt",
                   header= T)
mock <- otus[which(otus[,194] > 0), c(1,194)]
mock
