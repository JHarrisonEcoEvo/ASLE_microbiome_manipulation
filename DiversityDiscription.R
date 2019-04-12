#read and wrangle a bit
dat <- read.csv("./data/fungiOTUtableZotus.csv")
dat$samps <- as.character(dat$samps)
dim(dat)

#Import treatment information
samps <- read.csv("./data/planting_order.csv")
dim(dat)

#Merge datasets
dat <- merge(dat, 
      samps,
      by.x = "samps",
      by.y = "rns",
      all.x = T)

dat$vol = dat$height*dat$width*dat$length
dat <- dat[order(dat$height*dat$width*dat$length),]

table(cbind(dat$treament, dat$vol)[143:172,])
dat <- dat[order(dat$treatment),]

dim(dat)

#clean up a bit
#Remove the control, tech replicates, inoculum mixture, culture
#etc. 

#Also remove samples with no reads, if there are any
rmv <- which(rowSums(dat[,2:22]) == 0)
if(length(rmv) > 0){
  dat <- dat[-rmv,]
}
dat <- dat[order(dat$treament),]
treatment <- dat$treament

#Remove the non ESV count columns
dat <- dat[,-c(1, 23:length(dat))]
dim(dat) #21 taxa, as expected, however some of these have zero counts in the data.
#That is because these were in inoculum mixture, cultures, or one fo the negative/positive controls
#We remove those here
dat <- dat[,-which(colSums(dat) == 0)]
sum(colSums(dat))
dat <- t(dat)
dat$otus <- row.names(dat)

tax <- read.table("./data/ITSdata/combinedTaxonomy_zotus_ITS (1).txt",
                  fill = T)
tax$V1 <-  as.character(tax$V1)

tax[which(tax$V1 %in% dat$otus),]

#################
#Do for bacteria#
#################

#read and wrangle a bit
dat <- read.csv("./data/bacteriaOTUtableZotus.csv")
dat$samps <- as.character(dat$samps)

#Import treatment information
samps <- read.csv("./data/planting_order.csv")
dim(dat)

#Merge datasets
dat <- merge(dat, 
      samps,
      by.x = "samps",
      by.y = "rns",
      all.x = T)

#Also remove the odd sample that had no reads
rmv <- which(rowSums(dat[,2:113]) == 0)

if(length(rmv) > 0){
  dat <- dat[-rmv,]
}

dat <- dat[order(dat$treament),]
dim(dat)

#clean up a little
dat <- dat[-(174:200),]
dat <- dat[,-1]
treatment <- dat$treament

dat <- dat[,1:112] #DOUBLECHECK THIS DURING EACH RUN

#Remove taxa that are not present in samples
dat <- dat[,-which(colSums(dat) == 0)]
sum(colSums(dat))

tax <- read.table("./data/16s_data/combinedTaxonomy.txt",
                  fill = T)

tax$V1 <- as.character(tax$V1)

dat <- tax[which(tax$V1 %in% names(dat)),]

tax <- strsplit(x = as.character(dat$V3), split = ",")

table(sapply(tax,"[", 2))

##################################################
#Determine what from inoculum made it into plants#
##################################################

tax <- read.table("./data/ITSdata/combinedTaxonomy_zotus_ITS (1).txt",
                  fill = T)
tax$V1 <-  as.character(tax$V1)

dat <- read.csv("./data/fungiOTUtableZotus.csv")
dat$samps <- as.character(dat$samps)
dim(dat)

#Import treatment information
samps <- read.csv("./data/planting_order.csv")
dim(dat)

#Merge datasets
dat <- merge(dat, 
      samps,
      by.x = "samps",
      by.y = "rns",
      all.x = T)

dat <- dat[order(dat$treament),]
dim(dat)

dat$samps 

#inoculum samples
dat[187:190,]

colSums(dat[187:190,2:22])
length(which(colSums(dat[187:190,2:22]) > 0 ))

colSums(dat[1:186,2:22])

tax[which(tax$V1 == "Otu11"),]

#for bacteria

#read and wrangle a bit
dat <- read.csv("./data/bacteriaOTUtableZotusALLSEQS.csv")
dat$samps <- as.character(dat$samps)

germies <- dat[grep("Inoculum", dat$samps),]

colSums(germies[,2:length(germies)])
length(which(colSums(germies[,2:length(germies)]) > 0 ))

insamples <- dat[-grep("[A-Za-z]", dat$samps),]
insamples <- insamples[1:172,]

colSums(insamples[,2:length(insamples)]) 

x <- cbind(colSums(germies[,2:length(germies)]),
      colSums(insamples[,2:length(insamples)]))

x[which(x[,1] > 0 ), 1] <- 1
x[which(x[,2] > 0 ), 2] <- 1
length(which(rowSums(x) == 2 ))

#I did a sanity check and reloaded data to make sure this was correct.

row.names(x)[which(rowSums(x) == 2 )]

tax <- read.table("./data/16s_data/combinedTaxonomy.txt",
                  fill = T)

tax$V1 <- as.character(tax$V1)

inboth <- tax[tax$V1 %in% row.names(x)[which(rowSums(x) == 2 )],]

inboth2 <- strsplit(x = as.character(inboth$V3), split = ",")

table(sapply(inboth2,"[", 2))
