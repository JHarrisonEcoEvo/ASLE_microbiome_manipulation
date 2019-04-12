plantorder = read.csv("./data/planting_order.csv")
table(plantorder$treament)

#Add sla and leaf area
sla = read.csv("./data/sla.csv")

alldat <- merge(sla, plantorder, by.x = "Plant", by.y="rns", all=T)
dim(alldat)
dim(plantorder)

alldat$leafarea_avg = rowMeans(data.frame(alldat[,5:7], na.rm=T))

alldat$sla_1 = alldat$a1/alldat$m1
alldat$sla_2 = alldat$a2/alldat$m2
alldat$sla_3 = alldat$a3/alldat$m3
alldat$sla_4 = alldat$a4/alldat$m4

alldat$sla_avg = rowMeans(alldat[,27:30],na.rm=T)

#Add CN data
cn <- read.csv("./data/CNmeasurements.csv")
key <- read.csv("./data/CNkey.csv")
names(key)[which(names(key)=="Sample")] <- "Sample2"

#first replace duplicate samples with their mean (tech.replicates)
key$well <- paste("T", key$Tray, "-", key$Row., sep="")
dat <- merge(cn, key, by.x="Sample", by.y = "well")
dat$percN <- dat$wt...N/dat$Weight.mg.
dat$percC <- dat$wt...C/dat$Weight.mg.

dupes <- dat[grep("*[ab]", dat$Sample2),]

dupes$samp <- gsub("(\\d+)[ab]","\\1", dupes$Sample2)

d15N<- NA
whichDupe <- NA
l<- 1
for(i in unique(dupes$samp)){
  d15N[l] <- mean(dupes[which(dupes$samp==i),2])
  whichDupe[l]<- i
  l<-l+1
}
c13<- NA
l<- 1
for(i in unique(dupes$samp)){
  c13[l] <- mean(dupes[which(dupes$samp==i),2])
  l<-l+1
}

percN<- NA
l<- 1
for(i in unique(dupes$samp)){
  percN[l] <- mean(dupes[which(dupes$samp==i),13])
  l<-l+1
}
percC<- NA
l<- 1
for(i in unique(dupes$samp)){
  percC[l] <- mean(dupes[which(dupes$samp==i),14])
  l<-l+1
}
derepDupes <- data.frame(whichDupe, NA, d15N, NA, c13, NA, percN, percC )
#Remove dupes from orig data frame and add the derep data back

dat = dat[-grep("*[ab]", dat$Sample2),]
dat = subset(dat, dat$SampleType != "Extra")
dat2 = dat[,c(9,10, 2:5, 13, 14)]
names(derepDupes) <- names(dat2)
cndat <- rbind(dat2, derepDupes)

cndat<- cndat[which(cndat$Sample2 !="other"),]

#merge data
str(alldat$Plant)
cndat$Sample2 <- as.numeric(as.character(cndat$Sample2))

alldat2 <- merge(cndat, alldat, by.x = "Sample2", by.y = "Plant", all=T)

#Calculate NDFA
for(i in 1:length(alldat2[,1])){
  alldat2$ndfa[i] <- 100*((5.14 - alldat2$d15N....vs..air.[i])/(5.14 - min(na.omit(as.numeric(alldat2$d15N....vs..air.)))))
}

dim(cndat)
dim(alldat2)
dim(alldat)

#Sanity check
grep("294", alldat2$Plant.Number)

#End of season measurements
eos <- read.csv("./data/AsleBoxMeasurementsEndofSeason.csv")
alldat2 <- merge(eos, alldat2, by.x = "Plant.Number", by.y = "Sample2", all=T)

#Sanity check
grep("294", alldat2$Plant.Number)

#Add swainsonine results
swain <- read.csv("./data/SwainsonineMg_JHarrison.csv")
alldat2 <- merge(alldat2, swain, by.x = "Plant.Number", by.y = "sample_num", all=T)
dim(alldat2)

#add sequenced categorical
alldat2$sequenced <- "no"
sequenced <- read.csv("./data/bacteriaOTUtableZotus.csv")
seqd <- c(na.omit(as.numeric(as.character(sequenced[,1]))))

alldat2 <- alldat2[-which(duplicated(alldat2$Plant.Number)), ]
grep("294", alldat2$Plant.Number)

alldat2$sequenced[which(alldat2$Plant.Number %in% seqd)] <- "yes"

table(alldat2$sequenced)

#Add culturing information

dat <- read.csv("data/Culturing_xp_results.csv")

alldat2 <- merge(dat, alldat2, 
                by.x = "plant", 
                by.y = "Plant.Number",
                all = T)
dim(alldat2)

#Add sequencing information (these data have already been cleaned)
otus <- read.csv("./data/fungiOTUtableZotus.csv")
alldat2 <- merge(alldat2, otus, by.x = "plant", by.y = "samps", all = T)

otus <- read.csv("./data/bacteriaOTUtableZotus.csv")
names(otus) <- c("samps", paste("bact",names(otus)[2:length(otus)], sep = ""))
alldat2 <- merge(alldat2, otus, by.x = "plant", by.y = "samps", all = T)

table(alldat2$treament)
#################
#Remove those plants where treatment appeared to be unsuccessful. 
#Use this to remove swainsonine positive plants with reduced a. fulva from consideration
#these are likely those plants for which treatment was unsuccessful

#Recode the #Value! data as zero as these had very small decimal values associated with them. 
alldat2$X..Swainsonine <- as.character(alldat2$X..Swainsonine)
alldat2$X..Swainsonine[alldat2$X..Swainsonine == "#VALUE!"] <- 0
dim(alldat2)
alldat2$X..Swainsonine <- as.numeric(as.character(alldat2$X..Swainsonine))
table(alldat2$treament)

#thees had swainsonine when they were not supposed to
failed <- which(alldat2$X..Swainsonine > 0.001 & alldat2$treament %in% c("untreated_neg","treated_neg"))

alldat2 <- alldat2[-failed,]
dim(alldat2) #omitted 5 plants
table(alldat2$treament)

#Check out any plants that A. fulva was cultured from but that had their seed coats removed and didn't have swainsonine
which(alldat2$Undifilum > 1 & alldat2$treament %in% c("untreated_neg","treated_neg"))
alldat2[300,]
alldat2 <- alldat2[-300,]

# #use this to remove swainsonine negative plants that were supposed to have A. fulva
suspect <- which(alldat2$X..Swainsonine < 0.001 & alldat2$treament %in% c("treated_plus","untreated_plus"))

#retain those samples that had Undifilum or either OTU
good <- which(alldat2$Undifilum > 0 %in% c("treated_plus","untreated_plus"))
good <- c(good, which(alldat2$Otu6 > 0 %in% c("treated_plus","untreated_plus")))
good <- c(good, which(alldat2$Otu9 > 0 %in% c("treated_plus","untreated_plus")))

suspect <- suspect[!(suspect %in% good)]

#Do not consider those plants that were not sequenced or cultured 

notSurveyed <- c(which(is.na(alldat2$Undifilum)), 
                 which(is.na(alldat2$Otu9)), 
                 which(is.na(alldat2$Otu6)))

suspect <- suspect[!(suspect %in% notSurveyed)]

alldat2 <- alldat2[-suspect,]
table(alldat2$treament)

alldat2 <- alldat2[-which(alldat2$status %in% c("dead","dontuse")),]
  
#remove technical replicates, inoculum, and cultures
alldat2 <- alldat2[grep("^\\d+$",alldat2$plant),]

sum(colSums(alldat2[,grep("^Otu",names(alldat2))], na.rm = T), na.rm=T)
sum(colSums(alldat2[,grep("^bact",names(alldat2))], na.rm = T), na.rm=T)

write.csv(alldat2, file="./data/alldata2.csv", row.names = F)

