#J. Harrison
set.seed(666)

#input otu table, sample-well key, and taxonomy file
tax <- read.delim("data/ITSdata/combinedTaxonomy_zotus_ITS (1).txt", header=F)
#Note that it appears the second columnbn in tax is more accurate so I am using that one

otus <- read.delim("data/ITSdata/otuTableZotus_rawITS.txt", header=T)
sampNames <- read.csv("data/Sample_Well_Key_includesQuants (1).csv")
sampNames$pl_well <- paste(sampNames$plate, sampNames$well, sep=".")

#now we need to pick which OTUs to keep
#I am first removing any OTUs that the UNITE database said were plants.
#Note I also remove blank values here, spot checks showed those were likely plant

tax[,4] <- as.character(tax[,4])
fungiOnly <- tax[-which(tax[,4] %in% c("d:Plantae", "")),]

#overwrite our OTU table so that it includes only the taxon of interest
otus <- otus[which(as.character(otus$OTU_ID) %in% as.character(fungiOnly$V1)),]

#remove empty rows (no taxa of interest in a sample)
otus <- otus[which(rowSums(otus[,2:length(otus)]) != 0),]
dim(otus)

#Transpose and rename
nms <- otus$OTU_ID 
otus <- data.frame(t(otus[,2:length(otus)]))
names(otus) <- nms

row.names(otus) <- gsub("X(\\d+\\.[a-h]\\d+)_.*","\\1",row.names(otus))

#Check to see if anything in the row.names is not in the key
which(!(row.names(otus) %in% sampNames$pl_well))
#The missing sample is the control 

otus <- data.frame(otus)

#replace well based row names with the sample names
otus$samps <- row.names(otus)
newotus <- merge(sampNames,
                 otus,
                 by.y = "samps", 
                 by.x = "pl_well")

#Remove contaminants here. These were identifed via my cleaning bash script.
theBadBatch <- read.table("./data/ITSdata/OtusToRmv.txt")
newotus <- newotus[,-which(names(newotus) %in% as.character(theBadBatch$V1))]

#Reorder
newotus <- newotus[order(as.character(newotus$sampleName)),]

#Determine which contaminants are present in the control at 1% or greater
tormv <- which(newotus[which(newotus$sampleName == "neg_control"), 10:length(newotus)] / 
  colSums(na.omit(newotus[-184,10:length(newotus)])) > 0.01)

samps <- newotus$sampleName
newotus <- newotus[,10:length(newotus)]
newotus <- newotus[,-tormv]
dim(newotus)

#See retained taxa
names(newotus)
tax[which(tax[,1] %in% names(newotus)),]
#Hand blasted the ones that were just listed as fungi.
#To see what I thought. Decided to keep all of em.
dim(newotus)

write.csv(data.frame(samps, 
                     newotus),
          file="./data/fungiOTUtableZotus.csv", row.names = F)
