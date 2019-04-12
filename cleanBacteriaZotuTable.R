#J. Harrison
set.seed(666)

#input otu table, sample-well key, and taxonomy file
tax <- read.delim("./data/16s_data/combinedTaxonomy.txt", header = F)

otus <- read.delim("data/16s_data/otuTableZotus_raw (1).txt", header = T)
sampNames <- read.csv("data/Sample_Well_Key_includesQuants (1).csv")
sampNames$pl_well <- paste(sampNames$plate, sampNames$well, sep = ".")

#These are likely contaminants that I identified via my bash script.
theBadBatch <-
  read.table("./data/16s_data/possible_mt_cpDNA/allhits.b6")

#get otus in same order as taxonomy file
otus <- otus[match(tax[, 1], as.character(otus[, 1])), ]
tax <- tax[match(as.character(otus[, 1]), tax[, 1]), ]

#Sanity check
cbind(as.character(otus[, 1]),
      as.character(tax[, 1]))

#I am first removing any OTUs that are cpDNA or mitochondria.
contams <- NA

contams <-c(contams, grep("c:[Cc]hloroplast", tax[, 4]))
contams <-c(contams, grep("c:[Cc]hloroplast", tax[, 3]))

contams <-c(contams, grep("[Mm]itochondria", tax[, 4]))
contams <-c(contams, grep("[Mm]itochondria", tax[, 3]))

contams_names <- tax[contams,1]
#sanity check
tax[contams,]

length(unique(contams))

#overwrite our OTU table so that it includes only good data
#sanity check
# otus$OTU_ID
# contams_names
# which(otus$OTU_ID %in% contams_names)

otus <-  otus[-which(otus$OTU_ID %in% contams_names), ]

#Transpose and assign meaningful names
otus_t <- data.frame(t(otus[,2:length(otus)]))
names(otus_t) <- as.character(otus$OTU_ID)

#Check to see if anything in the row.names is not in the key
#The positive controls are missing, which is fine
otus_t[-which(gsub("X(\\d+\\.\\w+\\d+)\\.16S_S\\d+_L001_R1_001.fastq", "\\1", row.names(otus_t)) %in% sampNames$pl_well),
       ]

#replace well based row names with the sample names
otus_t$samps <- gsub("X(\\d+\\.\\w+\\d+)\\.16S_S\\d+_L001_R1_001.fastq", 
                     "\\1", 
                     row.names(otus_t))

newotus <- merge(sampNames,
                otus_t,
                by.y = "samps",
                by.x = "pl_well",
                all.y = T)

newotus <- newotus[order(as.character(newotus$sampleName)),]
newotus <- newotus[-1,] #Remove a pesky NA line

#Remove contaminants identified through the negative control 

#Determine which contaminants are present in the control at 1% or greater
tormv <- which((newotus[which(newotus$sampleName == "neg_control"), 10:length(newotus)] / 
  colSums(na.omit(newotus[1:188,10:length(newotus)]))) > 0.01) #Not summing the cultures and controls etc. 
samps <- newotus$sampleName
newotus <- newotus[,10:length(newotus)]
newotus <- newotus[,-tormv]
dim(newotus)

#Doublecheck that we don't have cpDNA
test <- tax[which(as.character(tax[,1]) %in% names(newotus)),]

grep("Chloroplast", test$V4)

dim(newotus)

#Calculate sum reads remaining
newotus <- data.frame(samps, 
                     newotus)
newotus <- newotus[order(as.character(newotus[,1])),]
test <- newotus[-grep("[A-Za-z]",as.character(newotus[,1])),]
test$samps <- as.character(test$samps)

sum(colSums(test[1:172,2:length(test)])) #doesnt include positive controls

write.csv(#test,
          newotus,
          file="./data/bacteriaOTUtableZotusALLSEQS.csv", 
          row.names = F)
