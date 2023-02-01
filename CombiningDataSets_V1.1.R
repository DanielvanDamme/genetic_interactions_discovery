# Libraries
# Novel method to determine candidate genetic interactions
library(discover)
# Permutation of tables, vectors arrays and matrices
library(dplyr)
# Chi-square tests
library(MASS)

rm(list=ls())
gc()

# FUNCTIONS
# Simple convert number of patients and amount of genes   
# to obtain alteriation percentage of that gene
percentage <- function(x, output){
  return((x[1]/nrow(uniquepatients))*100)
}
cb <- function(df, sep="\t", dec=",", max.size=(200*1000)){
  write.table(df, 
              paste0("clipboard-", 
                     formatC(max.size, 
                             format="f", 
                             digits=0)), 
              sep=sep, 
              row.names=FALSE, 
              dec=dec)
}

# Known cancer driver genes
oncoKbGenes <- read.csv("cancerGeneList.tsv",
                        header = TRUE,
                        sep = "\t")
# Preparing multiple combined datasets
combinedsets = data.frame(Hugo_Symbol = c(),
                          Tumor_Sample_Barcode = c())
combinedclinical = data.frame(patientid=c(),
                              Sample.ID=c())


# Standard filename
filename = "data_mutations_"

# Loading multiple datasets taking everything interesting
for(i in 2:2){
  file <- paste(filename,i,".txt",sep="")
  clinicalfile <- paste(filename,i,".tsv",sep="")
  
  temp <- read.csv(file,
                   header = TRUE,
                   sep = "\t")
  clinicaltemp <- read.csv(clinicalfile,
                           header = TRUE,
                           sep = "\t")
  
  newrows <- data.frame(Hugo_Symbol = temp$Hugo_Symbol, 
                        Tumor_Sample_Barcode = temp$Tumor_Sample_Barcode)
  newclinical <- data.frame(patientid = clinicaltemp$Patient.ID,
                            Sample.ID = clinicaltemp$Sample.ID,
                            tmb = clinicaltemp$TMB..nonsynonymous.)
  
  
  combinedsets = rbind(combinedsets, newrows)
  combinedclinical = rbind(combinedclinical, newclinical)
  
  rm(newrows, temp, file, newclinical, clinicaltemp, clinicalfile)
}

# Get rid of hypermutators
toohighTMB <- (dplyr::filter(combinedclinical, tmb > 20))$Sample.ID
combinedsets = dplyr::filter(combinedsets, !(Tumor_Sample_Barcode %in% toohighTMB))
combinedclinical = dplyr::filter(combinedclinical, !(Sample.ID %in% toohighTMB))

uniquepatients <- data.frame(Tumor_Sample_Barcode = unique(combinedsets$Tumor_Sample_Barcode))
uniquegenes <- data.frame(Hugo_Symbol = unique(combinedsets$Hugo_Symbol))

uniquegenes = dplyr::filter(uniquegenes, Hugo_Symbol %in% oncoKbGenes$Hugo.Symbol)

emptymatrix <- matrix(0,
                      nrow = nrow(uniquegenes),
                      nrow(uniquepatients))

rownames(emptymatrix) <- pull(uniquegenes)
colnames(emptymatrix) <- pull(uniquepatients)

combinedsets = dplyr::filter(combinedsets, Hugo_Symbol %in% oncoKbGenes$Hugo.Symbol)

for(i in 1:nrow(combinedsets)){
  whichrow <- combinedsets[i, "Hugo_Symbol"]
  whichcol <- combinedsets[i, "Tumor_Sample_Barcode"]
  emptymatrix[whichrow, whichcol] = 1
}
rm(uniquegenes)
gc()

# summing all rows for a frequency count
listofsums <- data.frame(rowSums(emptymatrix))
percentages <- data.frame(apply(listofsums, 1, percentage))

# Stratification of the datasets
list1 <- c(rep(1,109))
list2 <- c(rep(2,175))
list3 <- c(rep(3,381))
list4 <- c(rep(4,140))

# stratification of which experiment
whichstudy <- data.frame(experiment = uniquepatients,
                         strata = c(rep(0,nrow(uniquepatients))))
whichstudy[1] <- lapply(whichstudy[1], substr, 1, 3)
for (i in 1:nrow(uniquepatients)) {
  if(whichstudy[i,1] == "PDA"){
    whichstudy[i,2] = 1
  }  
  else if(whichstudy[i,1] == "TCG"){
    whichstudy[i,2] = 2
  }
  else if(whichstudy[i,1] == "ICG"){
    whichstudy[i,2] = 3
  }
  else if(whichstudy[i,1] == "GAR"){
    whichstudy[i,2] = 4
  }
  else if(whichstudy[i,1] == "C3N"){
    whichstudy[i,2] = 5
  }
  else if(whichstudy[i,1] == "C3L"){
    whichstudy[i,2] = 6
  }
}
strataexperiment <- whichstudy$strata

# estimation of the background matrix
stratatest <- array(c(list1, list2, list3, list4))
events <- discover.matrix(emptymatrix)
threshold <- 10
subset <- rowSums(emptymatrix) > threshold

result.mutex <- pairwise.discover.test(events[subset, ],
                                       alternative = c("greater"))
print(result.mutex, fdr.threshold = 0.2)

signME <- as.data.frame(result.mutex, q.threshold = 0.2)

# Groupwise comparisons
# genes <- c("TGFBR2","TP53","GNAS")
# groupwise.discover.test(events[genes, ],
#                         method = "exclusivity")
# plot(events[genes, ])

signME = dplyr::filter(signME, p.value < 0.1)

#write.csv(signME, "C:\\Users\\Daan\\Desktop\\forstringdb\\combined.csv")

cb(signME)

# allgenes <- data.frame(genes = c(signME$gene1, signME$gene2))
# freqtable <- table(allgenes$genes)
# cb(freqtable)

# Copy the percentages per gene
# data1 <- data.frame(percentage = percentages$apply.listofsums..1..percentage.)
# row_names <- data.frame(row_names = row.names(percentages))
# data1.2 <- data.frame(percentage = data1,
#                       row_names1 = row_names)
# data_1.1 <- data1.2[with(data1.2, order(-percentage)),]
# cb(data_1.1)

# ______________________________________________________________________________

# Get other clinicaldata to attach to gene pairs (bipartition)
# clinical1 <- read.csv("data_mutations_1.tsv",
#                       header = TRUE,
#                       sep = "\t")
# clinical2 <- read.csv("data_mutations_2.tsv",
#                       header = TRUE,
#                       sep = "\t")
# clinical3 <- read.csv("data_mutations_3.tsv",
#                       header = TRUE,
#                       sep = "\t")
# clinical4 <- read.csv("data_mutations_4.tsv",
#                       header = TRUE,
#                       sep = "\t")
# 
# 
# toadd1 <- dplyr::select(clinical1, #dataset
#                  Sample.ID = Sample.ID,
#                  race = Race.Category,
#                  sex = Sex)
# toadd2 <- dplyr::select(clinical2, #dataset
#                  Sample.ID = Sample.ID,
#                  race = Race.Category,
#                  cancerstate = Person.Neoplasm.Cancer.Status,
#                  sex = Sex)
# toadd3 <- dplyr::select(clinical3, #dataset
#                  Sample.ID = Sample.ID,
#                  country = Country,
#                  race = Ethnicity.Category,
#                  sex = Sex)
# toadd4 <- dplyr::select(clinical4, #dataset
#                  Sample.ID = Sample.ID,
#                  country = Participant.Country,
#                  race = Race,
#                  sex = Sex)
# 
# # sets that have country tag
# #countrycompare <- rbind(toadd3, toadd4)
# # sets that have sex and race (everything)
# compareadd <- bind_rows(toadd1, toadd2, toadd3, toadd4)
# 
# rm(clinical1, clinical2, clinical3,clinical4, toadd1, toadd2, toadd3, toadd4)
# 
# # Bipartition for all genes in result of DISCOVER output (signME)
# for(i in 1:nrow(signME)){
#   gene1 <- signME[i,1]
#   gene2 <- signME[i,2]
#   gene1vals <- c()
#   gene2vals <- c()
# 
#   # Every tumor barcode that has mutation in gene add to list
#   for(j in 1:nrow(combinedsets)){
#     if(combinedsets[j,1] == gene1){
#       gene1vals = append(gene1vals, combinedsets[j,2])
#     }
#     else if(combinedsets[j,1] == gene2){
#       gene2vals = append(gene2vals, combinedsets[j,2])
#     }
#   }
#   # list of pairs of genes with their tumor barcode
#   res <- data.frame(origin = c(rep(gene1, length(gene1vals)),
#                                rep(gene2, length(gene2vals))),
#                     Sample.ID = c(gene1vals, gene2vals))
#   #Save location
#   filename <- paste("C:\\Users\\Daan\\Desktop\\resultscombined\\", "row", i,"-", threshold, ".csv", sep="")
# 
#   # Join clinical data
#   res = inner_join(res, compareadd,
#                    by = "Sample.ID")
#   res = res[!duplicated(res),]
#   
#   # Remove from the group those that are mutated in both
#   resduplicates <- res %>% group_by(Sample.ID) %>% filter(n()>1)
#   resfortest <- anti_join(res, resduplicates)
#   # Add which experiment they come from
#   resfortest[7] <- lapply(resfortest[2], substr, 1, 3)
#   
#   # Add which experiment they come from
#   res[7] <- lapply(res[2], substr, 1, 3)
#   
#   # chi square test based on clinical data for current pair
#   print(paste("Row",i))
#   
#   # !!!! OPTIONS !!!!
#   # Options are country, sex, race, Sample.ID.1 (which experiment)
#   testdata = table(resfortest$origin,
#                    resfortest$race)
#   print(chisq.test(testdata, simulate.p.value = TRUE))
#   
#   # Write every row into file for use with Cytoscape for example
#   #write.csv(res, filename)
# }
