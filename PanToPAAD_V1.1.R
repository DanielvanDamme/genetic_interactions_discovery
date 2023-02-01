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
percentage = function(x, output){
  return((x[1]/nrow(uniquepatients)) * 100)
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

# Known cancer driver genes + CDKN2Ap
oncoKbGenes <- read.csv("cancerGeneList.tsv",
                        header = TRUE,
                        sep = "\t")
oncoKbGenes = data.frame(Hugo_Symbol = oncoKbGenes$Hugo.Symbol)
toadd = data.frame(Hugo_Symbol = c("CDKN2Ap14ARF", "CDKN2Ap16INK4A"))
oncoKbGenes = rbind(oncoKbGenes, toadd)

# Loading the mutations matrix
pancreasdata1 <- read.csv("data_mutations_pan1.txt",
                          header = TRUE,
                          sep = "\t")
clinicaldata <- read.csv("data_mutations_pan1.tsv",
                         header = TRUE,
                         sep = "\t")

# Get rid of hypermutators
toohighTMB <- (dplyr::filter(clinicaldata, TMB..nonsynonymous. > 20))$Sample.ID
pancreasdata1 = dplyr::filter(pancreasdata1, !(Tumor_Sample_Barcode %in% toohighTMB))
clinicaldata = dplyr::filter(clinicaldata, !(Sample.ID %in% toohighTMB))

# Only driver genes from the oncoKbGene database
pancreasdata1 = dplyr::filter(pancreasdata1, Hugo_Symbol %in% oncoKbGenes$Hugo_Symbol)

# Remove duplicate rows
pancreasdata1 = pancreasdata1[!duplicated(pancreasdata1),]

# Only the rows that differed except the first column
# filter those rows from pancreasdata1
# choose one from duplicates
duplicates <- pancreasdata1 %>% group_by_at(vars(-Hugo_Symbol)) %>% filter(n()>1)
pancreasdata1 <- anti_join(pancreasdata1, duplicates)
duplicates = duplicates[!duplicated(duplicates[, c("Start_Position")]), ]

pancreasdata1 = rbind(pancreasdata1, duplicates)

panpancandata <- dplyr::filter(clinicaldata,
                               Cancer.Type.Detailed == "Pancreatic Adenocarcinoma")

sampleIdsPA <- data.frame(Tumor_Sample_Barcode = panpancandata$Sample.ID)

temppancreasdata4 <- dplyr::filter(pancreasdata1,
                                   Tumor_Sample_Barcode %in% pull(sampleIdsPA))

pancreasdata4 <- data.frame(Hugo_Symbol = temppancreasdata4$Hugo_Symbol,
                            Tumor_Sample_Barcode = temppancreasdata4$Tumor_Sample_Barcode)

uniquepatients <- data.frame(Tumor_Sample_Barcode = unique(sampleIdsPA))
uniquegenes <- data.frame(Hugo_Symbol = unique(pancreasdata4$Hugo_Symbol))

emptymatrix <- matrix(0,
                      nrow = nrow(uniquegenes),
                      nrow(uniquepatients))

rownames(emptymatrix) <- pull(uniquegenes)
colnames(emptymatrix) <- pull(uniquepatients)

gc()
for(i in 1:nrow(pancreasdata4)){
  whichrow <- pancreasdata4[i, "Hugo_Symbol"]
  whichcol <- pancreasdata4[i, "Tumor_Sample_Barcode"]
  emptymatrix[whichrow, whichcol] = 1
}
gc()

# summing all rows for a frequency count
listofsums <- data.frame(rowSums(emptymatrix))
percentages <- data.frame(apply(listofsums, 1, percentage))

events <- discover.matrix(emptymatrix)
threshold <- 14
mode <- "less"
subset <- rowSums(emptymatrix) > threshold

result.mutex <- pairwise.discover.test(events[subset, ],
                                       alternative = c(mode))
print(result.mutex, fdr.threshold = 0.20)

signME <- as.data.frame(result.mutex, q.threshold = 0.20)

# genes <- c("IGF2BP2","REL")
# groupwise.discover.test(events[genes, ], method = "exclusivity")
# plot(events[genes, ])

signME = dplyr::filter(signME, p.value < 0.1)

# genes <- c("KRAS", "CDKN2A", "ARID1A")
# groupwise.discover.test(events[genes,], method = c("exclusivity"))
# plot(events[genes,])
#write.csv(signME, "C:\\Users\\Daan\\Desktop\\forstringdb\\mskmet.csv")
#cb(signME)

# count how many times gene is in mutual exclusive
# allgenes <- data.frame(genes = c(signME$gene1, signME$gene2))
# freqtable <- table(allgenes$genes)
# cb(freqtable)

# Methods for copying the percentages to clipboard (WINDOWS specific probably)
# #____________________________________________________________________________________________________________________________________

# data1 <- data.frame(percentage = percentages$apply.listofsums..1..percentage.)
# row_names <- data.frame(row_names = row.names(percentages))
# data1.2 <- data.frame(percentage = data1,
#                       row_names1 = row_names)
# data_1.1 <- data1.2[with(data1.2, order(-percentage)),]
# cb(data_1.1)

# # Subsetting the data to test subsampling results
# # ____________________________________________________________________________________________________________________________________

# rm(list=ls())
# gc()
# # FUNCTIONS
# percentage = function(x, output){
#   return((x[1]/nrow(uniquepatients)) * 100)
# }
# cb <- function(df, sep="\t", dec=",", max.size=(200*1000)){
#   write.table(df,
#               paste0("clipboard-",
#                      formatC(max.size,
#                              format="f",
#                              digits=0)),
#               sep=sep,
#               row.names=FALSE,
#               dec=dec)
# }
#
# # oncokbgene list to filter + CDKN2A variants
# oncoKbGenes <- read.csv("cancerGeneList.tsv",
#                         header = TRUE,
#                         sep = "\t")
# oncoKbGenes = data.frame(Hugo_Symbol = oncoKbGenes$Hugo.Symbol)
# toadd = data.frame(Hugo_Symbol = c("CDKN2Ap14ARF", "CDKN2Ap16INK4A"))
# oncoKbGenes = rbind(oncoKbGenes, toadd)
#
# # Loading the original sets
# pancreasdata1 <- read.csv("data_mutations_pan1.txt",
#                           header = TRUE,
#                           sep = "\t")
# clinicaldata <- read.csv("data_mutations_pan1.tsv",
#                          header = TRUE,
#                          sep = "\t")
#
# # Get rid of hypermutators
# toohighTMB <- (dplyr::filter(clinicaldata, TMB..nonsynonymous. > 20))$Sample.ID
# pancreasdata1 = dplyr::filter(pancreasdata1, !(Tumor_Sample_Barcode %in% toohighTMB))
# clinicaldata = dplyr::filter(clinicaldata, !(Sample.ID %in% toohighTMB))
#
# # Only driver genes from the oncoKbGene database
# pancreasdata1 = dplyr::filter(pancreasdata1, Hugo_Symbol %in% oncoKbGenes$Hugo_Symbol)
#
# # Remove duplicate rows
# pancreasdata1 = pancreasdata1[!duplicated(pancreasdata1),]
#
# # Only the rows that differed except the first column
# # filter those rows from pancreasdata1
# # choose one from duplicates
# duplicates <- pancreasdata1 %>% group_by_at(vars(-Hugo_Symbol)) %>% filter(n()>1)
# pancreasdata1 <- anti_join(pancreasdata1, duplicates)
# duplicates = duplicates[!duplicated(duplicates[, c("Start_Position")]), ]
# pancreasdata1 = rbind(pancreasdata1, duplicates)
#
# # Only PAAD patients and a subset of them
# panpancandata <- dplyr::filter(clinicaldata,
#                                Cancer.Type.Detailed == "Pancreatic Adenocarcinoma")
# panpancandata <- panpancandata[sample(nrow(panpancandata),
#                                       1700, # how many should be sampled?
#                                       replace = FALSE,
#                                       prob = NULL),]
# sampleIdsPA <- data.frame(Tumor_Sample_Barcode = panpancandata$Sample.ID)
#
# # Retrieve the data on mutations of the PAAD patients
# temppancreasdata4 <- dplyr::filter(pancreasdata1,
#                                    Tumor_Sample_Barcode %in% pull(sampleIdsPA))
# pancreasdata4 <- data.frame(Hugo_Symbol = temppancreasdata4$Hugo_Symbol,
#                             Tumor_Sample_Barcode = temppancreasdata4$Tumor_Sample_Barcode)
#
# uniquepatients <- data.frame(Tumor_Sample_Barcode = unique(sampleIdsPA))
# uniquegenes <- data.frame(Hugo_Symbol = unique(pancreasdata4$Hugo_Symbol))

# # Fill in the matrix for DISCOVER
# emptymatrix <- matrix(0,
#                       nrow = nrow(uniquegenes),
#                       nrow(uniquepatients))
# rownames(emptymatrix) <- pull(uniquegenes)
# colnames(emptymatrix) <- pull(uniquepatients)
# for(i in 1:nrow(pancreasdata4)){
#   whichrow <- pancreasdata4[i, "Hugo_Symbol"]
#   whichcol <- pancreasdata4[i, "Tumor_Sample_Barcode"]
#   emptymatrix[whichrow, whichcol] = 1
# }
#
# events <- discover.matrix(emptymatrix)
#
# counts <- data.frame(countsME = c())
#
# for(i in 2:90){
#   threshold = i
#   mode <- "less"
#   subset <- rowSums(emptymatrix) > threshold
#   result.mutex <- pairwise.discover.test(events[subset, ],
#                                          alternative = c(mode))
#   print(result.mutex, fdr.threshold = 0.20)
#   signME <- as.data.frame(result.mutex, q.threshold = 0.20)
#   signME = dplyr::filter(signME, p.value < 0.1)
#
#   toadd <- data.frame(countsME = nrow(signME))
#   counts = rbind(counts, toadd)
# }
#
# cb(counts)
# threshold <- 2
# mode <- "less"
# subset <- rowSums(emptymatrix) > threshold
# result.mutex <- pairwise.discover.test(events[subset, ],
#                                        alternative = c(mode))
# print(result.mutex, fdr.threshold = 0.20)
# signME <- as.data.frame(result.mutex, q.threshold = 0.20)
# signME = dplyr::filter(signME, p.value < 0.1)

# cb(signME)

signME = head(signME,-1)

# Bipartition for all genes in result of DISCOVER output (signME)
# __________________________________________________________________
for(i in 1:nrow(signME)){
  gene1 <- signME[i,1]
  gene2 <- signME[i,2]
  gene1vals <- c()
  gene2vals <- c()
  for(j in 1:nrow(pancreasdata4)){
    if(pancreasdata4[j,1] == gene1){
      gene1vals = append(gene1vals, pancreasdata4[j,2])
    }
    else if(pancreasdata4[j,1] == gene2){
      gene2vals = append(gene2vals, pancreasdata4[j,2])
    }
  }
  res <- data.frame(origin = c(rep(gene1, length(gene1vals)),
                               rep(gene2, length(gene2vals))),
                    Sample.ID = c(gene1vals, gene2vals))

  res = inner_join(res, clinicaldata,
                   by = "Sample.ID")

  # MSKMET categories
  # res = dplyr::select(res,
  #                     origin,
  #                     Sample.ID,
  #                     Fraction.Genome.Altered, # Continous, maybe turn into buckets
  #                     Sample.Type,
  #                     Race.Category,
  #                     Mutation.Count) # Continous, maybe turn into buckets

  # CHINA categories
  res = dplyr::select(res,
                      origin,
                      Sample.ID,
                      Tumor.Stage,
                      Sample.Type,
                      Sex,
                      Tumor.Purity,
                      Treatment)

  res = res[!duplicated(res),]

  # Those individuals that have mutations in both genes
  resduplicates <- res %>% group_by(Sample.ID) %>% filter(n()>1)
  resfortest <- anti_join(res, resduplicates)
  
  if(nrow(resduplicates)>2){
    # Ones that have mutations in both are new category
    resduplicates2 = head(resduplicates, - (nrow(resduplicates)/2))
    for(x in 1:nrow(resduplicates2)){
      resduplicates2[x,1] = "BOTH"
    }
    res = rbind(res,resduplicates2)
  }
  

  filename <- paste("C:\\Users\\Daan\\Desktop\\pantopaad\\", "row", i,"-", threshold, ".csv", sep="")

  # chi square test based on clinical data for current pair
  print(paste(gene1," ",gene2, " Row",i))
  testdata = table(res$origin,
                   res$Tumor.Stage)
  print(chisq.test(testdata, simulate.p.value = FALSE))

  #write.csv(res, filename)
}

