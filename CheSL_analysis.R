#----------------------------------------------------------- Packages -----------------------------------------------------------------
#Some packages did not end up being used in the final analysis, but these should cover all depdencies

library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
source("scripts/ancom_v2.1.R")
library(readr)
library(tidyverse)
library(igraph)
library(visNetwork)
library('tidyverse')
library('vegan')
library('reshape2')
library('seqtime')
library('phyloseq')
library('dplyr')
library('evobiR')
library('slider')
library('runner')
library('ggVennDiagram')
library('ggallin')
library('limma')
library(phyloseqCompanion)
library(microViz)
#---------------------------------------------------------- Data Import ---------------------------------------------------------------

#Read in Count table and Taxonomy file, use full join to add taxonomy to the count table
chesl_asvs <- read.csv(file = "chesl_asvs.csv")
tax_table <- chesl_asvs[,c(694:700,694)]
otu_table <- chesl_asvs[,c(1:693)]
metadata <- read.csv(file = "metacomb.csv")

MCC_time <- read.csv("EarlyMidLate.csv")
rich_time <- read.csv(file = "EarlyMidLate.csv")
sub_list<- read.csv(file="exp1_carbohydrates.csv", header = TRUE)

media_summary <- read.csv(file="diversityComparisons.csv", header = TRUE)

chung_chesl <- read.csv("Chung_CheSL.csv")
Lindemann_chesl <- read.csv("Lindemann_CheSL.csv")
Loss_simp <- read.csv(file = "Loss_Inv_Simpson.csv")

node_table <- read.csv(file = "node_list.csv", header = TRUE, row.names = FALSE)
AP_network <- png::readPNG("AR_Genus_best.csv.png")
AM_network <- png::readPNG("ARM_Genus_best.csv.png")
B_network <- png::readPNG("ES_Genus_best.csv.png")
GM_network <- png::readPNG("HG_Genus_best.csv.png")
ILC_network <- png::readPNG("BRM4IS_Genus_best.csv.png")

# ========================================================== CheSL Calculation Function =================================================

#This version does not currently aggregate Subunit-Linkage conformations that are shared across carbohydrates 
#If there are overlaps in Subunits and linkages, it is recommended to perform these calculations stepwise
#Then aggregate fractions prior to CheSL diversity measure calculations

CheSL <- function(MediaMatrix, SubunitLinkageMatrix) {
  MediaMatrix <- t(MediaMatrix)
  SubunitLinkageMatrix <- t(SubunitLinkageMatrix)
  if (ncol(MediaMatrix) != nrow(SubunitLinkageMatrix)) {
    stop("The number of columns in the Media matrix must match the number of rows in SubunitLinkage matrix for multiplication.")
  }
  
  product_matrix <- MediaMatrix %*% SubunitLinkageMatrix
  product_matrix <- t(product_matrix)
  
  normalized_matrix <- apply(product_matrix, 2, function(col) col / sum(col))
  CheSL_Richness <- apply(normalized_matrix, 2, function(col) sum(col != 0))
  

  CheSL_Shannon <- apply(normalized_matrix, 2, function(col) {
    p <- col[col > 0] 
    -sum(p * log(p))
  })
  
  CheSL_Evenness <- CheSL_Shannon / log(CheSL_Richness)

  return(list(
    CheSL_Richness = CheSL_Richness,
    CheSL_Shannon = CheSL_Shannon,
    CheSL_Evenness = CheSL_Evenness
  ))
}

# ---------------------------------------------------------- Other Functions ----------------------------------------------------------
limitsOutput <- function(data) {
  data <- t(data)
  limitsData <- limits(data)$Aest
  return(limitsData)
}

ColVar <- function(x, ...) {
  vegdist(x, method="bray")
}

limitsAnalysis <- function(data){
  for (i in windowSize){
    #get the number of windows
    nWindowLength <- 11 - i
    nWindow <- 1:nWindowLength
    #get name of input variable to use for saving files
    dataName <- deparse(substitute(data))
    data_loop <- data
    
    windowList <- vector(mode="list", length = j)
    for (j in nWindow){
      #create a vector of the given time points for window j
      jfin <- j
      print(jfin)
      print(i)
      print("initial data")
      print(sample_data(data))
      timeSeries <- sample_data(data)$Time
      jwindow <- timeSeries[(jfin):(jfin+i-1)]
      print("Jwindow")
      print(jwindow)
      print(sample_data(data)$Time)
      data_loop <- ps_filter(data, Time %in% jwindow)
      print("initial data loop")
      print(sample_data(data_loop))
      
      
      otu_table(data_loop) <- normalize(otu_table(data_loop))
      nonzero_data <- genefilter_sample(data_loop, filterfun_sample(function(x) x > 1e-7), A=0.95*nsamples(data_loop))
      data_loop <- prune_taxa(nonzero_data, data_loop)
      
      #normalize the dataset
      #data_otu <- otu_table(data_loop)
      #data_otu <- normalize(data_otu)
      #data_table_otu <- data_otu
      #print otu tables
      #print(data_table_otu)
      
      #If filtering in first place
      data_otu <- as.data.frame(otu_table(data_loop))
      data_table_otu <- data_otu
      
      #determine mean cross correlation
      data_otu_limits <- tryCatch(limits(data_table_otu, bagging.iter=50), error = function(e) { cat('NA'); print(e); e }) #runs to here, need to stop it "stopping" when error occurs
      data_limits <- data_otu_limits$Aest
      tryCatch(rownames(data_limits)<- rownames(data_table_otu), error = function(e) { cat('NA'); print(e); e }) 
      tryCatch(colnames(data_limits)<-rownames(data_table_otu), error = function(e) { cat('NA'); print(e); e })
      #print(data_limits)
      
      #take timeseries after start
      kwindow <- timeSeries[(jfin):length(timeSeries)]
      print("Kwindow")
      print(kwindow)
      print("total frame")
      print(sample_data(data)$Time)
      data_loop_check <- ps_filter(data, Time %in% kwindow)
      my_subset <- subset(otu_table(data_loop_check), rownames(otu_table(data_loop_check)) %in% row.names(otu_table(data_loop)))
      data_loop_check <- merge_phyloseq(my_subset, tax_table(data_loop), sample_data(data_loop))
      data_otu_check <- as.data.frame(otu_table(data_loop_check))
      data_table_otu_check <- cbind(row.names(data_otu_check), data_otu_check)
      colnames(data_table_otu_check)[1]<- "OTU"
      genera <- tax[,c(1,7)]
      genera$Genus <- make.unique(genera$Genus, sep="-")
      colnames(genera)[1] <- "OTU"
      data_table_otu_Genus_check <- merge(x=data_table_otu_check, y=genera, by.x="OTU", by.y="OTU", all.x=TRUE)
      rownames(data_table_otu_Genus_check) <- data_table_otu_Genus_check$Genus
      data_table_otu_Genus_check <- data_table_otu_Genus_check[,-grep("OTU",colnames(data_table_otu_Genus_check))]
      data_table_otu_Genus_check <- data_table_otu_Genus_check[,-grep("Genus",colnames(data_table_otu_Genus_check))]
      data_table_otu_check <- data_table_otu_Genus_check
      print("data_table_otu_check")
      print(data_table_otu_check)
      #data_otu_slidequal1 <- tryCatch(limitsQuality(data_table_otu, A=data_limits, predict.stepwise = FALSE, plot=TRUE), error = function(e) { cat('NA'); print(e); e })
      data_otu_slidequal1 <- tryCatch(limitsQuality(data_table_otu, A=data_limits, plot=TRUE), error = function(e) { cat('NA'); print(e); e })
      
      #take final value of mcc
      data_mcc <- data_otu_slidequal1$meancrosscor[length(data_otu_slidequal1$meancrosscor)]
      cat("MCC VALUE:", data_mcc)
      #make list of the mean cross correlation and the interaction matrix
      data_list <- list(data_mcc, data_limits)
      
      #assign this list to a variable denoting which time window is represented
      
      windowList[[j]] <- data_list
      #print(windowList)
    }
    print(windowList)
    windowAnalysis[[i]] <-  windowList
    #print(windowAnalysis)
  }
  #return list of list of list of matrices and mcc
  return(windowAnalysis)
}

# ---------------------------------------------------------- Analysis ----------------------------------------------------------------

otu_table <- counts[,3:304]
row.names(otu_table) <- counts[,1]
tax_table <- tax[,1:7]
row.names(tax_table) <- tax[,1]
tax_table <- as.matrix(tax_table)
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

colnames(tax_table) <- c("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

rownames(otu_table) <- otu_table[,1]
otu_table <- otu_table[,-1]

rownames(tax_table) <- tax_table[,1]
tax_table <- tax_table[,-1]

OTU = otu_table(as.matrix(otu_table), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_table))
physeq = phyloseq(OTU,TAX)
meta = sample_data(metadata)
physeq1 = merge_phyloseq(physeq,meta)

OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
physeq = phyloseq(OTU,TAX)
meta = sample_data(br35meta)
physeq1 = merge_phyloseq(physeq,meta)
nonzero = prune_taxa(taxa_sums(physeq1)>5,physeq1)

windowSize <- 6:10

ES_1 <- subset_samples(ES, Replicate =="1")
ES_2 <-subset_samples(ES, Replicate == "2")
ES_3 <-subset_samples(ES, Replicate == "3")
ES_4 <-subset_samples(ES, Replicate == "4")

AR_1 <- subset_samples(AR, Replicate == "1")
AR_2 <- subset_samples(AR, Replicate == "2")
AR_3 <- subset_samples(AR, Replicate == "3")
AR_4 <- subset_samples(AR, Replicate == "4")

ARM_1 <- subset_samples(ARM, Replicate == "1")
ARM_2 <- subset_samples(ARM, Replicate == "2")
ARM_3 <- subset_samples(ARM, Replicate == "3")
ARM_4 <- subset_samples(ARM, Replicate == "4")

HS_1 <- subset_samples(HS, Replicate == "1")
HS_2 <- subset_samples(HS, Replicate == "2")
HS_3 <- subset_samples(HS, Replicate == "3")
HS_4 <- subset_samples(HS, Replicate == "4")

HG_1 <- subset_samples(HG, Replicate == "1")
HG_2 <- subset_samples(HG, Replicate == "2")
HG_3 <- subset_samples(HG, Replicate == "3")
HG_4 <- subset_samples(HG, Replicate == "4")

BRM4IS_1 <- subset_samples(BRM4IS, Replicate == "1")
BRM4IS_2 <- subset_samples(BRM4IS, Replicate == "2")
BRM4IS_3 <- subset_samples(BRM4IS, Replicate == "3")
BRM4IS_4 <- subset_samples(BRM4IS, Replicate == "4")

ES_1_limitsAnalysis <- limitsAnalysis(ES_1)
ES_2_limitsAnalysis <- limitsAnalysis(ES_2)
ES_3_limitsAnalysis <- limitsAnalysis(ES_3)
ES_4_limitsAnalysis <- limitsAnalysis(ES_4)

AR_1_limitsAnalysis <- limitsAnalysis(AR_1)
AR_2_limitsAnalysis <- limitsAnalysis(AR_2)
AR_3_limitsAnalysis <- limitsAnalysis(AR_3)
AR_4_limitsAnalysis <- limitsAnalysis(AR_4)

ARM_1_limitsAnalysis <- limitsAnalysis(ARM_1)
ARM_2_limitsAnalysis <- limitsAnalysis(ARM_2)
ARM_3_limitsAnalysis <- limitsAnalysis(ARM_3)
ARM_4_limitsAnalysis <- limitsAnalysis(ARM_4)

HS_1_limitsAnalysis <- limitsAnalysis(HS_1)
HS_2_limitsAnalysis <- limitsAnalysis(HS_2)
HS_3_limitsAnalysis <- limitsAnalysis(HS_3)
HS_4_limitsAnalysis <- limitsAnalysis(HS_4)

HG_1_limitsAnalysis <- limitsAnalysis(HG_1)
HG_2_limitsAnalysis <- limitsAnalysis(HG_2)
HG_3_limitsAnalysis <- limitsAnalysis(HG_3)
HG_4_limitsAnalysis <- limitsAnalysis(HG_4)

BRM4IS_1_limitsAnalysis <- limitsAnalysis(BRM4IS_1)
BRM4IS_2_limitsAnalysis <- limitsAnalysis(BRM4IS_2)
BRM4IS_3_limitsAnalysis <- limitsAnalysis(BRM4IS_3)
BRM4IS_4_limitsAnalysis <- limitsAnalysis(BRM4IS_4)

FS17 <- subset_samples(nonzeroEven1, FS == "FS17" & Time <= 240 & Time != 192)
FS18 <- subset_samples(nonzeroEven1, FS == "FS18" & Time <= 240 & Time != 192)

FS17_BRM3 <- subset_samples(FS17, Media == "BRM3")
FS17_BRM4 <- subset_samples(FS17, Media == "BRM4")
FS17_BRM4I <- subset_samples(FS17, Media == "BRM4I")

FS18_BRM3 <- subset_samples(FS18, Media == "BRM3")
FS18_BRM4 <- subset_samples(FS18, Media == "BRM4")
FS18_BRM4I <- subset_samples(FS18, Media == "BRM4I")
FS18_BRM5 <- subset_samples(FS18, Media == "BRM5")

FS17_BRM3A <- subset_samples(FS17_BRM3, Replicate =="A")
FS17_BRM3B <-subset_samples(FS17_BRM3, Replicate == "B")

FS17_BRM3A_limitsAnalysis <- limitsAnalysis(FS17_BRM3A)
FS17_BRM3B_limitsAnalysis <- limitsAnalysis(FS17_BRM3B)

FS17_BRM4A <- subset_samples(FS17_BRM4, Replicate =="A")
FS17_BRM4B <-subset_samples(FS17_BRM4, Replicate == "B")
FS17_BRM4C <-subset_samples(FS17_BRM4, Replicate == "C")

FS17_BRM4A_limitsAnalysis <- limitsAnalysis(FS17_BRM4A)
FS17_BRM4B_limitsAnalysis <- limitsAnalysis(FS17_BRM4B)
FS17_BRM4C_limitsAnalysis <- limitsAnalysis(FS17_BRM4C)

FS17_BRM4IA <- subset_samples(FS17_BRM4I, Replicate =="A")
FS17_BRM4IB <-subset_samples(FS17_BRM4I, Replicate == "B")
FS17_BRM4IC <-subset_samples(FS17_BRM4I, Replicate == "C")

FS17_BRM4IA_limitsAnalysis <- limitsAnalysis(FS17_BRM4IA)
FS17_BRM4IB_limitsAnalysis <- limitsAnalysis(FS17_BRM4IB)
FS17_BRM4IC_limitsAnalysis <- limitsAnalysis(FS17_BRM4IC)

FS18_BRM3A <- subset_samples(FS18_BRM3, Replicate =="A")
FS18_BRM3B <-subset_samples(FS18_BRM3, Replicate == "B")
FS18_BRM3C <-subset_samples(FS18_BRM3, Replicate == "C")

FS18_BRM3A_limitsAnalysis <- limitsAnalysis(FS18_BRM3A)
FS18_BRM3B_limitsAnalysis <- limitsAnalysis(FS18_BRM3B)
FS18_BRM3C_limitsAnalysis <- limitsAnalysis(FS18_BRM3C)

FS18_BRM4A <- subset_samples(FS18_BRM4, Replicate =="A")
FS18_BRM4B <-subset_samples(FS18_BRM4, Replicate == "B")
FS18_BRM4C <-subset_samples(FS18_BRM4, Replicate == "C")

FS18_BRM4A_limitsAnalysis <- limitsAnalysis(FS18_BRM4A)
FS18_BRM4B_limitsAnalysis <- limitsAnalysis(FS18_BRM4B)
FS18_BRM4C_limitsAnalysis <- limitsAnalysis(FS18_BRM4C)

FS18_BRM4IA <- subset_samples(FS18_BRM4I, Replicate =="A")
FS18_BRM4IB <-subset_samples(FS18_BRM4I, Replicate == "B")
FS18_BRM4IC <-subset_samples(FS18_BRM4I, Replicate == "C")

FS18_BRM4IA_limitsAnalysis <- limitsAnalysis(FS18_BRM4IA)
FS18_BRM4IB_limitsAnalysis <- limitsAnalysis(FS18_BRM4IB)
FS18_BRM4IC_limitsAnalysis <- limitsAnalysis(FS18_BRM4IC)

FS18_BRM5A <- subset_samples(FS18_BRM5, Replicate =="A")
FS18_BRM5B <-subset_samples(FS18_BRM5, Replicate == "B")
FS18_BRM5C <-subset_samples(FS18_BRM5, Replicate == "C")

FS18_BRM5A_limitsAnalysis <- limitsAnalysis(FS18_BRM5A)
FS18_BRM5B_limitsAnalysis <- limitsAnalysis(FS18_BRM5B)
FS18_BRM5C_limitsAnalysis <- limitsAnalysis(FS18_BRM5C)

fs17_brm3A_best <- FS17_BRM3A_limitsAnalysis[[6]][[2]][[2]]
fs17_brm3B_best <- FS17_BRM3B_limitsAnalysis[[6]][[2]][[2]]

fs18_brm3A_best <- FS18_BRM3A_limitsAnalysis[[6]][[2]][[2]]
fs18_brm3B_best <- FS18_BRM3B_limitsAnalysis[[6]][[2]][[2]]
fs18_brm3C_best <- FS18_BRM3C_limitsAnalysis[[6]][[3]][[2]]

fs17_brm4A_best <- FS17_BRM4A_limitsAnalysis[[6]][[2]][[2]]
fs17_brm4B_best <- FS17_BRM4B_limitsAnalysis[[6]][[3]][[2]]
fs17_brm4C_best <- FS17_BRM4C_limitsAnalysis[[6]][[2]][[2]]

fs18_brm4A_best <- FS18_BRM4A_limitsAnalysis[[6]][[2]][[2]]
fs18_brm4B_best <- FS18_BRM4B_limitsAnalysis[[6]][[2]][[2]]
fs18_brm4C_best <- FS18_BRM4C_limitsAnalysis[[6]][[2]][[2]]

fs17_brm4IA_best <- FS17_BRM4IA_limitsAnalysis[[6]][[2]][[2]]
fs17_brm4IB_best <- FS17_BRM4IB_limitsAnalysis[[6]][[2]][[2]]
fs17_brm4IC_best <- FS17_BRM4IC_limitsAnalysis[[6]][[2]][[2]]

fs18_brm4IA_best <- FS18_BRM4IA_limitsAnalysis[[6]][[2]][[2]]
fs18_brm4IB_best <- FS18_BRM4IB_limitsAnalysis[[6]][[3]][[2]]
fs18_brm4IC_best <- FS18_BRM4IC_limitsAnalysis[[6]][[3]][[2]]

fs18_brm5A_best <- FS18_BRM5A_limitsAnalysis[[6]][[2]][[2]]
fs18_brm5B_best <- FS18_BRM5B_limitsAnalysis[[6]][[2]][[2]]
fs18_brm5C_best <- FS18_BRM5C_limitsAnalysis[[6]][[3]][[2]]

AR_1_best <- AR_1_limitsAnalysis[[6]][[9]][[2]]
AR_2_best <- AR_2_limitsAnalysis[[6]][[8]][[2]]
AR_3_best <- AR_3_limitsAnalysis[[6]][[10]][[2]]
AR_4_best <- AR_4_limitsAnalysis[[6]][[9]][[2]]

ARM_1_best <- ARM_1_limitsAnalysis[[6]][[7]][[2]]
ARM_2_best <- ARM_2_limitsAnalysis[[6]][[1]][[2]]
ARM_3_best <- ARM_3_limitsAnalysis[[6]][[8]][[2]]
ARM_4_best <- ARM_4_limitsAnalysis[[6]][[3]][[2]]

ES_1_best <- ES_1_limitsAnalysis[[6]][[8]][[2]]
ES_2_best <- ES_2_limitsAnalysis[[6]][[9]][[2]]
ES_3_best <- ES_3_limitsAnalysis[[6]][[9]][[2]]
ES_4_best <- ES_4_limitsAnalysis[[6]][[10]][[2]]

HG_1_best <- HG_1_limitsAnalysis[[6]][[6]][[2]]
HG_2_best <- HG_2_limitsAnalysis[[6]][[8]][[2]]
HG_3_best <- HG_3_limitsAnalysis[[6]][[9]][[2]]
HG_4_best <- HG_4_limitsAnalysis[[6]][[2]][[2]]

HS_1_best <- HS_1_limitsAnalysis[[6]][[1]][[2]]
HS_2_best <- HS_2_limitsAnalysis[[6]][[7]][[2]]
HS_3_best <- HS_3_limitsAnalysis[[6]][[1]][[2]]
HS_4_best <- HS_4_limitsAnalysis[[6]][[1]][[2]]

BRM4IS_1_best <- BRM4IS_1_limitsAnalysis[[6]][[8]][[2]]
BRM4IS_2_best <- BRM4IS_2_limitsAnalysis[[6]][[8]][[2]]
BRM4IS_3_best <- BRM4IS_3_limitsAnalysis[[6]][[9]][[2]]
BRM4IS_4_best <- BRM4IS_4_limitsAnalysis[[6]][[9]][[2]]

fs17_brm3A_best_melt <- melt(fs17_brm3A_best)
fs17_brm3B_best_melt <- melt(fs17_brm3B_best)
fs17_brm3AB_best_melt <- merge(fs17_brm3A_best_melt, fs17_brm3B_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs17_brm3AB_best_melt[is.na(fs17_brm3AB_best_melt)] <- 0


fs17_brm3AB_best_genus <- merge(fs17_brm3AB_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
fs17_brm3AB_best_genus$Var1 <- fs17_brm3AB_best_genus$Genus
colnames(fs17_brm3AB_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Genus1")

fs17_brm3AB_best_genus <- merge(fs17_brm3AB_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(fs17_brm3AB_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Genus1", "Genus2")
fs17_brm3AB_best_genus$Var2 <- fs17_brm3AB_best_genus$Genus
fs17_brm3AB_best_genus <- fs17_brm3AB_best_genus[,-c(1,2)]


fs17_brm3_Genus_best <- aggregate(cbind(Rep1, Rep2)~Genus1 + Genus2,fs17_brm3AB_best_genus,sum)
fs17_brm3_Genus_best <- cbind(fs17_brm3_Genus_best, rowSums(fs17_brm3_Genus_best[,3:4]))
colnames(fs17_brm3_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Sum")
fs17_brm3_Genus_best <- fs17_brm3_Genus_best[fs17_brm3_Genus_best$Sum != 0,]
fs17_brm3_Genus_bestNet <- fs17_brm3_Genus_best[,3:4]
fs17_brm3_Genus_bestNet[fs17_brm3_Genus_bestNet < 0] <- -1
fs17_brm3_Genus_bestNet[fs17_brm3_Genus_bestNet > 0] <- 1
fs17_brm3_Genus_best <- cbind(fs17_brm3_Genus_best, rowSums(fs17_brm3_Genus_bestNet)) 
fs17_brm3_Genus_best <- fs17_brm3_Genus_best[abs(fs17_brm3_Genus_best$`rowSums(fs17_brm3_Genus_bestNet)`)> 1,]
write.csv(fs17_brm3_Genus_best, file="fs17_brm3_Genus_best.csv")

fs17_brm3_best <- merge(fs17_brm3AB_best_melt, genera, by.x = "Var1", by.y = "OTU", all=TRUE)
fs17_brm3_best$Var1 <- fs17_brm3_best$Genus
colnames(fs17_brm3_best) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Genus1")

fs17_brm3_best <- merge(fs17_brm3_best, genera, by.x = "OTU2", by.y = "OTU", all=TRUE)
fs17_brm3_best$OTU2 <- fs17_brm3_best$Genus

fs17_brm3_best <- fs17_brm3_best[,-c("Genus1", "Genus")]
fs17_brm3_best <- drop_na(fs17_brm3_best)
fs17_brm3_best$Int <- rowSums(fs17_brm3_best[,c("Rep1", "Rep2")])
fs17_brm3_best$Rep1V <- fs17_brm3_best$Rep1
fs17_brm3_best$Rep2V <- fs17_brm3_best$Rep2
fs17_brm3_best$Rep1V[fs17_brm3_best$Rep1V >0] <- 1
fs17_brm3_best$Rep1V[fs17_brm3_best$Rep1V <0] <- -1
fs17_brm3_best$Rep2V[fs17_brm3_best$Rep2V >0] <- 1
fs17_brm3_best$Rep2V[fs17_brm3_best$Rep2V <0] <- -1
fs17_brm3_best$Con <- rowSums(fs17_brm3_best[,c("Rep1V","Rep2V")])
fs17_brm3_best <- fs17_brm3_best[fs17_brm3_best$Con != 0,]
write.csv(fs17_brm3_best, file ="fs17_brm3_best.csv")

fs17_brm4A_best_melt <- melt(fs17_brm4A_best)
fs17_brm4B_best_melt <- melt(fs17_brm4B_best)
fs17_brm4C_best_melt <- melt(fs17_brm4C_best)
fs17_brm4AB_best_melt <- merge(fs17_brm4A_best_melt, fs17_brm4B_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs17_brm4ABC_best_melt <- merge(fs17_brm4AB_best_melt, fs17_brm4C_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs17_brm4ABC_best_melt[is.na(fs17_brm4ABC_best_melt)] <- 0


fs17_brm4ABC_best_genus <- merge(fs17_brm4ABC_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
fs17_brm4ABC_best_genus$Var1 <- fs17_brm4ABC_best_genus$Genus
colnames(fs17_brm4ABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1")

fs17_brm4ABC_best_genus <- merge(fs17_brm4ABC_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(fs17_brm4ABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1", "Genus2")
fs17_brm4ABC_best_genus$Var2 <- fs17_brm4ABC_best_genus$Genus
fs17_brm4ABC_best_genus <- fs17_brm4ABC_best_genus[,-c(1,2)]


fs17_brm4_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3)~Genus1 + Genus2,fs17_brm4ABC_best_genus,sum)
fs17_brm4_Genus_best <- cbind(fs17_brm4_Genus_best, rowSums(fs17_brm4_Genus_best[,3:5]))
colnames(fs17_brm4_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Sum")
fs17_brm4_Genus_best <- fs17_brm4_Genus_best[fs17_brm4_Genus_best$Sum != 0,]
fs17_brm4_Genus_bestNet <- fs17_brm4_Genus_best[,3:5]
fs17_brm4_Genus_bestNet[fs17_brm4_Genus_bestNet < 0] <- -1
fs17_brm4_Genus_bestNet[fs17_brm4_Genus_bestNet > 0] <- 1
fs17_brm4_Genus_best <- cbind(fs17_brm4_Genus_best, rowSums(fs17_brm4_Genus_bestNet))
fs17_brm4_Genus_best <- fs17_brm4_Genus_best[abs(fs17_brm4_Genus_best$`rowSums(fs17_brm4_Genus_bestNet)`)> 1,]
write.csv(fs17_brm4_Genus_best, file="fs17_brm4_Genus_best.csv")

fs17_brm4IA_best_melt <- melt(fs17_brm4IA_best)
fs17_brm4IB_best_melt <- melt(fs17_brm4IB_best)
fs17_brm4IC_best_melt <- melt(fs17_brm4IC_best)
fs17_brm4IAB_best_melt <- merge(fs17_brm4IA_best_melt, fs17_brm4IB_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs17_brm4IABC_best_melt <- merge(fs17_brm4IAB_best_melt, fs17_brm4IC_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs17_brm4IABC_best_melt[is.na(fs17_brm4IABC_best_melt)] <- 0


fs17_brm4IABC_best_genus <- merge(fs17_brm4IABC_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
fs17_brm4IABC_best_genus$Var1 <- fs17_brm4IABC_best_genus$Genus
colnames(fs17_brm4IABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1")

fs17_brm4IABC_best_genus <- merge(fs17_brm4IABC_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(fs17_brm4IABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1", "Genus2")
fs17_brm4IABC_best_genus$Var2 <- fs17_brm4IABC_best_genus$Genus
fs17_brm4IABC_best_genus <- fs17_brm4IABC_best_genus[,-c(1,2)]


fs17_brm4I_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3)~Genus1 + Genus2,fs17_brm4IABC_best_genus,sum)
fs17_brm4I_Genus_best <- cbind(fs17_brm4I_Genus_best, rowSums(fs17_brm4I_Genus_best[,3:5]))
colnames(fs17_brm4I_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Sum")
fs17_brm4I_Genus_best <- fs17_brm4I_Genus_best[fs17_brm4I_Genus_best$Sum != 0,]
fs17_brm4I_Genus_bestNet <- fs17_brm4I_Genus_best[,3:5]
fs17_brm4I_Genus_bestNet[fs17_brm4I_Genus_bestNet < 0] <- -1
fs17_brm4I_Genus_bestNet[fs17_brm4I_Genus_bestNet > 0] <- 1
fs17_brm4I_Genus_best <- cbind(fs17_brm4I_Genus_best, rowSums(fs17_brm4I_Genus_bestNet))
fs17_brm4I_Genus_best <- fs17_brm4I_Genus_best[abs(fs17_brm4I_Genus_best$`rowSums(fs17_brm4I_Genus_bestNet)`)> 1,]
write.csv(fs17_brm4I_Genus_best, file="fs17_brm4I_Genus_best.csv")



fs18_brm3A_best_melt <- melt(fs18_brm3A_best)
fs18_brm3B_best_melt <- melt(fs18_brm3B_best)
fs18_brm3C_best_melt <- melt(fs18_brm3C_best)
fs18_brm3AB_best_melt <- merge(fs18_brm3A_best_melt, fs18_brm3B_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs18_brm3ABC_best_melt <- merge(fs18_brm3AB_best_melt, fs18_brm3C_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs18_brm3ABC_best_melt[is.na(fs18_brm3ABC_best_melt)] <- 0

fs18_brm3ABC_best_genus <- merge(fs18_brm3ABC_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
fs18_brm3ABC_best_genus$Var1 <- fs18_brm3ABC_best_genus$Genus
colnames(fs18_brm3ABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1")

fs18_brm3ABC_best_genus <- merge(fs18_brm3ABC_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(fs18_brm3ABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1", "Genus2")
fs18_brm3ABC_best_genus$Var2 <- fs18_brm3ABC_best_genus$Genus
fs18_brm3ABC_best_genus <- fs18_brm3ABC_best_genus[,-c(1,2)]


fs18_brm3_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3)~Genus1 + Genus2,fs18_brm3ABC_best_genus,sum)
fs18_brm3_Genus_best <- cbind(fs18_brm3_Genus_best, rowSums(fs18_brm3_Genus_best[,3:5]))
colnames(fs18_brm3_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Sum")
fs18_brm3_Genus_best <- fs18_brm3_Genus_best[fs18_brm3_Genus_best$Sum != 0,]
fs18_brm3_Genus_bestNet <- fs18_brm3_Genus_best[,3:5]
fs18_brm3_Genus_bestNet[fs18_brm3_Genus_bestNet < 0] <- -1
fs18_brm3_Genus_bestNet[fs18_brm3_Genus_bestNet > 0] <- 1
fs18_brm3_Genus_best <- cbind(fs18_brm3_Genus_best, rowSums(fs18_brm3_Genus_bestNet))
fs18_brm3_Genus_best <- fs18_brm3_Genus_best[abs(fs18_brm3_Genus_best$`rowSums(fs18_brm3_Genus_bestNet)`)> 1,]
write.csv(fs18_brm3_Genus_best, file="fs18_brm3_Genus_best.csv")

fs18_brm4A_best_melt <- melt(fs18_brm4A_best)
fs18_brm4B_best_melt <- melt(fs18_brm4B_best)
fs18_brm4C_best_melt <- melt(fs18_brm4C_best)
fs18_brm4AB_best_melt <- merge(fs18_brm4A_best_melt, fs18_brm4B_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs18_brm4ABC_best_melt <- merge(fs18_brm4AB_best_melt, fs18_brm4C_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs18_brm4ABC_best_melt[is.na(fs18_brm4ABC_best_melt)] <- 0


fs18_brm4ABC_best_genus <- merge(fs18_brm4ABC_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
fs18_brm4ABC_best_genus$Var1 <- fs18_brm4ABC_best_genus$Genus
colnames(fs18_brm4ABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1")

fs18_brm4ABC_best_genus <- merge(fs18_brm4ABC_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(fs18_brm4ABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1", "Genus2")
fs18_brm4ABC_best_genus$Var2 <- fs18_brm4ABC_best_genus$Genus
fs18_brm4ABC_best_genus <- fs18_brm4ABC_best_genus[,-c(1,2)]


fs18_brm4_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3)~Genus1 + Genus2,fs18_brm4ABC_best_genus,sum)
fs18_brm4_Genus_best <- cbind(fs18_brm4_Genus_best, rowSums(fs18_brm4_Genus_best[,3:5]))
colnames(fs18_brm4_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Sum")
fs18_brm4_Genus_best <- fs18_brm4_Genus_best[fs18_brm4_Genus_best$Sum != 0,]
fs18_brm4_Genus_bestNet <- fs18_brm4_Genus_best[,3:5]
fs18_brm4_Genus_bestNet[fs18_brm4_Genus_bestNet < 0] <- -1
fs18_brm4_Genus_bestNet[fs18_brm4_Genus_bestNet > 0] <- 1
fs18_brm4_Genus_best <- cbind(fs18_brm4_Genus_best, rowSums(fs18_brm4_Genus_bestNet))
fs18_brm4_Genus_best <- fs18_brm4_Genus_best[abs(fs18_brm4_Genus_best$`rowSums(fs18_brm4_Genus_bestNet)`)> 1,]
write.csv(fs18_brm4_Genus_best, file="fs18_brm4_Genus_best.csv")

fs18_brm4IA_best_melt <- melt(fs18_brm4IA_best)
fs18_brm4IB_best_melt <- melt(fs18_brm4IB_best)
fs18_brm4IC_best_melt <- melt(fs18_brm4IC_best)
fs18_brm4IAB_best_melt <- merge(fs18_brm4IA_best_melt, fs18_brm4IB_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs18_brm4IABC_best_melt <- merge(fs18_brm4IAB_best_melt, fs18_brm4IC_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs18_brm4IABC_best_melt[is.na(fs18_brm4IABC_best_melt)] <- 0

fs18_brm4IABC_best_genus <- merge(fs18_brm4IABC_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
fs18_brm4IABC_best_genus$Var1 <- fs18_brm4IABC_best_genus$Genus
colnames(fs18_brm4IABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1")

fs18_brm4IABC_best_genus <- merge(fs18_brm4IABC_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(fs18_brm4IABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1", "Genus2")
fs18_brm4IABC_best_genus$Var2 <- fs18_brm4IABC_best_genus$Genus
fs18_brm4IABC_best_genus <- fs18_brm4IABC_best_genus[,-c(1,2)]


fs18_brm4I_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3)~Genus1 + Genus2,fs18_brm4IABC_best_genus,sum)
fs18_brm4I_Genus_best <- cbind(fs18_brm4I_Genus_best, rowSums(fs18_brm4I_Genus_best[,3:5]))
colnames(fs18_brm4I_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Sum")
fs18_brm4I_Genus_best <- fs18_brm4I_Genus_best[fs18_brm4I_Genus_best$Sum != 0,]
fs18_brm4I_Genus_bestNet <- fs18_brm4I_Genus_best[,3:5]
fs18_brm4I_Genus_bestNet[fs18_brm4I_Genus_bestNet < 0] <- -1
fs18_brm4I_Genus_bestNet[fs18_brm4I_Genus_bestNet > 0] <- 1
fs18_brm4I_Genus_best <- cbind(fs18_brm4I_Genus_best, rowSums(fs18_brm4I_Genus_bestNet))
fs18_brm4I_Genus_best <- fs18_brm4I_Genus_best[abs(fs18_brm4I_Genus_best$`rowSums(fs18_brm4I_Genus_bestNet)`)> 1,]
write.csv(fs18_brm4I_Genus_best, file="fs18_brm4I_Genus_best.csv")

fs18_brm5A_best_melt <- melt(fs18_brm5A_best)
fs18_brm5B_best_melt <- melt(fs18_brm5B_best)
fs18_brm5C_best_melt <- melt(fs18_brm5C_best)
fs18_brm5AB_best_melt <- merge(fs18_brm5A_best_melt, fs18_brm5B_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs18_brm5ABC_best_melt <- merge(fs18_brm5AB_best_melt, fs18_brm5C_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
fs18_brm5ABC_best_melt[is.na(fs18_brm5ABC_best_melt)] <- 0


fs18_brm5ABC_best_genus <- merge(fs18_brm5ABC_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
fs18_brm5ABC_best_genus$Var1 <- fs18_brm5ABC_best_genus$Genus
colnames(fs18_brm5ABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1")

fs18_brm5ABC_best_genus <- merge(fs18_brm5ABC_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(fs18_brm5ABC_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1", "Genus2")
fs18_brm5ABC_best_genus$Var2 <- fs18_brm5ABC_best_genus$Genus
fs18_brm5ABC_best_genus <- fs18_brm5ABC_best_genus[,-c(1,2)]


fs18_brm5_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3)~Genus1 + Genus2,fs18_brm5ABC_best_genus,sum)
fs18_brm5_Genus_best <- cbind(fs18_brm5_Genus_best, rowSums(fs18_brm5_Genus_best[,3:5]))
colnames(fs18_brm5_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Sum")
fs18_brm5_Genus_best <- fs18_brm5_Genus_best[fs18_brm5_Genus_best$Sum != 0,]
fs18_brm5_Genus_bestNet <- fs18_brm5_Genus_best[,3:5]
fs18_brm5_Genus_bestNet[fs18_brm5_Genus_bestNet < 0] <- -1
fs18_brm5_Genus_bestNet[fs18_brm5_Genus_bestNet > 0] <- 1
fs18_brm5_Genus_best <- cbind(fs18_brm5_Genus_best, rowSums(fs18_brm5_Genus_bestNet))
fs18_brm5_Genus_best <- fs18_brm5_Genus_best[abs(fs18_brm5_Genus_best$`rowSums(fs18_brm5_Genus_bestNet)`)> 1,]
write.csv(fs18_brm5_Genus_best, file="fs18_brm5_Genus_best.csv")


AR_1_best_melt <- melt(AR_1_best)
AR_2_best_melt <- melt(AR_2_best)
AR_3_best_melt <- melt(AR_3_best)
AR_4_best_melt <- melt(AR_4_best)
AR_12_best_melt <- merge(AR_1_best_melt, AR_2_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
AR_34_best_melt <- merge(AR_3_best_melt, AR_4_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
AR_1234_best_melt <- merge(AR_12_best_melt, AR_34_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
AR_1234_best_melt[is.na(AR_1234_best_melt)] <- 0

AR_1234_best_genus <- merge(AR_1234_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
AR_1234_best_genus$Var1 <- AR_1234_best_genus$Genus
colnames(AR_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Rep4", "Genus1")

AR_1234_best_genus <- merge(AR_1234_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(AR_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Rep4", "Genus1", "Genus2")
AR_1234_best_genus$Var2 <- AR_1234_best_genus$Genus
AR_1234_best_genus <- AR_1234_best_genus[,-c(1,2)]


AR_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3, Rep4)~Genus1 + Genus2,AR_1234_best_genus,sum)
AR_Genus_best <- cbind(AR_Genus_best, rowSums(AR_Genus_best[,3:6]))
colnames(AR_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Rep4", "Sum")
AR_Genus_best <- AR_Genus_best[AR_Genus_best$Sum != 0,]
AR_Genus_bestNet <- AR_Genus_best[,3:6]
AR_Genus_bestNet[AR_Genus_bestNet < 0] <- -1
AR_Genus_bestNet[AR_Genus_bestNet > 0] <- 1
AR_Genus_best <- cbind(AR_Genus_best, rowSums(AR_Genus_bestNet))
AR_Genus_best <- AR_Genus_best[abs(AR_Genus_best$`rowSums(AR_Genus_bestNet)`)> 1,]
write.csv(AR_Genus_best, file="AR_Genus_best.csv")


ARM_1_best_melt <- melt(ARM_1_best)
ARM_2_best_melt <- melt(ARM_2_best)
ARM_3_best_melt <- melt(ARM_3_best)
ARM_12_best_melt <- merge(ARM_1_best_melt, ARM_2_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
ARM_1234_best_melt <- merge(ARM_12_best_melt, ARM_3_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
ARM_1234_best_melt[is.na(ARM_1234_best_melt)] <- 0
#ARM_1234_best_melt <- ARM_1234_best_melt[,-5]

ARM_1234_best_genus <- merge(ARM_1234_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
ARM_1234_best_genus$Var1 <- ARM_1234_best_genus$Genus
colnames(ARM_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1")

ARM_1234_best_genus <- merge(ARM_1234_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(ARM_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Genus1", "Genus2")
ARM_1234_best_genus$Var2 <- ARM_1234_best_genus$Genus
ARM_1234_best_genus <- ARM_1234_best_genus[,-c(1,2)]


ARM_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3)~Genus1 + Genus2,ARM_1234_best_genus,sum)
ARM_Genus_best <- cbind(ARM_Genus_best, rowSums(ARM_Genus_best[,3:5]))
colnames(ARM_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Sum")
ARM_Genus_best <- ARM_Genus_best[ARM_Genus_best$Sum != 0,]
ARM_Genus_bestNet <- ARM_Genus_best[,3:5]
ARM_Genus_bestNet[ARM_Genus_bestNet < 0] <- -1
ARM_Genus_bestNet[ARM_Genus_bestNet > 0] <- 1
ARM_Genus_best <- cbind(ARM_Genus_best, rowSums(ARM_Genus_bestNet))  
ARM_Genus_best <- ARM_Genus_best[abs(ARM_Genus_best$`rowSums(ARM_Genus_bestNet)`)> 1,]
write.csv(ARM_Genus_best, file="ARM_Genus_best.csv")

HS_1234_best_melt <- melt(HS_2_best)
HS_1234_best_melt[is.na(HS_1234_best_melt)] <- 0

HS_1234_best_genus <- merge(HS_1234_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU")
HS_1234_best_genus$Var1 <- HS_1234_best_genus$Genus
colnames(HS_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Genus1")

HS_1234_best_genus <- merge(HS_1234_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU")
colnames(HS_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Genus1", "Genus2")
HS_1234_best_genus$Var2 <- HS_1234_best_genus$Genus
HS_1234_best_genus <- HS_1234_best_genus[,-c(1,2)]
write.csv(HS_1234_best_genus, file="HS_Genus_best.csv")


HG_1_best_melt <- melt(HG_1_best)
HG_2_best_melt <- melt(HG_2_best)
HG_3_best_melt <- melt(HG_3_best)
HG_4_best_melt <- melt(HG_4_best)
HG_12_best_melt <- merge(HG_1_best_melt, HG_2_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
HG_34_best_melt <- merge(HG_3_best_melt, HG_4_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
HG_1234_best_melt <- merge(HG_12_best_melt, HG_34_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
HG_1234_best_melt[is.na(HG_1234_best_melt)] <- 0

HG_1234_best_genus <- merge(HG_1234_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
HG_1234_best_genus$Var1 <- HG_1234_best_genus$Genus
colnames(HG_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Rep4", "Genus1")

HG_1234_best_genus <- merge(HG_1234_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(HG_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Rep4", "Genus1", "Genus2")
HG_1234_best_genus$Var2 <- HG_1234_best_genus$Genus
HG_1234_best_genus <- HG_1234_best_genus[,-c(1,2)]


HG_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3, Rep4)~Genus1 + Genus2,HG_1234_best_genus,sum)
HG_Genus_best <- cbind(HG_Genus_best, rowSums(HG_Genus_best[,3:6]))
colnames(HG_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Rep4", "Sum")
HG_Genus_best <- HG_Genus_best[HG_Genus_best$Sum != 0,]
HG_Genus_bestNet <- HG_Genus_best[,3:6]
HG_Genus_bestNet[HG_Genus_bestNet < 0] <- -1
HG_Genus_bestNet[HG_Genus_bestNet > 0] <- 1
HG_Genus_best <- cbind(HG_Genus_best, rowSums(HG_Genus_bestNet))
HG_Genus_best <- HG_Genus_best[abs(HG_Genus_best$`rowSums(HG_Genus_bestNet)`)> 1,]
write.csv(HG_Genus_best, file="HG_Genus_best.csv")

ES_1_best_melt <- melt(ES_1_best)
ES_2_best_melt <- melt(ES_2_best)
ES_3_best_melt <- melt(ES_3_best)
ES_4_best_melt <- melt(ES_4_best)
ES_12_best_melt <- merge(ES_1_best_melt, ES_2_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
ES_34_best_melt <- merge(ES_3_best_melt, ES_4_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
ES_1234_best_melt <- merge(ES_12_best_melt, ES_34_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
ES_1234_best_melt[is.na(ES_1234_best_melt)] <- 0

ES_1234_best_genus <- merge(ES_1234_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
ES_1234_best_genus$Var1 <- ES_1234_best_genus$Genus
colnames(ES_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Rep4", "Genus1")

ES_1234_best_genus <- merge(ES_1234_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(ES_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Rep4", "Genus1", "Genus2")
ES_1234_best_genus$Var2 <- ES_1234_best_genus$Genus
ES_1234_best_genus <- ES_1234_best_genus[,-c(1,2)]


ES_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3, Rep4)~Genus1 + Genus2,ES_1234_best_genus,sum)
ES_Genus_best <- cbind(ES_Genus_best, rowSums(ES_Genus_best[,3:6]))
colnames(ES_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Rep4", "Sum")
ES_Genus_best <- ES_Genus_best[ES_Genus_best$Sum != 0,]
ES_Genus_bestNet <- ES_Genus_best[,3:6]
ES_Genus_bestNet[ES_Genus_bestNet < 0] <- -1
ES_Genus_bestNet[ES_Genus_bestNet > 0] <- 1
ES_Genus_best <- cbind(ES_Genus_best, rowSums(ES_Genus_bestNet))
ES_Genus_best <- ES_Genus_best[abs(ES_Genus_best$`rowSums(ES_Genus_bestNet)`)> 1,]
write.csv(ES_Genus_best, file="ES_Genus_best.csv")


BRM4IS_1_best_melt <- melt(BRM4IS_1_best)
BRM4IS_2_best_melt <- melt(BRM4IS_2_best)
BRM4IS_3_best_melt <- melt(BRM4IS_3_best)
BRM4IS_4_best_melt <- melt(BRM4IS_4_best)
BRM4IS_12_best_melt <- merge(BRM4IS_1_best_melt, BRM4IS_2_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
BRM4IS_34_best_melt <- merge(BRM4IS_3_best_melt, BRM4IS_4_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
BRM4IS_1234_best_melt <- merge(BRM4IS_12_best_melt, BRM4IS_34_best_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), all=TRUE)
BRM4IS_1234_best_melt[is.na(BRM4IS_1234_best_melt)] <- 0

BRM4IS_1234_best_genus <- merge(BRM4IS_1234_best_melt, genera_nonun, by.x = "Var1", by.y = "OTU", all=TRUE)
BRM4IS_1234_best_genus$Var1 <- BRM4IS_1234_best_genus$Genus
colnames(BRM4IS_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Rep4", "Genus1")

BRM4IS_1234_best_genus <- merge(BRM4IS_1234_best_genus, genera_nonun, by.x = "OTU2", by.y = "OTU", all=TRUE)
colnames(BRM4IS_1234_best_genus) <- c("OTU1", "OTU2", "Rep1", "Rep2", "Rep3", "Rep4", "Genus1", "Genus2")
BRM4IS_1234_best_genus$Var2 <- BRM4IS_1234_best_genus$Genus
BRM4IS_1234_best_genus <- BRM4IS_1234_best_genus[,-c(1,2)]


BRM4IS_Genus_best <- aggregate(cbind(Rep1, Rep2, Rep3, Rep4)~Genus1 + Genus2,BRM4IS_1234_best_genus,sum)
BRM4IS_Genus_best <- cbind(BRM4IS_Genus_best, rowSums(BRM4IS_Genus_best[,3:6]))
colnames(BRM4IS_Genus_best) <- c("Genus1", "Genus2", "Rep1", "Rep2", "Rep3", "Rep4", "Sum")
BRM4IS_Genus_best <- BRM4IS_Genus_best[BRM4IS_Genus_best$Sum != 0,]
BRM4IS_Genus_bestNet <- BRM4IS_Genus_best[,3:6]
BRM4IS_Genus_bestNet[BRM4IS_Genus_bestNet < 0] <- -1
BRM4IS_Genus_bestNet[BRM4IS_Genus_bestNet > 0] <- 1
BRM4IS_Genus_best <- cbind(BRM4IS_Genus_best, rowSums(BRM4IS_Genus_bestNet))
BRM4IS_Genus_best <- BRM4IS_Genus_best[abs(BRM4IS_Genus_best$`rowSums(BRM4IS_Genus_bestNet)`)> 1,]
write.csv(BRM4IS_Genus_best, file="BRM4IS_Genus_best.csv")


#===================================== check this



MCC_time <- subset(MCC_time, FS == "A")

MCC_time <- subset(MCC_time, Experiment == "2")

MCC_time <- MCC_time[!is.na(MCC_time$MCC),]

anova <- aov(MCC ~ Time, data = MCC_time)
summary(anova)

tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, tukey)

Tk <- group_by(MCC_time, Time) %>%
  summarise(mean=mean(MCC), quant = quantile(MCC, probs = 0.75)) %>%
  arrange(desc(mean))

cld <- as.data.frame.list(cld$Time)
Tk$cld <- cld$Letters

MCC_time$Time <- factor(MCC_time$Time, levels = c("Early", "Mid", "Late"))

ggplot(MCC_time, aes(x=Time, y=MCC)) + geom_boxplot()+ 
  labs(x="Culture Time Frame", y="Mean Cross Correlation of OTU Abundances")+ theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data = Tk, aes(x = Time, y = quant, label = cld))



rich_time <- subset(rich_time, Experiment == "2")


anova_rich <- aov(Shared_Richness ~ Time + Media, data = rich_time)
summary(anova_rich)
tukey_rich <- TukeyHSD(anova_rich)
cld_rich <- multcompLetters4(anova_rich, tukey_rich)

Tk_rich <- group_by(rich_time, Time) %>%
  summarise(mean=mean(Shared_Richness), quant = quantile(Shared_Richness, probs = 0.75)) %>%
  arrange(desc(mean))

cld_rich <- as.data.frame.list(cld_rich$Time)
Tk_rich$cld <- cld_rich$Letters

rich_time$Time <- factor(rich_time$Time, levels = c("Early", "Mid", "Late"))
rich_time$Media <- factor(rich_time$Media, levels = c("GP", "AM", "5MM", "GM", "AP", "NS", "B", "MM", "ILC", "IL"))


late_time <- subset(rich_time, Time == "Late")

late_comp <- list(c("GP", "AM"),  c("AP", "B"))


nb_model_rich <- glm.nb(Shared_Richness ~ Media, data = late_time)
summary(nb_model_rich)

library(emmeans)
emmeans_results <- emmeans(nb_model_rich, pairwise ~ Media, type = "response")

late_plot <- ggplot(late_time, aes(x=Media, y=Shared_Richness)) + geom_dotplot( aes(x=Media, y=Shared_Richness, fill = Media), binaxis = "y", binwidth = 4)+ 
  labs(x="Medium", y="Shared Richness") +
  theme(axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=22), plot.title = element_text(size=22),
        legend.title=element_text(size=22), legend.text=element_text(size=22), strip.text.x = element_text(size = 22))+
  stat_compare_means(comparisons = late_comp, method = "t.test", label = "p.signif", size = 6)+
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)

ggsave("late_plot.svg", late_plot, "svg", width=14, height = 4.25, units = "in")

rich_plot <- ggplot(rich_time, aes(x=Time, y=Shared_Richness, fill=Time)) + geom_dotplot( aes(x=Time, y=Shared_Richness, fill=Time), binaxis = "y", binwidth = 2)+ 
  labs(x="Time", y="Shared Richness") + facet_grid(. ~ Media) + 
  theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=20), plot.title = element_text(size=20),
        legend.title=element_text(size=20), legend.text=element_text(size=15), strip.text.x = element_text(size = 20))+
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)

mcc_plot <- ggplot(rich_time, aes(x=Time, y=MCC, fill=Time)) + geom_dotplot(aes(x=Time, y=MCC, fill=Time), binaxis="y", binwidth = 0.02)+ labs(x="Time", y= "Mean Cross Correlation") +
  facet_grid(. ~ Media) + theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=20), plot.title = element_text(size=20),
                                legend.title=element_text(size=20), legend.text=element_text(size=15), strip.text.x = element_text(size = 20))+
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)

limitsperf <- ggarrange(rich_plot, mcc_plot, nrow=2, labels = c("A", "B"), font.label = list(size = 24))

ggsave("Fig.S5.svg", limitsperf, "svg", width=14, height = 8, units = "in")

SR_res_aov <- aov(Shared_Richness ~ Media + Time + Media*Time,
                  data = rich_time
)

MCC_res_aov <- aov(MCC ~ Media + Time + Media*Time,
                   data = rich_time
)


MCC_tukey <- TukeyHSD(MCC_res_aov)
SR_tukey <- TukeyHSD(SR_res_aov)

#edit
x <- which(names(rich_time) == "Media") # name of grouping variable
y <- which(
  names(rich_time) == "Shared_Richness" # names of variables to test
)
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("GP", "AM"), c("GP", "GM"), c("GP", "AP"), c("GP", "B"), c("GP", "ILC"), c("AM","GM"), c("AM", "AP"), c("AM", "B"), c("AM", "ILC"), c("GM", "AP"), c("GM", "B"), c("GM", "ILC"), c("AP", "B"), c("AP", "ILC"), c("B", "ILC")) # comparisons for post-hoc tests
# Edit until here

# Edit at your own risk
library(ggpubr)
for (i in y) {
  for (j in x) {
    p <- ggboxplot(rich_time,
                   x = colnames(rich_time[j]), y = colnames(rich_time[i]),
                   color = colnames(rich_time[j]),
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                             method = method1, label.y = max(rich_time[, i], na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}

ggplot(rich_time, aes(x=MCC, y=Shared_Richness, fill=Time)) + geom_point(binaxis="y", binwidth = 0.05)

# --------------------------------------------------------DESEQ Analysis ----------------------------------------------------------------

nonzeroEven_exp1 <- subset_samples(nonzeroEven1, Media %in% experiment1media)
nonzeroEven_exp1 <- subset_samples(nonzeroEven_exp1, Time <= 168)

fsA_exp1 <- subset_samples(nonzeroEven_exp1, FS == "FSA")
fsb_exp1 <- subset_samples(nonzeroEven_exp1, FS == "FSB")

exp1_FSA.ord <- ordinate(fsA_exp1, "PCoA", "bray")
FSA_bray = plot_ordination(fsA_exp1, exp1_FSA.ord, type="samples", color="Media", title="Beta-Diversity of FSA by Media")
exp1_FSB.ord <- ordinate(fsb_exp1, "PCoA","bray")
FSB_bray = plot_ordination(fsb_exp1, exp1_FSB.ord, type="samples", color="Media", title="Beta-Diversity of FSB by Media")

group.colors <- c(IL = "#7CAE00", MM = "#00BFC4", NS ="#C77CFF" )

FSA_bray <- FSA_bray + geom_point(size = 6, aes(color=Media))+ scale_colour_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF")) +theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                       panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=16), plot.title = element_text(size=16),
                                                                                                                                       legend.title=element_text(size=16), legend.text=element_text(size=16)) 

FSB_bray <- FSB_bray + geom_point(size = 6)+theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=16), plot.title = element_text(size=16),
                                                  legend.title=element_text(size=16), legend.text=element_text(size=16)) 

sample_data(nonzeroEven1)$Time <- factor(sample_data(nonzeroEven1)$Time)
nonzeroEven1 <- subset_samples(nonzeroEven1, Media %in% experiment2media_ILC)

nonzeroEven1 = subset_samples(nonzeroEven1, Time != "NA")
otu_table(nonzeroEven1) <- otu_table(nonzeroEven1) + 1
deseq2_Media <- phyloseq_to_deseq2(nonzeroEven1, ~Media + Time + Media:Time)
deseq2_Media$Media <- factor(deseq2_Media$Media, ordered= FALSE)
deseq2_Media$Media <- relevel(deseq2_Media$Media, ref = "B")

peptidelevels <- c("Low", "Med", "High")
nonzeroEven1_peptide <- subset_samples(nonzeroEven1, Peptide %in% peptidelevels)
deseq2_peptide <- phyloseq_to_deseq2(nonzeroEven1_peptide, ~Peptide + Peptide:Time)
deseq2_peptide$Peptide <- factor(deseq2_peptide$Peptide, ordered= FALSE)
deseq2_peptide$Peptide <- relevel(deseq2_peptide$Peptide, ref = "Med")
deseq2_Peptide <- DESeq2::DESeq(deseq2_peptide, test="LRT", modelMatrixType = "standard", reduced = ~Peptide)
DESeq2::resultsNames(deseq2_Peptide)

HighPeptide_Results <- DESeq2::results(deseq2_Peptide, name = "Peptide_High_vs_Med", cooksCutoff = FALSE)
LowPeptide_Results <- DESeq2::results(deseq2_Peptide, name = "Peptide_Low_vs_Med", cooksCutoff = FALSE)

alpha=0.01
High_sigtab = HighPeptide_Results[which(HighPeptide_Results$padj < alpha), ]
High_sigtab = cbind(as(High_sigtab, "data.frame"), as(tax_table(nonzeroEven1)[rownames(High_sigtab), ], "matrix"))

High_DeseqOut = tapply(High_sigtab$log2FoldChange, High_sigtab$Genus, function(x) max(x))
High_DeseqOut = sort(High_DeseqOut, TRUE)
High_sigtab$Genus = factor(as.character(High_sigtab$Genus), levels=names(High_DeseqOut))
High_sigtab1 <- High_sigtab %>% filter(abs(log2FoldChange) > 2 & pvalue < 0.05)
ggplot(High_sigtab1, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("Differential Abundance with High Peptide")

alpha=0.01
Low_sigtab = LowPeptide_Results[which(LowPeptide_Results$padj < alpha), ]
Low_sigtab = cbind(as(Low_sigtab, "data.frame"), as(tax_table(nonzeroEven1)[rownames(Low_sigtab), ], "matrix"))

Low_DeseqOut = tapply(Low_sigtab$log2FoldChange, Low_sigtab$Genus, function(x) max(x))
Low_DeseqOut = sort(Low_DeseqOut, TRUE)
Low_sigtab$Genus = factor(as.character(Low_sigtab$Genus), levels=names(Low_DeseqOut))
Low_sigtab1 <- Low_sigtab %>% filter(abs(log2FoldChange) > 2 & pvalue < 0.05)
ggplot(Low_sigtab1, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("Differential Abundance with Low Peptide")

deseq2_Media <- DESeq2::DESeq(deseq2_Media, test="LRT", fitType = "local", reduced = ~Media)
residual_Media <- mcols(deseq2_Media)$dispGeneEst - mcols(deseq2_Media)$dispFit

DESeq2::resultsNames(deseq2_Media)

ARM_HG_Results <- DESeq2::results(deseq2_Media, name = "Media_AM_vs_GM", cooksCutoff = TRUE)

ARMHG_sigtab = ARM_HG_Results[which(ARM_HG_Results$padj < 0.05), ]
ARMHG_sigtab = cbind(as(ARMHG_sigtab, "data.frame"), as(tax_table(nonzeroEven1)[rownames(ARMHG_sigtab), ], "matrix"))
head(ARMHG_sigtab)
dim(ARMHG_sigtab)

ARMHG_DeseqOut = tapply(ARMHG_sigtab$log2FoldChange, ARMHG_sigtab$Genus, function(x) max(x))
ARMHG_DeseqOut = sort(ARMHG_DeseqOut, TRUE)
ARMHG_sigtab$Genus = factor(as.character(ARMHG_sigtab$Genus), levels=names(ARMHG_DeseqOut))
ARMHG_sigtab1 <- ARMHG_sigtab %>% filter(abs(log2FoldChange) > 2 & pvalue < 0.05)
ggplot(ARMHG_sigtab1, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("Differential Abundance with Arabinogalactan Monomer Overrepresentation Compared to Glucose")

AR_Results <- DESeq2::results(deseq2_Media, name = "Media_AP_vs_B", cooksCutoff = FALSE)
ARM_Results <- DESeq2::results(deseq2_Media, name = "Media_AM_vs_B", cooksCutoff = FALSE)
HS_Results <- DESeq2::results(deseq2_Media, name = "Media_GP_vs_B", cooksCutoff = FALSE)
HG_Results <- DESeq2::results(deseq2_Media, name = "Media_GM_vs_B", cooksCutoff = FALSE)
BRM4IS_Results <- DESeq2::results(deseq2_Media, name = "Media_ILC_vs_B", cooksCutoff = FALSE)

alpha = 0.01
AR_sigtab = AR_Results[which(AR_Results$padj < alpha), ]
AR_sigtab = cbind(as(AR_sigtab, "data.frame"), as(tax_table(nonzeroEven1)[rownames(AR_sigtab), ], "matrix"))
head(AR_sigtab)
dim(AR_sigtab)

AR_DeseqOut = tapply(AR_sigtab$log2FoldChange, AR_sigtab$Genus, function(x) max(x))
AR_DeseqOut = sort(AR_DeseqOut, TRUE)
AR_sigtab$Genus = factor(as.character(AR_sigtab$Genus), levels=names(AR_DeseqOut))
AR_sigtab1 <- AR_sigtab %>% filter(abs(log2FoldChange) > 2 & padj < 0.05)
ggplot(AR_sigtab1, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("Differential Abundance with Arabinogalactan Overrepresentation")

ARM_sigtab = ARM_Results[which(ARM_Results$padj < alpha), ]
ARM_sigtab = cbind(as(ARM_sigtab, "data.frame"), as(tax_table(nonzeroEven1)[rownames(ARM_sigtab), ], "matrix"))
ARM_DeseqOut = tapply(ARM_sigtab$log2FoldChange, ARM_sigtab$Genus, function(x) max(x))
ARM_DeseqOut = sort(ARM_DeseqOut, TRUE)
ARM_sigtab$Genus = factor(as.character(ARM_sigtab$Genus), levels=names(ARM_DeseqOut))
ARM_sigtab1 <- ARM_sigtab %>% filter(abs(log2FoldChange) > 2 & padj < 0.05)
ggplot(ARM_sigtab1, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("Differential Abundance with Arabinogalactan Monomer Overrepresentation")

HS_sigtab = HS_Results[which(HS_Results$padj < alpha), ]
HS_sigtab = cbind(as(HS_sigtab, "data.frame"), as(tax_table(nonzeroEven1)[rownames(HS_sigtab), ], "matrix"))
HS_DeseqOut = tapply(HS_sigtab$log2FoldChange, HS_sigtab$Genus, function(x) max(x))
HS_DeseqOut = sort(HS_DeseqOut, TRUE)
HS_sigtab$Genus = factor(as.character(HS_sigtab$Genus), levels=names(HS_DeseqOut))
HS_sigtab1 <- HS_sigtab %>% filter(abs(log2FoldChange) > 2 & padj < 0.05)
ggplot(HS_sigtab1, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("Differential Abundance with Soluble Starch Overrepresentation")

HG_sigtab = HG_Results[which(HG_Results$padj < alpha), ]
HG_sigtab = cbind(as(HG_sigtab, "data.frame"), as(tax_table(nonzeroEven1)[rownames(HG_sigtab), ], "matrix"))
HG_DeseqOut = tapply(HG_sigtab$log2FoldChange, HG_sigtab$Genus, function(x) max(x))
HG_DeseqOut = sort(HG_DeseqOut, TRUE)
HG_sigtab$Genus = factor(as.character(HG_sigtab$Genus), levels=names(HG_DeseqOut))
HG_sigtab1 <- HG_sigtab %>% filter(abs(log2FoldChange) > 2 & padj < 0.05)
ggplot(HG_sigtab1, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("Differential Abundance with Glucose Overrepresentation")

BRM4IS_sigtab = BRM4IS_Results[which(BRM4IS_Results$padj < alpha), ]
BRM4IS_sigtab = cbind(as(BRM4IS_sigtab, "data.frame"), as(tax_table(nonzeroEven1)[rownames(BRM4IS_sigtab), ], "matrix"))
BRM4IS_DeseqOut = tapply(BRM4IS_sigtab$log2FoldChange, BRM4IS_sigtab$Genus, function(x) max(x))
BRM4IS_DeseqOut = sort(BRM4IS_DeseqOut, TRUE)
BRM4IS_sigtab$Genus = factor(as.character(BRM4IS_sigtab$Genus), levels=names(BRM4IS_DeseqOut))
BRM4IS_sigtab1 <- BRM4IS_sigtab %>% filter(abs(log2FoldChange) > 2 & padj < 0.05)
ggplot(BRM4IS_sigtab1, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("Differential Abundance with greater carbohydrate richness with BRM4IS")

Basal_OTU_ranking <- subset_samples(nonzeroEven1, Media == "B")
Basal_OTU_order <- rowMeans(otu_table(Basal_OTU_ranking))
Basal_OTU_order1 <- cbind(rownames(otu_table(Basal_OTU_ranking)), Basal_OTU_order)
colnames(Basal_OTU_order1) <- c("Taxonomy", "Count")
Basal_OTU_order1 <- as.data.frame(Basal_OTU_order1)
Basal_OTU_order1$Count <- as.numeric(Basal_OTU_order1$Count) 
Basal_OTU_order <- Basal_OTU_order1 %>% arrange(desc(Count))

ARM_sigtab2 <- ARM_sigtab1[,c(2,6,13)]
AR_sigtab2 <- AR_sigtab1[,c(2,6,13)]
HG_sigtab2 <- HG_sigtab1[,c(2,6,13)]
HS_sigtab2 <- HS_sigtab1[,c(2,6,13)]
BRM4IS_sigtab2 <- BRM4IS_sigtab1[,c(2,6,13)]

AR_ARM_sigtab <- merge(AR_sigtab2, ARM_sigtab2, by.x = "unique", by.y = "unique", all =TRUE)
colnames(AR_ARM_sigtab) <- c("unique", "log2FC_AR", "padj_AR", "log2FC_ARM", "padj_ARM")
AR_ARM_HG <- merge(AR_ARM_sigtab, HG_sigtab2, by.x = "unique", by.y = "unique", all = TRUE)
colnames(AR_ARM_HG) <- c("unique", "log2FC_AR", "padj_AR", "log2FC_ARM", "padj_ARM", "log2FC_HG", "padj_HG")
AR_ARM_HG_HS <- merge(AR_ARM_HG, HS_sigtab2, by.x = "unique", by.y = "unique", all = TRUE)
colnames(AR_ARM_HG_HS) <- c("unique", "log2FC_AR", "padj_AR", "log2FC_ARM", "padj_ARM", "log2FC_HG", "padj_HG", "log2FC_HS", "padj_HS")
All_Comp_B <- merge(AR_ARM_HG_HS, BRM4IS_sigtab2, by.x = "unique", by.y= "unique", all = TRUE)
colnames(All_Comp_B) <- c("unique", "log2FC_AR", "padj_AR", "log2FC_ARM", "padj_ARM", "log2FC_HG", "padj_HG", "log2FC_HS", "padj_HS", "log2FC_ILC", "padj_ILC")

All_Comp_B_log2FC <- All_Comp_B[,c(1,2,4,6,8,10)] 
colnames(All_Comp_B_log2FC) <- c("unique", "AP", "AM", "GM", "GP", "ILC")
All_Comp_B_log2FC_melt <- melt(All_Comp_B_log2FC, id = "unique")
All_Comp_B_log2FC_melt <- merge(All_Comp_B_log2FC_melt , Basal_OTU_order, by.x = "unique", by.y = "Taxonomy")
All_Comp_B_log2FC_melt <- All_Comp_B_log2FC_melt %>% arrange(desc(Count))

All_Comp_B_log2FC_melt[All_Comp_B_log2FC_melt == "k__Bacteria_p__Firmicutes_c__Clostridia_o__Clostridiales_f__Lachnospiraceae_g__Ruminococcus"] <- "f__Lachnospiraceae_g__Ruminococcus"
All_Comp_B_log2FC_melt[All_Comp_B_log2FC_melt == "k__Bacteria_p__Firmicutes_c__Clostridia_o__Clostridiales_f__Ruminococcaceae_g__Ruminococcus"] <- "f__Ruminococcaceae_g__Ruminococcus"

DESEQ_Media_v_Basal <- ggplot(All_Comp_B_log2FC_melt, aes(x= variable, y= reorder(unique, Count, decreasing=T), fill= value)) + 
  geom_tile()+ scale_fill_gradientn("value", colours=c("blue", "white", "red"), limits = c(-6,6), na.value = "transparent") + scale_x_discrete(limits = c("ILC", "AP", "GM", "GP", "AM"))+
  xlab("Medium") + ylab("Differentially Abundant Genera")+
  theme(axis.text.x = element_text(size=14), legend.key.size= unit(1, "cm"), axis.text.y = element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=14), legend.text=element_text(size=12))

deseq2_Media_HS <- phyloseq_to_deseq2(nonzeroEven1, ~Media + Time + Media:Time)
deseq2_Media_HS$Media <- factor(deseq2_Media_HS$Media, ordered= FALSE)
deseq2_Media_HS$Media <- relevel(deseq2_Media_HS$Media, ref = "GM")
deseq2_Media_HS <- DESeq2::DESeq(deseq2_Media_HS, test="LRT", fitType = "local", reduced = ~Media)

HG_v_HS_Results <- DESeq2::results(deseq2_Media_HS, name = "Media_GP_vs_GM", cooksCutoff = FALSE)

HGHS_sigtab = HG_v_HS_Results[which(HG_v_HS_Results$padj < alpha), ]
HGHS_sigtab = cbind(as(HGHS_sigtab, "data.frame"), as(tax_table(nonzeroEven1)[rownames(HGHS_sigtab), ], "matrix"))
HGHS_DeseqOut = tapply(HGHS_sigtab$log2FoldChange, HGHS_sigtab$Genus, function(x) max(x))
HGHS_DeseqOut = sort(HGHS_DeseqOut, TRUE)
HGHS_sigtab$Genus = factor(as.character(HGHS_sigtab$Genus), levels=names(HGHS_DeseqOut))
HGHS_sigtab1 <- HGHS_sigtab %>% filter(abs(log2FoldChange) > 2 & padj < 0.05)
HGHS_sigtab1 <- cbind(HGHS_sigtab1, rep("GP", times=nrow(HGHS_sigtab1)))

HGHS_sigtab2 <- HGHS_sigtab1[,c(2,6,13,14)]
colnames(HGHS_sigtab2) <- c("log2FC", "padj", "unique", "Media")
HGHS_sigtab2$Media <- rep("GM vs GP", nrow(HGHS_sigtab2))

ggplot(HGHS_sigtab2, aes(x=Media, y=unique, fill=log2FC)) + geom_tile() + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("Differential Abundance Against ")

deseq2_Media_AR <- phyloseq_to_deseq2(nonzeroEven1, ~Media + Time + Media:Time)
deseq2_Media_AR$Media <- factor(deseq2_Media_AR$Media, ordered= FALSE)
deseq2_Media_AR$Media <- relevel(deseq2_Media_AR$Media, ref = "AM")
deseq2_Media_AR <- DESeq2::DESeq(deseq2_Media_AR, test="LRT", fitType = "local", reduced = ~Media)

AR_v_ARM_Results <- DESeq2::results(deseq2_Media_AR, name = "Media_AP_vs_AM", cooksCutoff = FALSE)

ARARM_sigtab = AR_v_ARM_Results[which(AR_v_ARM_Results$padj < alpha), ]
ARARM_sigtab = cbind(as(ARARM_sigtab, "data.frame"), as(tax_table(nonzeroEven1)[rownames(ARARM_sigtab), ], "matrix"))
ARARM_DeseqOut = tapply(ARARM_sigtab$log2FoldChange, ARARM_sigtab$Genus, function(x) max(x))
ARARM_DeseqOut = sort(ARARM_DeseqOut, TRUE)
ARARM_sigtab$Genus = factor(as.character(ARARM_sigtab$Genus), levels=names(ARARM_DeseqOut))
ARARM_sigtab1 <- ARARM_sigtab %>% filter(abs(log2FoldChange) > 2 & padj < 0.05)
ARARM_sigtab1 <- cbind(ARARM_sigtab1, rep("AP vs AM", times=nrow(ARARM_sigtab1)))

ARARM_sigtab2 <- ARARM_sigtab1[,c(2,6,13,14)]
colnames(ARARM_sigtab2) <- c("log2FC", "padj", "unique", "Media")

ggplot(ARARM_sigtab1, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("Differential Abundance with between Arabinogalactan Monomer and Arabinogalactan Overrepresentation")

All_Comp <- merge(All_Comp_B, HGHS_sigtab2, by.x = "unique", by.y= "unique", all = TRUE)
colnames(All_Comp) <- c("unique", "log2FC_AR", "padj_AR", "log2FC_ARM", "padj_ARM", "log2FC_HG", "padj_HG", "log2FC_HS", "padj_HS", "log2FC_ILC", "padj_ILC", "log2FC_GMvsGP", "padj_GMvsGP", "Media")

All_Comparisons <- merge(All_Comp, ARARM_sigtab2, by.x =c("unique"), by.y = c("unique"), all = TRUE)

All_Comparisons <- All_Comparisons[,-c(14,17)]

colnames(All_Comparisons) <- c("unique", "log2FC_AR", "padj_AR", "log2FC_ARM", "padj_ARM", "log2FC_HG", "padj_HG", "log2FC_HS", "padj_HS", "log2FC_ILC", "padj_ILC", "log2FC_GMvsGP", "padj_GMvsGP", "log2FC_APvsAM", "padj_APvsAM")

All_Comp_log2FC <- All_Comparisons[,c(1,2,4,6,8,10,12,14)]
All_Comp_padj <- All_Comparisons[,c(1,3,5,7,9,11,13,15)]
colnames(All_Comp_log2FC) <- c("unique", "AP", "AM", "GM", "GP", "ILC", "GP vs GM", "AP vs AM")
colnames(All_Comp_padj) <- c("unique", "AP", "AM", "GM", "GP", "ILC", "GP vs GM", "AP vs AM")

All_Comp_log2FC_melt <- melt(All_Comp_log2FC, id = "unique")
All_Comp_padj_melt <- melt(All_Comp_padj, id = "unique")
colnames(All_Comp_log2FC_melt) <- c("unique", "Comparison", "log2FC")
colnames(All_Comp_padj_melt) <- c("unique", "Comparison", "padj")
All_Comp_melt <- merge(All_Comp_log2FC_melt , All_Comp_padj_melt, by = c("unique","Comparison"))

All_Comp_melt <- merge(All_Comp_melt , Basal_OTU_order, by.x = "unique", by.y = "Taxonomy")

All_Comp_melt <- All_Comp_melt %>% arrange(desc(Count))

All_Comp_melt[All_Comp_melt == "k__Bacteria_p__Firmicutes_c__Clostridia_o__Clostridiales_f__Lachnospiraceae_g__Ruminococcus"] <- "f__Lachnospiraceae_g__Ruminococcus"
All_Comp_melt[All_Comp_melt == "k__Bacteria_p__Firmicutes_c__Clostridia_o__Clostridiales_f__Ruminococcaceae_g__Ruminococcus"] <- "f__Ruminococcaceae_g__Ruminococcus"

All_Comp_melt$facet <- ifelse(All_Comp_melt$Comparison %in% c("AM", "AP", "GM", "GP", "ILC"), "vs Basal", All_Comp_melt$Comparison)

All_Comp_melt <- filter(All_Comp_melt, padj < 1e-5 )

DESEQ_plot <- ggplot(All_Comp_padj, aes(x=Comparison, y= reorder(unique, Count, decreasing=T)))
DESEQ_plot <- DESEQ_plot + geom_count(aes(color = log2FC, size = padj)) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
  theme_bw()+ scale_size_area(trans = 'log', max_size = 4, breaks = c(1e-85, 1e-65, 1e-45, 1e-25, 1e-10), labels = expression(1e-85, 1e-65, 1e-45, 1e-25, 1e-10))+
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=30),
        axis.text.y = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.title=element_text(size=20), legend.text=element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        line = element_blank())

DESEQ_plot$labels$colour <- "Log 2 Fold Change"
DESEQ_plot$labels$size <- "Adjusted P-value"

DESEQ_plot <- DESEQ_plot + facet_grid(~factor(facet,  levels = c("vs Basal", "6", "7")), space = "free", scales = "free_x") + theme(strip.text.x = element_blank())

DESEQ_Media <- ggplot(All_Comp_B_log2FC_melt, aes(x= variable, y= reorder(unique, Count, decreasing=T), fill= value)) + 
  geom_tile()+ scale_fill_gradientn("value", colours=c("blue", "white", "red"), limits = c(-6,6), na.value = "transparent") + scale_x_discrete(limits = c("ILC", "AP", "GM", "GP", "AM"))+
  xlab("Medium") + ylab("Differentially Abundant Genera")+
  theme(axis.text.x = element_text(size=14), legend.key.size= unit(1, "cm"), axis.text.y = element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=14), legend.text=element_text(size=12))

ggsave("DESEQ_5.pdf", DESEQ_plot, "pdf", width=14, height = 8, units = "in")

# =================================================== mccplots ==========================================================================

library(ggplot2)
library(grid)
library(gridExtra)
library(png)

late_table <- grid.table(late_comp)

node_table <- data.table::as.data.table(node_table, rownames = FALSE)

colnames(node_table) <- c("Node Number", "Taxonomy")

nodes <- tableGrob(node_table, rows = NULL)

fig5_lay <- rbind(
  c(1, 1, 1, 1),
  c(2, 3, 7, 7),
  c(2, 3, 7, 7),
  c(4, 5, 7, 7),
  c(4, 5, 7, 7),
  c(6, 6, 7, 7),
  c(6, 6, 7, 7)
)

fig5_grobs <- list(rasterGrob(late_plot), rasterGrob(AM_network), rasterGrob(GM_network), rasterGrob(AP_network), rasterGrob(B_network), rasterGrob(ILC_network), nodes)

# Specify the scale factor to shrink the second PNG
scale_factor <- 0.6

# Wrap the second grob in a gtable object
scaled_raster_grob <- rasterGrob(AM_network, 
                                 width = unit(scale_factor, "npc"), 
                                 height = unit(scale_factor, "npc"))



# Replace the modified gtable in the fig5_grobs list
fig5_grobs[[2]] <- scaled_raster_grob

figure5 <- grid.arrange(grobs = fig5_grobs, ncol = 4,
                        layout_matrix = fig5_lay,
                        labels = c("A", "B", "C", "D", "E", "F", "G"),
                        widths = c(1, 1, 1, 0.5),
                        heights = unit(c(2, 1, 1, 1, 1, 1, 1), c("null", rep("null", 5))))

layout_matrix <- matrix(c(
  1, 1, 1, 1,
  1, 1, 1, 1,# Row 1: Grob A spans 4 columns
  2, 3, 7, 7,       # Row 2: Grob B in Col 1, Grob C in Col 2, Table Grob in Col 3 and 4
  2, 3, 7, 7,       # Row 3: Same as Row 2
  4, 5, 7, 7,     # Row 4: Grob D spans 2 columns
  4, 5, 7, 7,     # Row 5: Grob E spans 2 columns
  6, 6, 7, 7,     # Row 6: Grob F spans 2 columns
  6, 6, 7, 7      # Row 7: Same as Row 6
), nrow = 8, byrow = TRUE)

labelA <- textGrob("A", x = unit(0.05, "npc"), y = unit(0.95, "npc"), just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
labelB <- textGrob("B", x = unit(0.05, "npc"), y = unit(0.95, "npc"), just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
labelC <- textGrob("C", x = unit(0.05, "npc"), y = unit(0.95, "npc"), just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
labelD <- textGrob("D", x = unit(0.05, "npc"), y = unit(0.95, "npc"), just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
labelE <- textGrob("E", x = unit(0.05, "npc"), y = unit(0.95, "npc"), just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
labelF <- textGrob("F", x = unit(0.05, "npc"), y = unit(0.95, "npc"), just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))

# Combine labels and Grobs into a new list
combinedGrobs <- list(
  grid.arrange(labelA, rasterGrob(late_plot), ncol = 4),  # Grob A with label
  grid.arrange(labelB, rasterGrob(AM_network,
                                  width = unit(scale_factor, "npc"), 
                                  height = unit(scale_factor, "npc")), ncol = 1),  # Grob B with label
  grid.arrange(labelC, rasterGrob(GM_network), ncol = 1),  # Grob C with label
  grid.arrange(labelD, rasterGrob(AP_network), ncol = 1),  # Grob D with label
  grid.arrange(labelE, rasterGrob(B_network), ncol = 1),  # Grob E with label
  grid.arrange(labelF, rasterGrob(ILC_network), ncol = 1),  # Grob F with label
  nodes                             # Table Grob (no label)
)

combinedGrobs <- list(
  ggplot() + 
    annotation_custom(rasterGrob(late_plot)) + 
    annotation_custom(labelA)+ theme_minimal(),  # Grob A with label
  ggplot() + 
    annotation_custom(rasterGrob(AM_network,
                                 width = unit(scale_factor, "npc"), 
                                 height = unit(scale_factor, "npc")))+
    annotation_custom(labelB)+ theme_minimal(),  # Grob B with label
  ggplot() + 
    annotation_custom(rasterGrob(GM_network)) + 
    annotation_custom(labelC) + theme_minimal(),  # Grob C with label
  ggplot() + 
    annotation_custom(rasterGrob(AP_network)) + 
    annotation_custom(labelD) + theme_minimal(),  # Grob D with label
  ggplot() + 
    annotation_custom(rasterGrob(B_network)) + 
    annotation_custom(labelE)+ theme_minimal(),  # Grob E with label
  ggplot() + 
    annotation_custom(rasterGrob(ILC_network)) + 
    annotation_custom(labelF)+ theme_minimal(),  # Grob F with label
  nodes                             # Table Grob (no label)
)

# Arrange the Grobs using ggarrange
final_plot <- grid.arrange(
  grobs = combinedGrobs,
  layout_matrix = layout_matrix
)

width <- 12
height <- 14

# Set the desired resolution in dots per inch
dpi <- 600

# Save the plot while maintaining relative proportions
ggsave("fig5_update2.pdf", plot = final_plot, 'pdf', width = width, height = height, dpi = dpi)

#plots for windows and the test non-parametrically
#mtc between windows
#switch stool to individual
#substrate plot



nonzeroEven_stable <- subset_samples(nonzeroEven, Time >= 72 & Time <=240 & Time != 192)
nonzeroEven_exp1_stable <- subset_samples(nonzeroEven_stable, Media %in% experiment1media)

exp1.ord <- ordinate(nonzeroEven_exp1_stable, "PCoA", "bray")
bray_exp1 = plot_ordination(nonzeroEven_exp1_stable, exp1.ord, type="samples", color="Media", shape = "FS")

exp1_bray <- bray_exp1 + geom_point(size=6)+ ggtitle("Sample Community Composition for Experiment 1")+ theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),plot.title = element_text(size=20, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), 
                                                                                                             legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 18))

nonzeroEven_FSA <- subset_samples(nonzeroEven_stable, FS == "FSA")
nonzeroEven_FSB<- subset_samples(nonzeroEven_stable, FS == "FSB")
nonzeroEven_FSA_exp2 <- subset_samples(nonzeroEven_FSA, Media %in% experiment2media_ILC)
nonzeroEven_FSA_exp1 <- subset_samples(nonzeroEven_FSA, Media %in% experiment1media)

FSA_exp1.ord <- ordinate(nonzeroEven_FSA_exp1, "PCoA", "bray")
FSAbray_exp1 = plot_ordination(nonzeroEven_FSA_exp1, FSA_exp1.ord, type="samples", color="Media")

FSA_exp1_bray <- FSAbray_exp1 + geom_point(size=6)+ ggtitle("FSA Sample Community Composition by Media type")+ theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),plot.title = element_text(size=20, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), 
                                                                                                                     legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 18))

FSA_exp2.ord <- ordinate(nonzeroEven_FSA_exp2, "PCoA", "bray")
FSA.ord <- ordinate(nonzeroEven_FSA, "PCoA", "bray")
FSB.ord <- ordinate(nonzeroEven_FSB, "PCoA", "bray")
FSAbray = plot_ordination(nonzeroEven_FSA, FSA.ord, type="samples", color="Media")
FSBbray = plot_ordination(nonzeroEven_FSB, FSB.ord, type="samples", color="Media")

FSA_exp2_bray <- plot_ordination(nonzeroEven_FSA_exp2, FSA_exp2.ord, type = "samples", color = "Media")

FSAbray <- FSAbray + geom_point(size=6)+ ggtitle("FSA Sample Community Composition by Media type")+ theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),plot.title = element_text(size=20, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), 
                                                                                                          legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 18))

FSA_exp2_bray <- FSA_exp2_bray + geom_point(size=6)+ ggtitle("FSA Sample Community Composition by Media type")+ theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),plot.title = element_text(size=20, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), 
                                                                                                                      legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 18))


FSBbray <- FSBbray + geom_point(size=6) + ggtitle("FSB Sample Community Composition by Media type")+ theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18), plot.title = element_text(size=20, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), 
                                                                                                           legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 18))


#Initial analysis of Carbohydrate Diversity (not CheSL)
row.names(sub_list) <- sub_list[,1]
sub_list<- sub_list[,-1]
sub_list<- t(sub_list)
shannon <- diversity(sub_list, "shannon")
evenness <- shannon/log(ncol(sub_list))

media_comp_exp1 <- data.frame( col1 = c("Media", "Concentration", "Richness", "Evenness", "Shannon Index", "Microbiota Richness"),
                               col2 = c("NS", 0.64, 5, 0.556, 1.510, 64.2),
                               col3 = c("MM", 0.90, 7, 0.584, 1.580, 66.5),
                               col4 = c("IL", 2.44, 14, 0.604, 1.640, 92.833),
                               col5 = c("5MM", 3.30, 7, .503, 1.360, 32.666))

media_comp_FSA <- data.frame( col1 = c("Media", "Concentration", "Richness", "Evenness", "Shannon Index", "Microbiota Richness"),
                              col2 = c("NS", 0.64, 5, 0.556, 1.510, 75),
                              col3 = c("MM", 0.90, 7, 0.584, 1.580, 74),
                              col4 = c("IL", 2.44, 14, 0.604, 1.640, 92),
                              col5 = c("B", 0.84, 6, 0.947, 1.7, 71),
                              col6 = c("AP", 2.74, 6, .553, 0.992, 62.5),
                              col7 = c("AM", 2.74, 7, .665, 1.290, 14.5),
                              col8 = c("GP", 2.64, 6, 0.513, 0.919, 13.75),
                              col9 = c("GM", 2.80, 6, 0.586, 1.050, 51.5),
                              col10 = c("ILC", 0.50, 11, .988, 2.370, 79.5)
)

media_comp_FSB <- data.frame( col1 = c("Media", "Concentration", "Richness", "Evenness", "Shannon Index", "Microbiota Richness"),
                              col2 = c("NS", 0.64, 5, 0.556, 1.510, 57),
                              col3 = c("MM", 0.90, 7, 0.584, 1.580, 59),
                              col4 = c("IL", 2.44, 14, 0.604, 1.640, 93.66),
                              col5 = c("5MM", 3.30, 7, .503, 1.360, 32.666))


media_comp_FSA_exp2 <- data.frame( col1 = c("Media", "Concentration", "Richness", "Evenness", "Shannon Index", "Microbiota Richness"),
                                   col2 = c("B", 0.84, 6, 0.947, 1.7, 71),
                                   col3 = c("AP", 2.74, 6, .553, 0.992, 62.5),
                                   col4 = c("AM", 2.74, 7, .665, 1.290, 14.5),
                                   col5 = c("GP", 2.64, 6, 0.513, 0.919, 13.75),
                                   col6 = c("GM", 2.80, 6, 0.586, 1.050, 51.5),
                                   col7 = c("ILC", 0.50, 11, .988, 2.370, 79.5)
)

media_comp_FSA_exp1 <- data.frame( col1 = c("Media", "Concentration", "Richness", "Evenness", "Shannon Index", "Microbiota Richness"),
                                   col2 = c("NS", 0.64, 5, 0.556, 1.510, 75),
                                   col3 = c("MM", 0.90, 7, 0.584, 1.580, 74),
                                   col4 = c("IL", 2.44, 14, 0.604, 1.640, 92))


colnames(media_comp_FSA_exp1) <- media_comp_FSA_exp1[1,]
media_comp_FSA_exp1 <- media_comp_FSA_exp1[-1,]
rownames(media_comp_FSA_exp1) <- media_comp_FSA_exp1[,1]
media_comp_FSA_exp1 <- media_comp_FSA_exp1[,-1]

media_comp_FSA_exp1_long <- media_comp_FSA_exp1 %>%
  rownames_to_column(var = "Metric") %>%
  pivot_longer(cols = -Metric, names_to = "Condition", values_to = "Value")

media_comp_FSA_exp1_long$Value <- as.numeric(media_comp_FSA_exp1_long$Value)

media_comp_FSA_exp1_long$Metric = factor(media_comp_FSA_exp1_long$Metric, levels=c('Concentration','Richness','Evenness','Shannon Index', "Microbiota Richness"))

media_comp_plot_FSA_exp1 <- ggplot(media_comp_FSA_exp1_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  theme_minimal() +
  ggtitle("Media Carbohydrate Diversity Measures and Microbiota Richness")+  xlab("Medium") + ylab("")+
  theme(axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=20, hjust = 0.5),
        legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 16))

colnames(media_comp_FSA_exp2) <- media_comp_FSA_exp2[1,]
media_comp_FSA_exp2 <- media_comp_FSA_exp2[-1,]
rownames(media_comp_FSA_exp2) <- media_comp_FSA_exp2[,1]
media_comp_FSA_exp2 <- media_comp_FSA_exp2[,-1]

media_comp_FSA_exp2_long <- media_comp_FSA_exp2 %>%
  rownames_to_column(var = "Metric") %>%
  pivot_longer(cols = -Metric, names_to = "Condition", values_to = "Value")

media_comp_FSA_exp2_long$Value <- as.numeric(media_comp_FSA_exp2_long$Value)

media_comp_FSA_exp2_long$Metric = factor(media_comp_FSA_exp2_long$Metric, levels=c('Concentration','Richness','Evenness','Shannon Index', "Microbiota Richness"))

media_comp_plot_FSA_exp2 <- ggplot(media_comp_FSA_exp2_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  theme_minimal() +
  ggtitle("Media Carbohydrate Diversity Measures and Microbiota Richness")+  xlab("Medium") + ylab("")+
  theme(axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=20, hjust = 0.5),
        legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 16))

Concentration_FSA <- ggplot(media_summary_FSA_exp2, aes(x=Concentration, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Concentration")+  xlab("Carbohydrate Concentration") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(0,4,0.5))+
  stat_cor(method="pearson", label.x = 1, label.y=30, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

DoubleRichness_FSA <- ggplot(media_summary_FSA_exp2, aes(x=Carbohydrate.Richness, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Carbohydrate Richness")+  xlab("Carbohydrate Richness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(2,18,2))+
  stat_cor(method="pearson", label.x = 8, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Even_FSA <- ggplot(media_summary_FSA_exp2, aes(x=Carbohydrate.Evenness, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Evenness")+  xlab("Carbohydrate Evenness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = .8, label.y=30, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Shannon_FSA <- ggplot(media_summary_FSA_exp2, aes(x=Carbohydrate.Shannon, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Shannon Index")+  xlab("Carbohydrate Shannon Index") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.4, label.y=30, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

FSAfirst_row <- ggarrange(FSA_exp2_bray, labels = c("A"), font.label = list(size = 24))

# Add the label "C" to the middle row
FSAsecond_row <- annotate_figure(media_comp_plot_FSA_exp2, top = text_grob("B", size = 24, hjust = 0, x=0.01,face="bold"))

# Arrange the last two plots with labels D and E
FSAthird_row <- ggarrange(Concentration_FSA, DoubleRichness_FSA, ncol = 2, labels = c("C", "D"), font.label = list(size = 24))

#Fourth
FSAfourth_row<- ggarrange(Rich_Even_FSA, Rich_Shannon_FSA, ncol = 2, labels = c("E", "F"), font.label = list(size = 24))

# Combine all three rows into the final figure
Figure_3_FSA <- ggarrange(FSAfirst_row, FSAsecond_row, FSAthird_row, FSAfourth_row, nrow = 4)

ggsave("figure3evenFSA_new.pdf", Figure_3_FSA, "pdf", width=14, height = 14, units = "in")

Linkage_shannon_FSA <- ggplot(media_summary_FSA_exp2, aes(x=LS, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by CheSL Shannon Index")+  xlab("CheSL Shannon Index") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.95, label.y=20, size=8)+
  theme(axis.text.x = element_text(size=24), axis.text.y = element_text(size=24),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=24), plot.title = element_text(size=24),
        legend.title=element_text(size=24), legend.text=element_text(size=24))


# Experiment 1 FSB ----------------------------------------------------------------------------------------------------------------------


Concentration_FSB <- ggplot(media_summary_FSB, aes(x=Concentration, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Concentration")+  xlab("Carbohydrate Concentration") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(0,4,0.5))+
  stat_cor(method="pearson", label.x = 2, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 



DoubleRichness_FSB <- ggplot(media_summary_FSB, aes(x=Carbohydrate.Richness, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Carbohydrate Richness")+  xlab("Carbohydrate Richness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(2,18,2))+
  stat_cor(method="pearson", label.x = 8, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 


Rich_Even_FSB <- ggplot(media_summary_FSB, aes(x=Carbohydrate.Evenness, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Evenness")+  xlab("Carbohydrate Evenness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = .55, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Shannon_FSB <- ggplot(media_summary_FSB, aes(x=Carbohydrate.Shannon, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Shannon Index")+  xlab("Carbohydrate Shannon Index") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.2, label.y=80, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Concentration_FSA <- ggplot(media_summary_FSA, aes(x=Concentration, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Concentration")+  xlab("Carbohydrate Concentration") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(0,4,0.5))+
  stat_cor(method="pearson", label.x = 1, label.y=30, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

DoubleRichness_FSA <- ggplot(media_summary_FSA, aes(x=Carbohydrate.Richness, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Carbohydrate Richness")+  xlab("Carbohydrate Richness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(2,18,2))+
  stat_cor(method="pearson", label.x = 8, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Even_FSA <- ggplot(media_summary_FSA, aes(x=Carbohydrate.Evenness, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Evenness")+  xlab("Carbohydrate Evenness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = .8, label.y=30, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Shannon_FSA <- ggplot(media_summary_FSA, aes(x=Carbohydrate.Shannon, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Shannon Index")+  xlab("Carbohydrate Shannon Index") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.4, label.y=30, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

# --------------------------------------------------- Experiment 1 -----------------------------------------------------------------------------
colnames(media_comp_FSA) <- media_comp_FSA[1,]
media_comp_FSA <- media_comp_FSA[-1,]
rownames(media_comp_FSA) <- media_comp_FSA[,1]
media_comp_FSA <- media_comp_FSA[,-1]

media_comp_FSA_long <- media_comp_FSA %>%
  rownames_to_column(var = "Metric") %>%
  pivot_longer(cols = -Metric, names_to = "Condition", values_to = "Value")

media_comp_FSA_long$Value <- as.numeric(media_comp_FSA_long$Value)

media_comp_FSA_long$Metric = factor(media_comp_FSA_long$Metric, levels=c('Concentration','Richness','Evenness','Shannon Index', "Microbiota Richness"))

media_comp_plot_exp1 <- ggplot(media_comp_exp1_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  theme_minimal() +
  ggtitle("Media Carbohydrate Diversity Measures and Microbiota Richness")+  xlab("Medium") + ylab("")+
  theme(axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=20, hjust = 0.5), strip.text = element_text(size = 14),
        legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 14))

Concentration_exp1 <- ggplot(media_summary_exp1, aes(x=Concentration, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Concentration")+  xlab("Carbohydrate Concentration") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(0,4,0.5))+
  stat_cor(method="pearson", label.x = 2, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 


DoubleRichness_exp1 <- ggplot(media_summary_exp1, aes(x=Carbohydrate.Richness, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Carbohydrate Richness")+  xlab("Carbohydrate Richness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(2,18,2))+
  stat_cor(method="pearson", label.x = 8, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 


Rich_Even_exp1 <- ggplot(media_summary_exp1, aes(x=Carbohydrate.Evenness, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Evenness")+  xlab("Carbohydrate Evenness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = .55, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Shannon_exp1 <- ggplot(media_summary_exp1, aes(x=Carbohydrate.Shannon, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Shannon Index")+  xlab("Carbohydrate Shannon Index") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.2, label.y=80, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Linkage_richness_exp1 <- ggplot(media_summary_exp1, aes(x=LR, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by CheSL Richness")+  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.25)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 7, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_shannon_exp1 <- ggplot(media_summary_exp1, aes(x=LS, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by CheSL Shannon Index")+  xlab("CheSL Shannon Index") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.6, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_evenness_exp1 <- ggplot(media_summary_exp1, aes(x=LE, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by CheSLEvenness")+  xlab("CheSL Evenness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.55, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

exp1first_row <- ggarrange(FSA_bray, labels = c("A"), font.label = list(size = 24))

# Add the label "C" to the middle row
exp1second_row <- annotate_figure(media_comp_plot_exp1, top = text_grob("B", size = 24, hjust = 0, x=0.01,face="bold"))

# Arrange the last two plots with labels D and E
exp1third_row <- ggarrange(Concentration_exp1, DoubleRichness_exp1, ncol = 2, labels = c("C", "D"), font.label = list(size = 24))

#Fourth
exp1fourth_row<- ggarrange(Rich_Even_exp1, Rich_Shannon_exp1, ncol = 2, labels = c("E", "F"), font.label = list(size = 24))

exp1fifth_row <- ggarrange(Linkage_richness_exp1, Linkage_shannon_exp1, Linkage_evenness_exp1, ncol = 3, labels = c("G", "H", "I"), font.label = list(size = 24))

# Combine all three rows into the final figure
Figure_S2_FSA <- ggarrange(FSAfirst_row, FSAsecond_row, FSAthird_row, FSAfourth_row, FSAfifth_row, nrow = 5)

ggsave("figureS2evenFSA_new1.svg", Figure_S2_FSA, "svg", width=14, height = 20, units = "in")

#---------------------------------------------------- Just correlations FSA and FSB Experiment 1--------------------------------


media_comp_plot_exp1 <- ggplot(media_comp_exp1_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  theme_minimal() +
  ggtitle("Media Carbohydrate Diversity Measures and Microbiota Richness")+  xlab("Medium") + ylab("")+
  theme(axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=20, hjust = 0.5), strip.text = element_text(size = 14),
        legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 14))

Concentration_exp1 <- ggplot(media_summary_exp1, aes(x=Concentration, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Concentration")+  xlab("Carbohydrate Concentration") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(0,4,0.5))+
  stat_cor(method="pearson", label.x = 2, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 


DoubleRichness_exp1 <- ggplot(media_summary_exp1, aes(x=Carbohydrate.Richness, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Carbohydrate Richness")+  xlab("Carbohydrate Richness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(2,18,2))+
  stat_cor(method="pearson", label.x = 8, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 


Rich_Even_exp1 <- ggplot(media_summary_exp1, aes(x=Carbohydrate.Evenness, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Evenness")+  xlab("Carbohydrate Evenness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = .55, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Shannon_exp1 <- ggplot(media_summary_exp1, aes(x=Carbohydrate.Shannon, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Shannon Index")+  xlab("Carbohydrate Shannon Index") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.2, label.y=80, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Linkage_richness_exp1 <- ggplot(media_summary_exp1, aes(x=LR, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Richness")+  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.25)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 7, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_shannon_exp1 <- ggplot(media_summary_exp1, aes(x=LS, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Shannon Index")+  xlab("CheSL Shannon Index") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.5, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_evenness_exp1 <- ggplot(media_summary_exp1, aes(x=LE, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Evenness")+  xlab("CheSL Evenness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.55, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

exp1first_row <- ggarrange(exp1_bray, labels = c("A"), font.label = list(size = 24))

# Add the label "C" to the middle row
exp1second_row <- annotate_figure(media_comp_plot_exp1, top = text_grob("B", size = 24, hjust = 0, x=0.01,face="bold"))

# Arrange the last two plots with labels D and E
exp1third_row <- ggarrange(Concentration_exp1, DoubleRichness_exp1, ncol = 2, labels = c("C", "D"), font.label = list(size = 24))

#Fourth
exp1fourth_row<- ggarrange(Rich_Even_exp1, Rich_Shannon_exp1, ncol = 2, labels = c("E", "F"), font.label = list(size = 24))

exp1fifth_row <- ggarrange(Linkage_richness_exp1, Linkage_shannon_exp1, Linkage_evenness_exp1, ncol = 3, labels = c("G", "H", "I"), font.label = list(size = 24))

# Combine all three rows into the final figure
Figure_S2_exp1 <- ggarrange(exp1first_row, exp1second_row, exp1third_row, exp1fourth_row, exp1fifth_row, nrow = 5)

ggsave("figureS2even_new1.pdf", Figure_S2_exp1, "pdf", width=14, height = 20, units = "in")

# -----------------------------------------------------FSB full --------------------------------------------------------------------------------
media_comp_FSB <- data.frame( col1 = c("Media", "Concentration", "Richness", "Evenness", "Shannon Index", "Microbiota Richness"),
                              col2 = c("NS", 0.64, 5, 0.556, 1.510, 57),
                              col3 = c("MM", 0.90, 7, 0.584, 1.580, 59),
                              col4 = c("IL", 2.44, 14, 0.604, 1.640, 93.66),
                              col5 = c("5MM", 3.30, 7, .503, 1.360, 32.666))

colnames(media_comp_FSB) <- media_comp_FSB[1,]
media_comp_FSB <- media_comp_FSB[-1,]
rownames(media_comp_FSB) <- media_comp_FSB[,1]
media_comp_FSB <- media_comp_FSB[,-1]

media_comp_FSB_long <- media_comp_FSB %>%
  rownames_to_column(var = "Metric") %>%
  pivot_longer(cols = -Metric, names_to = "Condition", values_to = "Value")

media_comp_FSB_long$Value <- as.numeric(media_comp_FSB_long$Value)

media_comp_FSB_long$Metric = factor(media_comp_FSB_long$Metric, levels=c('Concentration','Richness','Evenness','Shannon Index', "Microbiota Richness"))

media_comp_plot_FSB <- ggplot(media_comp_FSB_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  theme_minimal() +
  ggtitle("Media Carbohydrate Diversity Measures and Microbiota Richness")+  xlab("Medium") + ylab("")+
  theme(axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=20, hjust = 0.5), strip.text = element_text(size = 14),
        legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 18))

media_comp_plot_FSB <- ggplot(media_comp_FSB_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  theme_minimal() +
  ggtitle("Media Carbohydrate Diversity Measures and Microbiota Richness")+  xlab("Medium") + ylab("")+
  theme(axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=20, hjust = 0.5), strip.text = element_text(size = 14),
        legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 14))

Concentration_FSB <- ggplot(media_summary_FSB, aes(x=Concentration, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Concentration")+  xlab("Carbohydrate Concentration") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(0,4,0.5))+
  stat_cor(method="pearson", label.x = 2, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

DoubleRichness_FSB <- ggplot(media_summary_FSB, aes(x=Carbohydrate.Richness, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Carbohydrate Richness")+  xlab("Carbohydrate Richness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(2,18,2))+
  stat_cor(method="pearson", label.x = 8, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Even_FSB <- ggplot(media_summary_FSB, aes(x=Carbohydrate.Evenness, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Evenness")+  xlab("Carbohydrate Evenness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = .55, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Shannon_FSB <- ggplot(media_summary_FSB, aes(x=Carbohydrate.Shannon, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Shannon Index")+  xlab("Carbohydrate Shannon Index") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.2, label.y=80, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Linkage_richness_FSB <- ggplot(media_summary_FSB, aes(x=LR, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Richness")+  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.25)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 9, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_shannon_FSB <- ggplot(media_summary_FSB, aes(x=LS, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Shannon Index")+  xlab("CheSL Shannon Index") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.4, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_evenness_FSB <- ggplot(media_summary_FSB, aes(x=LE, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Evenness")+  xlab("CheSL Evenness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.55, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

FSB_first_row <- ggarrange(FSBbray, labels = c("A"), font.label = list(size = 24))

FSB_second_row <- annotate_figure(media_comp_plot_FSB, top = text_grob("B", size = 24, hjust = 0, x=0.01,face="bold"))

# Arrange the last two plots with labels D and E
FSB_third_row <- ggarrange(Concentration_FSB, DoubleRichness_FSB, ncol = 2, labels = c("C", "D"), font.label = list(size = 24))

#Fourth
FSB_fourth_row<- ggarrange(Rich_Even_FSB, Rich_Shannon_FSB, ncol = 2, labels = c("E", "F"), font.label = list(size = 24))

FSB_fifth_row <- ggarrange(Linkage_richness_FSB, Linkage_shannon_FSB, Linkage_evenness_FSB, ncol = 3, labels = c("G", "H", "I"), font.label = list(size = 24))

# Combine all three rows into the final figure
Figure_S1_FSB <- ggarrange(FSB_first_row, FSB_second_row, FSB_third_row, FSB_fourth_row, FSB_fifth_row, nrow = 5)

ggsave("figureS1_FSB_new.svg", Figure_S1_FSB, "svg", width=14, height = 20, units = "in")

# _--------------------------------------------------------- FSA Experiment 1 -------------------------------------------------------

media_comp_plot_exp1 <- ggplot(media_comp_FSA_exp1_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  theme_minimal() +
  ggtitle("Media Carbohydrate Diversity Measures and Microbiota Richness")+  xlab("Medium") + ylab("")+
  theme(axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=20, hjust = 0.5), strip.text = element_text(size = 14),
        legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 14))

Concentration_FSA_exp1 <- ggplot(media_summary_FSA_exp1, aes(x=Concentration, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Concentration")+  xlab("Carbohydrate Concentration") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(0,4,0.5))+
  stat_cor(method="pearson", label.x = 1, label.y=60, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

DoubleRichness_FSA_exp1 <- ggplot(media_summary_FSA_exp1, aes(x=Carbohydrate.Richness, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Carbohydrate Richness")+  xlab("Carbohydrate Richness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(2,18,2))+
  stat_cor(method="pearson", label.x = 8, label.y=60, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Even_FSA_exp1 <- ggplot(media_summary_FSA_exp1, aes(x=Carbohydrate.Evenness, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Evenness")+  xlab("Carbohydrate Evenness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = .55, label.y=60, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Shannon_FSA_exp1 <- ggplot(media_summary_FSA_exp1, aes(x=Carbohydrate.Shannon, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Shannon Index")+  xlab("Carbohydrate Shannon Index") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.5, label.y=90, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Linkage_richness_FSA_exp1 <- ggplot(media_summary_FSA_exp1, aes(x=LR, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Richness")+  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.25)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 10, label.y=50, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_shannon_FSA_exp1 <- ggplot(media_summary_FSA_exp1, aes(x=LS, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Shannon Index")+  xlab("CheSL Shannon Index") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.6, label.y=50, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_evenness_FSA_exp1 <- ggplot(media_summary_FSA_exp1, aes(x=LE, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSLEvenness")+  xlab("CheSL Evenness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.65, label.y=50, size=6)+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

FSAexp1first_row <- ggarrange(FSA_exp1_bray, labels = c("A"), font.label = list(size = 24))

# Add the label "C" to the middle row
FSAexp1second_row <- annotate_figure(media_comp_plot_FSA_exp1, top = text_grob("B", size = 24, hjust = 0, x=0.01,face="bold"))

# Arrange the last two plots with labels D and E
FSAexp1third_row <- ggarrange(Concentration_FSA_exp1, DoubleRichness_FSA_exp1, ncol = 2, labels = c("C", "D"), font.label = list(size = 24))

#Fourth
FSAexp1fourth_row<- ggarrange(Rich_Even_FSA_exp1, Rich_Shannon_FSA_exp1, ncol = 2, labels = c("E", "F"), font.label = list(size = 24))

FSAexp1fifth_row <- ggarrange(Linkage_richness_FSA_exp1, Linkage_shannon_FSA_exp1, Linkage_evenness_FSA_exp1, ncol = 3, labels = c("G", "H", "I"), font.label = list(size = 24))

# Combine all three rows into the final figure
Figure_S3_exp1 <- ggarrange(FSAexp1first_row, FSAexp1second_row, FSAexp1third_row, FSAexp1fourth_row, FSAexp1fifth_row, nrow = 5)

ggsave("figureS2_FSA_Exp1_11.pdf", Figure_S3_exp1, "pdf", width=14, height = 20, units = "in")

#------------------------------------------------------------- cheSL Final Plots -----------------------------------------------------------------


experiment2media_ILC <- c("AM", "GM", "GP", "AP", "ILC", "B")
experiment1media <- c("NS", "MM", "5MM", "IL")
experiment1media_noIL <- c("NS", "MM", "5MM")

media_summary_FSB <- subset(media_summary, FS == "FSB")
media_summary_FSA <- subset(media_summary, FS == "FSA")
media_summary_FSA_exp1 <- subset(media_summary_FSA, Media %in% experiment1media)
media_summary_FSB_exp1 <- subset(media_summary_FSB, Media %in% experiment1media)
media_summary_FSAnoIL <- subset(media_summary_FSA, Media != "IL")
media_summary_exp1 <- subset(media_summary, Media %in% experiment1media)
media_summary_FSA_exp2 <- subset(media_summary, Media %in% experiment2media_ILC)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Figure S1 ~~~~~~~~~~~~~~~~~~~~~~~~~~

Linkage_richness_FSA_exp2 <- ggplot(media_summary_FSA_exp2, aes(x=LR, y=Richness.Species)) + 
  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.25)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 9, label.y=25, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_shannon_FSA_exp2 <- ggplot(media_summary_FSA_exp2, aes(x=LS, y=Richness.Species)) + 
  xlab("CheSL Shannon Index") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.4, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_evenness_FSA_exp2 <- ggplot(media_summary_FSA_exp2, aes(x=LE, y=Richness.Species)) + 
  xlab("CheSL Evenness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.45, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

all_FSB <- ggarrange(Linkage_richness_FSB,  Linkage_evenness_FSB, Linkage_shannon_FSB, ncol = 3, labels = c("A", "B", "C"), font.label = list(size = 24))
exp1_FSA <- ggarrange(Linkage_richness_FSA_exp1, Linkage_evenness_FSA_exp1, Linkage_shannon_FSA_exp1,  ncol = 3, labels = c("D", "E", "F"), font.label = list(size = 24))
exp2_FSA <- ggarrange(Linkage_richness_FSA_exp2, Linkage_evenness_FSA_exp2, Linkage_shannon_FSA_exp2,  ncol = 3, labels = c("G", "H", "I"), font.label = list(size = 24))

Figure_S1_12_23 <- ggarrange(all_FSB, exp1_FSA, exp2_FSA, nrow = 3)

ggsave("FigureS1_12_23.pdf", Figure_S1_12_23, "pdf", width=16, height = 12, units = "in")

#---------------------------------------------------- Figure 3 (Just correlations FSA and FSB Experiment 1) --------------------------------

Linkage_richness_exp1 <- ggplot(media_summary_exp1, aes(x=LR, y=Richness.Species)) + 
  xlab("") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media, shape = FS), size = 4)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 14, label.y=43, size = 4)+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=12), plot.title = element_text(size=12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))

Linkage_shannon_exp1 <- ggplot(media_summary_exp1, aes(x=LS, y=Richness.Species)) + 
  xlab("") + ylab("")+
  geom_point(aes(color = Media, shape = FS), size = 4)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.025)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.56, label.y=40, size = 4)+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=12), plot.title = element_text(size=12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))

Linkage_evenness_exp1 <- ggplot(media_summary_exp1, aes(x=LE, y=Richness.Species)) + 
  xlab("") + ylab("")+
  geom_point(aes(color = Media, shape = FS), size = 4)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.0075)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.55, label.y=45, size = 4)+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=12), plot.title = element_text(size=12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))

Linkage_richness_FSAnoIL <- ggplot(media_summary_FSAnoIL, aes(x=LR, y=Richness.Species)) + 
  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 4)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 10, label.y=30, size = 4)+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=12), plot.title = element_text(size=12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))

Linkage_shannon_FSAnoIL <- ggplot(media_summary_FSAnoIL, aes(x=LS, y=Richness.Species)) + 
  xlab("CheSL Shannon Index") + ylab("")+
  geom_point(aes(color = Media), size = 4)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  scale_y_continuous(breaks=seq(0,80,20))+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.35, label.y=20, size = 4)+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=12), plot.title = element_text(size=12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))

Linkage_evenness_FSAnoIL <- ggplot(media_summary_FSAnoIL, aes(x=LE, y=Richness.Species)) + 
  xlab("CheSL Evenness") + ylab("")+
  geom_point(aes(color = Media), size = 4)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.45, label.y=20, size = 4)+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=12), plot.title = element_text(size=12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))

exp1_bothFS <- ggarrange(Linkage_richness_exp1, Linkage_evenness_exp1, Linkage_shannon_exp1,  ncol = 3, labels = c("A", "B", "C"), font.label = list(size = 24))
all_FSAnoIL <- ggarrange(Linkage_richness_FSAnoIL,  Linkage_evenness_FSAnoIL, Linkage_shannon_FSAnoIL, ncol = 3, labels = c("D", "E", "F"), font.label = list(size = 24))

# Combine all three rows into the final figure
Figure_3_12_20 <- ggarrange(exp1_bothFS, all_FSAnoIL, nrow = 2)

ggsave("Figure3.pdf", Figure_3_12_20, "pdf", width=16, height = 8, units = "in")

Linkage_richness_FSA <- ggplot(media_summary_FSA, aes(x=LR, y=Richness.Species)) + 
  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 4)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  scale_y_continuous(breaks=seq(0,110,30))+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 10, label.y=30, size = 4)+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=12), plot.title = element_text(size=12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))

Linkage_shannon_FSA <- ggplot(media_summary_FSA, aes(x=LS, y=Richness.Species)) + 
  xlab("CheSL Shannon Index") + ylab("")+
  geom_point(aes(color = Media), size = 4)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  scale_y_continuous(breaks=seq(0,110,30))+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.35, label.y=20, size = 4)+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=12), plot.title = element_text(size=12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))

Linkage_evenness_FSA <- ggplot(media_summary_FSA, aes(x=LE, y=Richness.Species)) + 
  xlab("CheSL Evenness") + ylab("")+
  geom_point(aes(color = Media), size = 4)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  scale_y_continuous(breaks=seq(0,110,30))+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.45, label.y=20, size = 4)+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=12), plot.title = element_text(size=12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))

exp1_bothFS <- ggarrange(Linkage_richness_exp1, Linkage_evenness_exp1, Linkage_shannon_exp1,  ncol = 3, labels = c("A", "B", "C"), font.label = list(size = 24))
all_FSA <- ggarrange(Linkage_richness_FSA,  Linkage_evenness_FSA, Linkage_shannon_FSA, ncol = 3, labels = c("D", "E", "F"), font.label = list(size = 24))

# Combine all three rows into the final figure
Figure_3_12_20wIL <- ggarrange(Linkage_richness_exp1, Linkage_evenness_exp1, Linkage_shannon_exp1,Linkage_richness_FSA,  Linkage_evenness_FSA, Linkage_shannon_FSA, nrow = 3, ncol=2, labels = c("A", "B", "C", "D", "E","F"), font.label = list(size = 18))

ggsave("Figure3_withIL.pdf", Figure_3_12_20wIL, "pdf", width=7, height = 9, units = "in")

# _--------------------------------------------------------- FSA Experiment 1 -------------------------------------------------------

Linkage_richness_FSA_exp1 <- ggplot(media_summary_FSA_exp1, aes(x=LR, y=Richness.Species)) + 
  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.25)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 10, label.y=90, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_shannon_FSA_exp1 <- ggplot(media_summary_FSA_exp1, aes(x=LS, y=Richness.Species)) + 
  xlab("CheSL Shannon Index") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.6, label.y=90, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_evenness_FSA_exp1 <- ggplot(media_summary_FSA_exp1, aes(x=LE, y=Richness.Species)) + 
  xlab("CheSL Evenness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.6, label.y=60, size=6)+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

#---------------------------------------------------- Supplemental Figure 1: FSA and FSB independent CheSL correlations December 19, 2024 --------------------------------
Linkage_richness_FSAnoILs1 <- ggplot(media_summary_FSAnoIL, aes(x=LR, y=Richness.Species)) + 
  xlab("") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 10, label.y=30, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_shannon_FSAnoILs1 <- ggplot(media_summary_FSAnoIL, aes(x=LS, y=Richness.Species)) + 
  xlab("") + ylab("")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  scale_y_continuous(breaks=seq(0,80,20))+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.35, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_evenness_FSAnoILs1 <- ggplot(media_summary_FSAnoIL, aes(x=LE, y=Richness.Species)) + 
  xlab("") + ylab("")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.45, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_richness_FSB <- ggplot(media_summary_FSB, aes(x=LR, y=Richness.Species)) + 
  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 14, label.y=42, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_shannon_FSB <- ggplot(media_summary_FSB, aes(x=LS, y=Richness.Species)) + 
  xlab("CheSL Shannon Index") + ylab("")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.025)+
  scale_y_continuous(breaks=seq(0,100,20))+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.62, label.y=38, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_evenness_FSB <- ggplot(media_summary_FSB, aes(x=LE, y=Richness.Species)) + 
  xlab("CheSL Evenness") + ylab("")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.005)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.6, label.y=40, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

all_FSAs1 <- ggarrange(Linkage_richness_FSAnoILs1,  Linkage_evenness_FSAnoILs1, Linkage_shannon_FSAnoILs1, ncol = 3, labels = c("A", "B", "C"), font.label = list(size = 24))
all_FSB <- ggarrange(Linkage_richness_FSB,  Linkage_evenness_FSB, Linkage_shannon_FSB, ncol = 3, labels = c("D", "E", "F"), font.label = list(size = 24))

# Combine all three rows into the final figure
Figure_S1_12_20 <- ggarrange(all_FSAs1, all_FSB, nrow = 2)

ggsave("FigureS1.pdf", Figure_S1_12_20, "pdf", width=16, height = 8, units = "in")

# -----------------------------------------------------FSB full (Supplemental Figure 1 DEF) --------------------------------------------------------------------------------
media_comp_FSB <- data.frame( col1 = c("Media", "Concentration", "Richness", "Evenness", "Shannon Index", "Microbiota Richness"),
                              col2 = c("NS", 0.64, 5, 0.556, 1.510, 57),
                              col3 = c("MM", 0.90, 7, 0.584, 1.580, 59),
                              col4 = c("IL", 2.44, 14, 0.604, 1.640, 93.66),
                              col5 = c("5MM", 3.30, 7, .503, 1.360, 32.666))

colnames(media_comp_FSB) <- media_comp_FSB[1,]
media_comp_FSB <- media_comp_FSB[-1,]
rownames(media_comp_FSB) <- media_comp_FSB[,1]
media_comp_FSB <- media_comp_FSB[,-1]

media_comp_FSB_long <- media_comp_FSB %>%
  rownames_to_column(var = "Metric") %>%
  pivot_longer(cols = -Metric, names_to = "Condition", values_to = "Value")

media_comp_FSB_long$Value <- as.numeric(media_comp_FSB_long$Value)

media_comp_FSB_long$Metric = factor(media_comp_FSB_long$Metric, levels=c('Concentration','Richness','Evenness','Shannon Index', "Microbiota Richness"))

media_comp_plot_FSB <- ggplot(media_comp_FSB_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  theme_minimal() +
  ggtitle("Media Carbohydrate Diversity Measures and Microbiota Richness")+  xlab("Medium") + ylab("")+
  theme(axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=20, hjust = 0.5), strip.text = element_text(size = 14),
        legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 18))

media_comp_plot_FSB <- ggplot(media_comp_FSB_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  theme_minimal() +
  ggtitle("Media Carbohydrate Diversity Measures and Microbiota Richness")+  xlab("Medium") + ylab("")+
  theme(axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=20, hjust = 0.5), strip.text = element_text(size = 14),
        legend.title=element_text(size=18), legend.text=element_text(size=18), strip.text.x = element_text(size = 14))

Concentration_FSB <- ggplot(media_summary_FSB, aes(x=Concentration, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Concentration")+  xlab("Carbohydrate Concentration") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(0,4,0.5))+
  stat_cor(method="pearson", label.x = 2, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

DoubleRichness_FSB <- ggplot(media_summary_FSB, aes(x=Carbohydrate.Richness, y=Richness.Species))+
  geom_point(aes(color = Media),  size =5)+
  ggtitle("Microbiota Richness by Carbohydrate Richness")+  xlab("Carbohydrate Richness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  scale_x_continuous(breaks = seq(2,18,2))+
  stat_cor(method="pearson", label.x = 8, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Even_FSB <- ggplot(media_summary_FSB, aes(x=Carbohydrate.Evenness, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Evenness")+  xlab("Carbohydrate Evenness") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = .55, label.y=40, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Rich_Shannon_FSB <- ggplot(media_summary_FSB, aes(x=Carbohydrate.Shannon, y=Richness.Species)) + 
  geom_point(aes(color = Media),position = position_dodge(width = .1), size = 5)+
  ggtitle("Microbiota Richness by Carbohydrate Shannon Index")+  xlab("Carbohydrate Shannon Index") + ylab("Microbiota Richness")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.2, label.y=80, size=5)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18)) 

Linkage_richness_FSB <- ggplot(media_summary_FSB, aes(x=LR, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Richness")+  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.25)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 9, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_shannon_FSB <- ggplot(media_summary_FSB, aes(x=LS, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Shannon Index")+  xlab("CheSL Shannon Index") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.4, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Linkage_evenness_FSB <- ggplot(media_summary_FSB, aes(x=LE, y=Richness.Species)) + 
  ggtitle("Microbiota Richness by \n CheSL Evenness")+  xlab("CheSL Evenness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6, position = position_dodge(width = .1))+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.55, label.y=20, size=6)+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

FSB_first_row <- ggarrange(FSBbray, labels = c("A"), font.label = list(size = 24))

FSB_second_row <- annotate_figure(media_comp_plot_FSB, top = text_grob("B", size = 24, hjust = 0, x=0.01,face="bold"))

# Arrange the last two plots with labels D and E
FSB_third_row <- ggarrange(Concentration_FSB, DoubleRichness_FSB, ncol = 2, labels = c("C", "D"), font.label = list(size = 24))

#Fourth
FSB_fourth_row<- ggarrange(Rich_Even_FSB, Rich_Shannon_FSB, ncol = 2, labels = c("E", "F"), font.label = list(size = 24))

FSB_fifth_row <- ggarrange(Linkage_richness_FSB, Linkage_shannon_FSB, Linkage_evenness_FSB, ncol = 3, labels = c("G", "H", "I"), font.label = list(size = 24))

# Combine all three rows into the final figure
Figure_S1_FSB <- ggarrange(FSB_first_row, FSB_second_row, FSB_third_row, FSB_fourth_row, FSB_fifth_row, nrow = 5)

ggsave("figureS1_FSB.svg", Figure_S1_FSB, "svg", width=14, height = 20, units = "in")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
chesl_asvs <- cbind(data.frame(otu_table(physeq)), data.frame(tax_table(physeq)))

write.csv(chesl_asvs, file = "chesl_asvs.csv")

#--------------------------------------------------Supplemental 2 -----------------------------------------------



lindemann_richness <- ggplot(Lindemann_chesl, aes(x=Richness, y=Microbiota.Richness)) + 
  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = .75)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 5, label.y=32, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

lindemann_shannon <- ggplot(Lindemann_chesl, aes(x=Shannon, y=Microbiota.Richness)) + 
  xlab("CheSL Shannon Index") + ylab("")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.125)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = .25, label.y=30, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

lindemann_evenness <- ggplot(Lindemann_chesl, aes(x=Evenness, y=Microbiota.Richness)) + 
  xlab("CheSL Evenness") + ylab("")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.05)+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.2, label.y=26, size=6)+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

chung_richness <- ggplot(chung_chesl, aes(x=Richness, y=Microbiota.Richness)) + 
  xlab("CheSL Richness") + ylab("Microbiota Richness")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 1)+
  scale_fill_brewer(palette="Dark2")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 5, label.y=40, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

chung_shannon <- ggplot(chung_chesl, aes(x=Shannon, y=Microbiota.Richness)) + 
  xlab("CheSL Shannon Index") + ylab("")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.1)+
  scale_fill_brewer(palette="Dark2")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.75, label.y=40, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

chung_evenness <- ggplot(chung_chesl, aes(x=Evenness, y=Microbiota.Richness)) + 
  xlab("CheSL Evenness") + ylab("")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.015)+
  scale_fill_brewer(palette="Dark2")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.65, label.y=40, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))


# Combine all three rows into the final figure
Figure_S2_12_22 <-ggarrange(lindemann_richness, lindemann_evenness, lindemann_shannon, chung_richness, chung_evenness, chung_shannon, ncol = 3, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"), font.label = list(size = 24))

ggsave("FigureS2.pdf", Figure_S2_12_22, "pdf", width=16, height = 8, units = "in")

#----------------------------------Supplemental Figure 6-------------------------------

ILC_1 <- subset_samples(ILC_1, Time >= 84)
ILC_inv_1 <- estimate_richness(ILC_1, measures = "InvSimpson")
ILC_mInv_1 <- mean(ILC_inv_1$InvSimpson)

ILC_2 <- subset_samples(ILC_2, Time >= 84)
ILC_inv_2 <- estimate_richness(ILC_2, measures = "InvSimpson")
ILC_mInv_2 <- mean(ILC_inv_2$InvSimpson)

ILC_3 <- subset_samples(ILC_3, Time >= 84)
ILC_inv_3 <- estimate_richness(ILC_3, measures = "InvSimpson")
ILC_mInv_3 <- mean(ILC_inv_3$InvSimpson)

ILC_4 <- subset_samples(ILC_4, Time >= 84)
ILC_inv_4 <- estimate_richness(ILC_4, measures = "InvSimpson")
ILC_mInv_4 <- mean(ILC_inv_4$InvSimpson)

B_1 <- subset_samples(B_1, Time >= 84)
B_inv_1 <- estimate_richness(B_1, measures = "InvSimpson")
B_mInv_1 <- mean(B_inv_1$InvSimpson)

B_2 <- subset_samples(B_2, Time >= 84)
B_inv_2 <- estimate_richness(B_2, measures = "InvSimpson")
B_mInv_2 <- mean(B_inv_2$InvSimpson)

B_3 <- subset_samples(B_3, Time >= 84)
B_inv_3 <- estimate_richness(B_3, measures = "InvSimpson")
B_mInv_3 <- mean(B_inv_3$InvSimpson)

B_4 <- subset_samples(B_4, Time >= 84)
B_inv_4 <- estimate_richness(B_4, measures = "InvSimpson")
B_mInv_4 <- mean(B_inv_4$InvSimpson)

AP_1 <- subset_samples(AP_1, Time >= 84)
AP_inv_1 <- estimate_richness(AP_1, measures = "InvSimpson")
AP_mInv_1 <- mean(AP_inv_1$InvSimpson)

AP_2 <- subset_samples(AP_2, Time >= 84)
AP_inv_2 <- estimate_richness(AP_2, measures = "InvSimpson")
AP_mInv_2 <- mean(AP_inv_2$InvSimpson)

AP_3 <- subset_samples(AP_3, Time >= 84)
AP_inv_3 <- estimate_richness(AP_3, measures = "InvSimpson")
AP_mInv_3 <- mean(AP_inv_3$InvSimpson)

AP_4 <- subset_samples(AP_4, Time >= 84)
AP_inv_4 <- estimate_richness(AP_4, measures = "InvSimpson")
AP_mInv_4 <- mean(AP_inv_4$InvSimpson)

AM_1 <- subset_samples(AM_1, Time >= 84)
AM_inv_1 <- estimate_richness(AM_1, measures = "InvSimpson")
AM_mInv_1 <- mean(AM_inv_1$InvSimpson)

AM_2 <- subset_samples(AM_2, Time >= 84)
AM_inv_2 <- estimate_richness(AM_2, measures = "InvSimpson")
AM_mInv_2 <- mean(AM_inv_2$InvSimpson)

AM_3 <- subset_samples(AM_3, Time >= 84)
AM_inv_3 <- estimate_richness(AM_3, measures = "InvSimpson")
AM_mInv_3 <- mean(AM_inv_3$InvSimpson)

AM_4 <- subset_samples(AM_4, Time >= 84)
AM_inv_4 <- estimate_richness(AM_4, measures = "InvSimpson")
AM_mInv_4 <- mean(AM_inv_4$InvSimpson)

GM_1 <- subset_samples(GM_1, Time >= 84)
GM_inv_1 <- estimate_richness(GM_1, measures = "InvSimpson")
GM_mInv_1 <- mean(GM_inv_1$InvSimpson)

GM_2 <- subset_samples(GM_2, Time >= 84)
GM_inv_2 <- estimate_richness(GM_2, measures = "InvSimpson")
GM_mInv_2 <- mean(GM_inv_2$InvSimpson)

GM_3 <- subset_samples(GM_3, Time >= 84)
GM_inv_3 <- estimate_richness(GM_3, measures = "InvSimpson")
GM_mInv_3 <- mean(GM_inv_3$InvSimpson)

GM_4 <- subset_samples(GM_4, Time >= 84)
GM_inv_4 <- estimate_richness(GM_4, measures = "InvSimpson")
GM_mInv_4 <- mean(GM_inv_4$InvSimpson)

GP_1 <- subset_samples(GP_1, Time >= 84)
GP_inv_1 <- estimate_richness(GP_1, measures = "InvSimpson")
GP_mInv_1 <- mean(GP_inv_1$InvSimpson)

GP_2 <- subset_samples(GP_2, Time >= 84)
GP_inv_2 <- estimate_richness(GP_2, measures = "InvSimpson")
GP_mInv_2 <- mean(GP_inv_2$InvSimpson)

GP_3 <- subset_samples(GP_3, Time >= 84)
GP_inv_3 <- estimate_richness(GP_3, measures = "InvSimpson")
GP_mInv_3 <- mean(GP_inv_3$InvSimpson)

GP_4 <- subset_samples(GP_4, Time >= 84)
GP_inv_4 <- estimate_richness(GP_4, measures = "InvSimpson")
GP_mInv_4 <- mean(GP_inv_4$InvSimpson)

Inverse_simpson <- c(
  ILC_mInv_1, ILC_mInv_2, ILC_mInv_3, ILC_mInv_4, B_mInv_1, B_mInv_2, B_mInv_3, B_mInv_4,
  AP_mInv_1, AP_mInv_2, AP_mInv_3, AP_mInv_4, AM_mInv_1, AM_mInv_2, AM_mInv_3, AM_mInv_4,
  GP_mInv_1, GP_mInv_2, GP_mInv_3, GP_mInv_4, GM_mInv_1, GM_mInv_2, GM_mInv_3, GM_mInv_4
)

media_list <- c(
  rep("ILC", 4), rep("B", 4), rep("AP",4), rep("AM", 4), rep("GP", 4), rep("GM",4)
)

#Values same as in Supplementary Data file

cheSL_shannon <- c(
  rep(2.52611181, 4), rep(1.636931906, 4), rep(2.123686579,4), rep(1.190280757, 4), rep(0.86473509, 4), rep(1.424182172,4)
)

cheSL_evenness <- c(
  rep(0.857926358, 4), rep(0.638192681, 4), rep(0.827964331,4), rep(0.611683308, 4), rep(0.337135346, 4), rep(0.593930097,4)
)

cheSL_richness <- c(
  rep(19, 4), rep(13, 4), rep(13,4), rep(7, 4), rep(13, 4), rep(11,4)
)

inverse_simpson_df <- data.frame(
  media_list, Inverse_simpson, cheSL_shannon, cheSL_evenness, cheSL_richness
) 

colnames(inverse_simpson_df) <- c("Media", "InverseSimpson", "Shannon", "Evenness", "Richness")

Loss_simp$Concentration <- as.factor(Loss_simp$Concentration)

Loss_simp <- ggplot(Loss_simp, aes(x=Carbohydrate.Complexity, y=InvSimpson)) + 
  xlab("Carbohydrate Complexity") + ylab("Inverse Simpson")+
  geom_point(aes(color = Carbohydrate, shape = Concentration), size = 6)+
  stat_summary(aes(group=Carbohydrate), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.015)+
  scale_fill_brewer(palette="Dark2")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = .05, label.y=4, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

inv_shannon <- ggplot(inverse_simpson_df, aes(x=cheSL_shannon, y=Inverse_simpson)) + 
  xlab("CheSL Shannon Index") + ylab("Inverse Simpson")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.075)+
  scale_fill_brewer(palette="Dark2")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 1.5, label.y=3, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

inv_evenness <- ggplot(inverse_simpson_df, aes(x=cheSL_evenness, y=Inverse_simpson)) + 
  xlab("CheSL Evenness") + ylab("Inverse Simpson")+
  geom_point(aes(color = Media), size = 6)+
  stat_summary(aes(group=Media), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.03)+
  scale_fill_brewer(palette="Dark2")+
  geom_smooth(method="lm", se=FALSE, color="black")+
  stat_cor(method="pearson", label.x = 0.3, label.y=11, size=6)+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size=18), plot.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=18))

Figure_S6_12_22 <-ggarrange(Loss_simp, inv_evenness, inv_shannon, ncol = 3, labels = c("A", "B", "C"), font.label = list(size = 24))

ggsave("FigureS6.pdf", Figure_S6_12_22, "pdf", width=18, height = 4, units = "in")
