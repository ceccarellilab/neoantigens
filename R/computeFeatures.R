#ComputeFeatures
#Compute peptide features
#

library(readxl)
library(dplyr)
pepTableb <- read.table("neoantigens_peptides.tsv", sep="\t", header=TRUE)#read_excel("PEPTIDE_IM.xlsx")
dim(pepTableb)

###  PropHydroAro ----

library(parallel)

# Function to calculate proportions
PropHydroAro <- function(peptide) {
    require(stringr)
    hydrophobic_residues <- c('A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y', 'C')
    aromatic_residues <- c('F', 'W', 'Y')
    # Convert the peptide sequence to uppercase
    peptide <- toupper(peptide)
    # Calculate the length of the peptide
    peptide_length <- nchar(peptide)
    # Calculate the count of hydrophobic residues
    hydrophobic_count <- sum(str_count(peptide, hydrophobic_residues))
    # Calculate the count of aromatic residues
    aromatic_count <- sum(str_count(peptide, aromatic_residues))
    # Calculate the proportions
    if(aromatic_count!=0)
        result <- hydrophobic_count/aromatic_count
    else
        result <- 9
    return(result)
}

# Test the function
PropHydroAroScore <- mclapply(pepTableb$pep, PropHydroAro, mc.cores = 64)
PropHydroAroScore <- unlist(PropHydroAroScore)

pepTableb$PropHydroAroScore <- PropHydroAroScore

#peptide <- "YTDQISKYA"
#proportion <- PropHydroAro(peptide)
#print(proportion)
#
t_test <- t.test(PropHydroAroScore ~ tgt, data = pepTableb)
p_value <- t_test$p.value
ggplot(pepTableb, aes(x = as.factor(tgt), y = PrimeScore)) +
    geom_boxplot() +
    labs(x = "Target (tgt)", y = "PropHydroAroScore", title = "Boxplot of PropHydroAroScore by Target") +
    theme_minimal()+
    annotate("text", x = 1.5, y = max(pepTableb$PropHydroAroScore), 
             label = paste("p-value =", signif(p_value, digits = 7)), size = 5)

### Prime score ----

# cd /home/mxc2982/sccc/TOOLS/MixMHCpred3.0/MixMHCpred
# conda activate PRIME

mhc_types <- unique(pepTableb$mhc)
for (mhct in mhc_types){
    filtered_data <- pepTableb %>% filter(mhc == mhct) 
    if (!dir.exists("temp")) {
        dir.create("temp")
    }
    fname <- paste0("temp/",mhct,"_peps.txt")
    write.table(filtered_data$pep,file=fname,quote = FALSE,row.names = FALSE,col.names = FALSE)
    ofname <- paste0("temp/",mhct,"_peps_primeout.txt")
    prime_cmd <- paste0("./PRIME/PRIME/PRIME -i ",fname, " -o ", ofname, " -a ", mhct, " -mix ./MixMHCpred/MixMHCpred/MixMHCpred" )
    
    #command <- paste("/home/ceccarelli/miniforge3/condabin/conda activate PRIME ",  ";", prime_cmd)
    #result <- system2("bash", args = c("-c", command), stdout = TRUE, stderr = TRUE)
    cat(prime_cmd, file = "commands.sh", sep = "\n",append = TRUE)
}
#read the PRIME output
pepTableb$PrimeScore <- NA
for (mhct in mhc_types){
    ofname <- paste0("temp/",mhct,"_peps_primeout.txt")
    primeout<-read.csv(ofname,sep="\t",skip = 11)
    primeout$Score_bestAllele
    filtered_data <- pepTableb %>% filter(mhc == mhct) 
    if(sum(match(filtered_data$pep,primeout$Peptide)!=seq(1:nrow(filtered_data)))!=0){
        print ("Something wrong")
    }
    positions <- which(pepTableb$mhc == mhct)
    pepTableb$PrimeScore[positions]<- primeout$Score_bestAllele
}
save(pepTableb,file="PEPTIDE_IM.RData")
library(ggplot2)
library(ggpubr)
t_test <- t.test(PrimeScore ~ tgt, data = pepTableb)

# Extract the p-value
p_value <- t_test$p.value
ggplot(pepTableb, aes(x = as.factor(tgt), y = PrimeScore)) +
    geom_boxplot() +
    labs(x = "Target (tgt)", y = "PrimeScore", title = "Boxplot of PrimeScore by Target") +
    theme_minimal()+
    annotate("text", x = 1.5, y = max(pepTableb$PrimeScore), 
             label = paste("p-value =", signif(p_value, digits = 7)), size = 5)

#then go in the terminal
#conda activate PRIME
#bash commands
### NetMHCpan ----
convertmhc <- function(vec){
    tmp <- gsub("HLA-","", vec )
    tmp <-  gsub(":","",tmp,fixed =T )
    tmp <- gsub("*","",tmp ,fixed =T )
}
library(data.table) 
netmhcres<- read.csv("~/home3/neoantigens/NetMHCpanmerged.csv")
netmhcres$mhc <- convertmhc(netmhcres$MHC)
netmhcres<- as.data.table(netmhcres)
netmhcres$mhc_pep <- paste0(netmhcres$mhc,"_",netmhcres$Peptide)
mergedTable <- inner_join(pepTableb,netmhcres,by="mhc_pep")

sum(mergedTable$mhc_pep!=pepTableb$mhc_pep)
pepTableb$Rank_EL <- mergedTable$Rank_EL
pepTableb$Rank_BA <- mergedTable$Rank_BA
save(pepTableb,file="PEPTIDE_IM.RData")


### 

### peptide properties ----
# The following  were calculated  using the ProteinAnalysis  module from
# BioPython: molecular  weight (mw);  molar extinction  coefficient; the
# relative frequency  of F, W, and  Y amino acids or  aromaticity (Aro);
# instability       index      (Inst);       and      the       relative
# frequencyo f V,I,Y,F,W,and L aminoacids or helix  (PropHydroAro)  (31).  The
# isoelectric  point  was  calculated  using EMBOSS  with  the  Peptides
# package (34) in R.Themean  hydrophobicity scale (33) and  the proportion
# of different physicochemical classes  of amino acids (small, aromatic,
# acidic, and basic) were calculated for the non-anchor subsequence.




######## COmpute peptide properties:
library(Biostrings)
library(Peptides)
library(seqinr)
library(stringr)
library(data.table)

# Define functions to calculate various parameters
calculate_aromaticity <- function(sequence) {
    amino_acids <- strsplit(sequence, NULL)[[1]]
    aromatic_count <- sum(amino_acids %in% c("F", "W", "Y"))
    aromaticity <- aromatic_count / length(amino_acids)
    return(aromaticity)
}

calculate_extinction_coefficient <- function(sequence) {
    n_W <- str_count(sequence, "W")
    n_Y <- str_count(sequence, "Y")
    n_C <- str_count(sequence, "C") / 2  # Each C forms half a cystine
    extinction_coefficient <- (n_W * 5500) + (n_Y * 1490) + (n_C * 125)
    return(extinction_coefficient)
}

get_non_anchor_subsequence <- function(sequence, anchor_positions) {
    seq_vec <- strsplit(sequence, NULL)[[1]]
    non_anchor_seq <- seq_vec[-anchor_positions]
    non_anchor_seq <- paste(non_anchor_seq, collapse = "")
    return(non_anchor_seq)
}

calculate_aa_composition <- function(sequence) {
    amino_acids <- strsplit(sequence, NULL)[[1]]
    total_aa <- length(amino_acids)
    
    small <- sum(amino_acids %in% c("A", "G", "C", "S", "P", "D", "N", "T")) / total_aa
    aromatic <- sum(amino_acids %in% c("F", "W", "Y")) / total_aa
    acidic <- sum(amino_acids %in% c("D", "E")) / total_aa
    basic <- sum(amino_acids %in% c("R", "H", "K")) / total_aa
    
    return(c(small = small, aromatic = aromatic, acidic = acidic, basic = basic))
}

calculate_proportion_aromatic <- function(sequence) {
    amino_acids <- strsplit(sequence, NULL)[[1]]
    total_aa <- length(amino_acids)
    aromatic_count <- sum(amino_acids %in% c("F", "W", "Y"))
    proportion_aromatic <- aromatic_count / total_aa
    return(proportion_aromatic)
}

calculate_proportion_cysteine <- function(sequence) {
    amino_acids <- strsplit(sequence, NULL)[[1]]
    total_aa <- length(amino_acids)
    cysteine_count <- sum(amino_acids %in% c("C"))
    proportion_cysteine <- cysteine_count / total_aa
    return(proportion_cysteine)
}

calculate_proportion_small <- function(sequence) {
    amino_acids <- strsplit(sequence, NULL)[[1]]
    total_aa <- length(amino_acids)
    small_count <- sum(amino_acids %in% c("A", "G", "C", "S", "P", "D", "N", "T"))
    proportion_small <- small_count / total_aa
    return(proportion_small)
}

calculate_proportion_acidic <- function(sequence) {
    amino_acids <- strsplit(sequence, NULL)[[1]]
    total_aa <- length(amino_acids)
    acidic_count <- sum(amino_acids %in% c("D", "E"))
    proportion_acidic <- acidic_count / total_aa
    return(proportion_acidic)
}

calculate_proportion_basic <- function(sequence) {
    amino_acids <- strsplit(sequence, NULL)[[1]]
    total_aa <- length(amino_acids)
    basic_count <- sum(amino_acids %in% c("R", "H", "K"))
    proportion_basic <- basic_count / total_aa
    return(proportion_basic)
}

# Compute parameters
anchor_positions <- c(2, 9)
pepTableb$HydroCore <-unlist(mclapply(pepTableb$pep, function(seq) hydrophobicity(get_non_anchor_subsequence(seq, anchor_positions)),
                                      mc.cores = 64)) 

pepTableb$HydroAll <-unlist(mclapply(pepTableb$pep, Peptides::hydrophobicity,
                                     mc.cores = 64)) 
pepTableb$Aro <-unlist(mclapply(pepTableb$pep, calculate_aromaticity,
                                mc.cores = 64)) 
pepTableb$PropAro <-unlist(mclapply(pepTableb$pep, function(seq) calculate_proportion_aromatic(get_non_anchor_subsequence(seq, anchor_positions)),
                                    mc.cores = 64)) 
pepTableb$CysRed <- unlist(mclapply(pepTableb$pep, calculate_proportion_cysteine,
                                    mc.cores = 64)) 
pepTableb$PropSmall <- unlist(mclapply(pepTableb$pep, function(seq) calculate_proportion_small(get_non_anchor_subsequence(seq, anchor_positions)),
                                       mc.cores = 64)) 
pepTableb$PropAcidic <-unlist(mclapply(pepTableb$pep, function(seq) calculate_proportion_acidic(get_non_anchor_subsequence(seq, anchor_positions)),
                                       mc.cores = 64)) 
pepTableb$Inst <- unlist(mclapply(pepTableb$pep, instaIndex,
                                  mc.cores = 64)) 
pepTableb$PropBasic <- unlist(mclapply(pepTableb$pep, function(seq) calculate_proportion_basic(get_non_anchor_subsequence(seq, anchor_positions)),
                                       mc.cores = 64)) 
pepTableb$pI <- unlist(mclapply(pepTableb$pep, pI,
                                mc.cores = 64)) 
pepTableb$mw <- unlist(mclapply(pepTableb$pep, Peptides::mw,
                                mc.cores = 64)) 
save(pepTableb,file="PEPTIDE_IM.RData")
library(writexl)
write_xlsx(pepTableb,"PEPTIDE_IM_20240613.xlsx")
