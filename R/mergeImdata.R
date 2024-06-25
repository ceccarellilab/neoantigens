#Merge multiple immunogenicity datasets

require(data.table)
setwd("~/neoantigens")


#data from #BigMHC Developed by the Karchin Lab (Nat Mach Intelligence 2023)
# Benjamin Alexander Albert, Deep neural networks predict class I major 
# histocompatibility complex epitope presentation
# and transfer learn neoepitope immunogenicity, Nature Machine Intelligence, 2023
# PMID: 37829001

author = "Alexander B. et al."
pmid = "37829001"
dataset="BigHMC"

# dsfolder <- "./Immunogenicity_datasets/datasets"

dsfolder <- "~/home3/neoantigens/BigMHC_Albert_etal/BigMHC Training and Evaluation Data/datasets"
im_test <- fread(file.path(dsfolder, "im_test.csv"))
im_train <- fread(file.path(dsfolder, "im_train.csv"))
im_val <- fread(file.path(dsfolder, "im_val.csv"))
iedb <- fread(file.path(dsfolder, "iedb.csv"))
table(iedb$tgt)
colnames(im_test)
colnames(im_train)
colnames(im_val)
im_test_s <- subset(im_test, select=c("mhc", "pep", "tgt"))
im_test_s$mhc <- convertmhc(im_test_s$mhc)
im_test_s$set <- "test"  

im_test_s$mhc_pep = paste0(im_test_s$mhc, "_", im_test_s$pep)
sum(duplicated(im_test_s$mhc_pep))

im_train <- subset(im_train, select=c("mhc","pep","tgt"))
im_train$mhc <- convertmhc(im_train$mhc)
im_train$set <- "train"
im_train$mhc_pep = paste0(im_train$mhc,"_",im_train$pep)
sum(duplicated(im_train$mhc_pep))

intersect(im_test_s$mhc_pep, im_train$mhc_pep)

im_val <- subset(im_val, select=c("mhc", "pep", "tgt"))
im_val$mhc <- convertmhc(im_val$mhc)
im_val$set <- "val"
im_val$mhc_pep = paste0(im_val$mhc, "_", im_val$pep)

intersect(im_test_s$mhc_pep, im_val$mhc_pep)
intersect(im_train$mhc_pep, im_val$mhc_pep)

iedb$mhc <- convertmhc(iedb$mhc)
iedb <- subset(iedb, select = c("mhc", "pep", "tgt"))
iedb$set <- "test_iedb"
iedb$mhc_pep <- paste0(iedb$mhc, "_", iedb$pep)
intersect(iedb$mhc_pep, im_train$mhc_pep)


PEPTIDE_IM <- rbind(im_train, im_test_s, im_val,iedb)
PEPTIDE_IM$source <- pmid
PEPTIDE_IM$author  <- author 
PEPTIDE_IM$dataset <- dataset

dim(PEPTIDE_IM)
table(PEPTIDE_IM$tgt, PEPTIDE_IM$set)
sum(duplicated(PEPTIDE_IM$mhc_pep))
dim(PEPTIDE_IM)


# Immunogenic peptides from Muller, Machine learning methods and harmonized 
# datasets improve immunogenic neoantigen prediction, Immunity 2023
# Datasets for mutations containing the feature values and immunogenicity 
# annotations

pmid=37816353
author = "Muller et al."
dsfolder <- "~/home3/neoantigens/Muller_etal_2023"
muts <- fread(file.path(dsfolder, "Mutation_data_org.txt"))
dim(muts)

# [1] 48306    59
colnames(muts)


# Datasets for neo-peptides containing the feature values and immunogenicity 
# annotations
peps <- fread(file.path(dsfolder, "Neopep_data_org.txt"))
dataset <- peps$dataset
dim(peps) # [1] 1787710      57
table(peps$Sample_Tissue)
# Adrenal Gland        Breast        Cervix Uteri 
# 3951                 73508         64070 
# Colon       Esophagus        Kidney 
# 303384          8077          3380 
# Liver          Lung      Pancreas 
# 283          197011         23580 
# Skin 
# 1110466 
table(peps$Cancer_Type)

# Breast Invasive Carcinoma         Cervical Cancer 
# 73508                             47985 
# Cholangiocarcinoma                  Colon Adenocarcinoma 
# 13896                                256329 
# Esophageal Carcinoma     Kidney Renal Clear Cell Carcinoma 
# 8077                      3380 
# Lung Adenocarcinoma          Lung Squamous Cell Carcinoma 
# 193354                                  3657 
# Pancreatic Adenoma and Adenocarcinoma                 Rectum Adenocarcinoma 
# 9684                                 47055 
# Sarcoma               Skin Cutaneous Melanoma 
# 4725                               1110466 
# Stomach Adenocarcinoma  Uterine Corpus Endometrial Carcinoma 
# 4234                                 11360 

names(peps)
sum(is.na(peps$response_type)) # 0
table(peps$response_type) 
# CD8   negative not_tested 
# 178     422907    1364625 

table(peps$train_test)
peps <- subset(peps, response_type %in% c("CD8","negative"))
dim(peps) 
table(peps$response_type)
#tmp<-subset(peps, select=c("mutant_best_alleles","mutant_seq","response_type"))

library(tidyr)
library(dplyr)


# Create a new data frame with unlisted alleles
pep1 <- peps %>%
    separate_rows(mutant_best_alleles, sep = ",") %>%
    mutate(response_type = ifelse(response_type == "negative"  | response_type == "negative ", 0, 1))

tmp<-subset(peps, select=c("mutant_best_alleles","mutant_seq","response_type","train_test","dataset"))
tmp$mhc_pep <- paste0(tmp$mutant_best_alleles,"_",tmp$mutant_seq)
#colnames(tmp[,1:5]) <- colnames(PEPTIDE_IM[1:5])
tmp$pmid = pmid
tmp$author = author

tmp <- tmp %>%
    select(mutant_best_alleles, mutant_seq, response_type, train_test, mhc_pep, pmid, author, dataset)

colnames(tmp) <- colnames(PEPTIDE_IM)
length(intersect(tmp$mhc_pep,PEPTIDE_IM$mhc_pep))
merged_df <- merge(tmp, PEPTIDE_IM, by = "mhc_pep")
#View(merged_df)
dim(merged_df)
table(tmp$dataset)
table(merged_df$dataset.x)
table(merged_df$tgt.x,merged_df$tgt.y)

#Omg!
# 
# 0   1
# 0 647  22
# 1   3  65
#I trust BigMHC (?)
tmp2 <- match(tmp$mhc_pep,PEPTIDE_IM$mhc_pep) 
sum(!is.na(tmp2))
tmp <- tmp[!!is.na(tmp2),]
dim(tmp)
table(tmp$tgt)
sum(duplicated((tmp$mhc_pep)))
#View(tmp[duplicated(tmp$mhc_pep),])
tmp <- tmp[!duplicated(tmp$mhc_pep),]
dim(tmp)
table(tmp$tgt)
PEPTIDE_IM <- rbind(PEPTIDE_IM,tmp)
dim(PEPTIDE_IM)
table(PEPTIDE_IM$tgt,PEPTIDE_IM$set)
sum(duplicated(PEPTIDE_IM$mhc_pep))


#data from IMPROVE

pmid="38633261"
author = "Borch et al."

dataImp <- fread("./IMPROVE_all_peptides.txt")
dim(dataImp)
#[1] 17520    98
dataImp$mhc <- convertmhc(dataImp$HLA_allele)
dataImp$Mut_peptide
dataImp$set <- "train"
dataImp$mhc_pep <- paste0(dataImp$mhc,"_",dataImp$Mut_peptide)
table(dataImp$response)
tmp <- dataImp %>%
    select(mhc,Mut_peptide,response,set,mhc_pep)
tmp$source <- pmid
tmp$author <- author
tmp$dataset <- "IMPROVE"
colnames(tmp) <- colnames(PEPTIDE_IM)
#check duplicates

sum(duplicated(tmp$mhc_pep))
#[1] 167
#check same tgt
duplicated_mhc_pep <- tmp[duplicated(tmp$mhc_pep), mhc_pep]
check_tgt <- tmp[mhc_pep %in% duplicated_mhc_pep, .(unique_tgt = unique(tgt)), by = mhc_pep]
#View(check_tgt)

tmp <- unique(tmp, by = "mhc_pep")


length(intersect(tmp$mhc_pep,PEPTIDE_IM$mhc_pep))
merged_df <- merge(tmp, PEPTIDE_IM, by = "mhc_pep")
#View(merged_df)
dim(merged_df)
table(tmp$dataset)
table(merged_df$dataset.x)
table(merged_df$tgt.x,merged_df$tgt.y)
#   
#0  1
#0 50  6
#1  1  3

tmp2 <- match(tmp$mhc_pep,PEPTIDE_IM$mhc_pep) 
sum(!is.na(tmp2))
tmp <- tmp[!!is.na(tmp2),]
dim(tmp)
table(tmp$tgt)

PEPTIDE_IM <- rbind(PEPTIDE_IM,tmp)
table(PEPTIDE_IM$tgt)
sum(duplicated(PEPTIDE_IM$mhc_pep))
save(PEPTIDE_IM, file="PEPTIDE_IM.RData")
library(openxlsx)

write.xlsx(PEPTIDE_IM, file = "PEPTIDE_IM.xlsx")
dim(PEPTIDE_IM)
