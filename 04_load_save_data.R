
# Functions for loading data and running association analysis
library(sas7bdat)
library(tidyverse)
library(readxl)
library(survival) 
library(mgcv)
library(openxlsx)
library(data.table)
library(Biobase)
library(dplyr)
library(glue)
library(gt)

manuscript_folder <- "/udd/nhast/keto-metabolomics/SchweickartAnnalise_HypocaloricDiets_022025"

# inputs
betas_file <- 
  paste0(manuscript_folder, 
         "/Processed_Data_and_Results/filtered_betas.xlsx")

nhs2_pheno_file <- # Contains covariates/metadata of NHS2 subjects
  "/udd/recpe/Metabolomics/KetoMetabolomics/nhs2final.sas7bdat"

nhs2_meno_er_file <- # Contains ER and menopausal status of cancer in NHS2
  "/udd/n2bch/request/Metab_2017/NHSII/data_mtbl_nhs2_18DEC18.csv"

# Load in metabolomics data

load(paste0(manuscript_folder, 
            "/Processed_Data_and_Results/NHS2_ExprSet_processed_maplet.RData"))

# Loads metabolomics data, covariates, and cancer case/control information
get_data_for_assoc <- function(
    ExprSet,
    beta_file,
    beta_sheet,
    pheno_data,
    meno_er_data){
  # Read in metabolomics data
  mData <- exprs(ExprSet)
  fData <- fData(ExprSet)
  pData <- pData(ExprSet)
  
  # Read in phenotype information, align with metabolomics pData
  pData$id6 <- substr(as.character(pData$id), 1, 6)
  pData$id7 <- pData$id
  pData <- subset(pData, select=-c(id))
  phenofinal <- merge(pheno_data,pData,by.x='id',by.y='id6')
  
  # Fill in missing metadata with median values
  phenofinal$bmi18[is.na(phenofinal$bmi18)]<-
    median(phenofinal$bmi18,na.rm=TRUE)
  phenofinal$weightchange[is.na(phenofinal$weightchange)]<-
    median(phenofinal$weightchange,na.rm=TRUE)
  phenofinal$act[is.na(phenofinal$act)]<-
    median(phenofinal$act,na.rm=TRUE)
  phenofinal$ahei_noal[is.na(phenofinal$ahei_noal)]<-
    median(phenofinal$ahei_noal,na.rm=TRUE)
  phenofinal$alco[is.na(phenofinal$alco)]<-
    median(phenofinal$alco,na.rm=TRUE)
  phenofinal$current.smk=ifelse(phenofinal$smoke=='current','1','0')
  phenofinal$menopause=ifelse(phenofinal$menopmh=='pre','pre','not_pre')
  phenofinal$fasting=ifelse(phenofinal$fast=='fasting','1','0')
  phenofinal$current.smk[is.na(phenofinal$current.smk)]<-'0'
  
  # Correct out batch, save residuals
  expr <- mData[,match(as.character(phenofinal$id7), colnames(mData))]
  expr <- as.matrix(expr)
  
  t.expr=t(expr)
  
  if(length(unique(phenofinal$labcode))>1){
    resExpr0=t.expr*NaN
    for (i in 1: dim(t.expr)[2]){
      residual=resid(lm(t.expr[,i] ~  as.factor(labcode), 
                        data=phenofinal, na.action=na.exclude)) 
      resExpr0[,i]=residual
      #print(i)
    }
    resExpr0.t=t(resExpr0)
  }
  else {resExpr0.t <- expr}
  
  
  # Load in metabolite betas
  betas <- read_excel(beta_file, sheet=beta_sheet)
  
  # Select metabolomics information for metabolites in beta sheet
  expression.resid=as.data.frame(resExpr0.t)
  expr.selected <- expression.resid[match(as.character(betas$HMDB_ID), 
                                          rownames(expression.resid)),]
  
  expr.selected <- expr.selected[complete.cases(expr.selected),]
  selected <- betas[match(rownames(expr.selected), as.character(betas$HMDB_ID)),]
  
  # Calculate diet score (dot product of betas and metabolite values)
  phenofinal.score <- phenofinal
  phenofinal.score$beta.score <- t(expr.selected) %*% selected$beta %>% .[,1]
  
  # Extract quartiles of scores from controls
  phenofinal.score$case=ifelse(phenofinal.score$caco=='case',1,0)
  
  # Add pre- or post- menopausal breast cancer
  # This removes samples (n=2) with missing ER and/or menopausal status at dx 
  meno_er_dat <- meno_er_data %>% 
    select(id, menopmh_at_dx, menopmh_dx, menopmh_bld, erpos, menarc, nulli, par_num, invas)
  phenofinal.score = merge(phenofinal.score, meno_er_dat, by = "id", all.X = T) 
  phenofinal.score
}

#### Make data with scores #####


##### NHS2 #####

NHS2_pheno_data <-read.sas7bdat(nhs2_pheno_file)
NHS2_meno_er_data <- read.csv(nhs2_meno_er_file)

# Keto analysis
NHS2_all_keto_dat <- get_data_for_assoc(
  ExprSet = NHS2_ExprSet_processed,
  beta_file = betas_file,
  beta_sheet = "NHS2_betas_all_keto",
  pheno_data = NHS2_pheno_data,
  meno_er_data = NHS2_meno_er_data)


# LFD analysis
NHS2_all_lfd_dat <- get_data_for_assoc(
  ExprSet = NHS2_ExprSet_processed,
  beta_file = betas_file,
  beta_sheet = "NHS2_betas_all_lfd",
  pheno_data = NHS2_pheno_data,
  meno_er_data = NHS2_meno_er_data)


save(NHS2_all_keto_dat,
     NHS2_all_lfd_dat,
     file = paste0(manuscript_folder, 
                   "/Processed_Data_and_Results/NHS_metadata_scores_df_maplet.Rdata"))
