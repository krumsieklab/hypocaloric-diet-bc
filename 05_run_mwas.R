# inputs
library(Biobase)
library(dplyr)

manuscript_folder <- 
  "/udd/nhast/keto-metabolomics/SchweickartAnnalise_HypocaloricDiets_022025"


betas_first <- paste0(manuscript_folder, 
                      "/Processed_Data_and_Results/filtered_betas.xlsx")
nhs2_pheno_file <- "/udd/recpe/Metabolomics/KetoMetabolomics/nhs2final.sas7bdat"
load(paste0(manuscript_folder, 
            "/Processed_Data_and_Results/NHS2_ExprSet_processed_maplet.RData"))

results_folder <- paste0(manuscript_folder, "/Processed_Data_and_Results/")

get_diet_mwas_res <- function(
    ExprSet,
    pheno_file,
    hmdbs_to_test,
    pre_post = "both"){
  
  # Read in metabolomics data
  mData <- exprs(ExprSet)
  fData <- fData(ExprSet)
  pData <- pData(ExprSet)
  
  # Read in phenotype information
  pheno_data <- read.sas7bdat(pheno_file)
  pData$id6 <- substr(as.character(pData$id), 1, 6)
  pData$id7 <- pData$id
  pData <- subset(pData, select=-c(id))
  phenofinal <- merge(pheno_data,pData,by.x='id',by.y='id6')
  
  expr <- mData[,match(as.character(phenofinal$id7), colnames(mData))]
  expr <- as.matrix(expr)
  
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
  
  # Subselect for desired metabolites to test
  expr.selected <- expr[match(as.character(hmdbs_to_test), rownames(expr)),]
  expr.selected <- expr.selected[complete.cases(expr.selected),]

  # Collect results
  expr.selected=as.matrix(expr.selected)
  m <- matrix(NA, ncol = 5, nrow = nrow(expr.selected))
  m = data.frame(m)
  rownames(m)=rownames(expr.selected)
  
  for(ii in 1:nrow(m)){
    
    results <- clogit(case ~ expr.selected[ii,] + ageyr + as.factor(menopause) + 
                        as.factor(fast) + act + as.factor(current.smk) + 
                        ahei_noal + alco + blddate + menage + as.factor(afbpar) + 
                        familyhx + bbddx + bmi18 + as.factor(changecat) + 
                        strata(matchid), phenofinal)
    temp=summary(results)
    m[ii,] <- temp$coefficients[1,]
  }
  
  m$id <- rownames(m)
  colnames(m)=c('coef','OR','se','z.stats','p.value','hmdb.id')
  results<-data.table(m)
  
  annotation = fData[c('metabolite_name','hmdb_id','biochemical_name',
                       'class_metabolon','sub_class_metabolon')]
  
  anno_matched <- annotation[match(as.character(results$hmdb.id), 
                                   rownames(annotation)),]
  
  res=cbind(results,anno_matched)
  
  # calculate FDR and Bonferroni vars
  res[, FDR := p.adjust(p.value, method = "BH")]
  res[, Bonferroni := p.adjust(p.value, method = "bonferroni")]
  
  res
}


#### Get MWAS Results

diet_hmdbs_nhs2 <-
  rbind(read.xlsx(betas_first, sheet ='NHS2_betas_all_lfd'),
        read.xlsx(betas_first, sheet ='NHS2_betas_all_keto')) %>% 
  .$HMDB_ID %>% unique()

nhs2_mwas_res <- 
  get_diet_mwas_res(NHS2_ExprSet_processed,
                    pheno_file = nhs2_pheno_file,
                    diet_hmdbs_nhs2)   

nhs2_mwas_file <- "NHS2.MWAS.withBMI.csv"

write.csv(nhs2_mwas_res, 
          file = paste0(results_folder, nhs2_mwas_file), row.names = T)
