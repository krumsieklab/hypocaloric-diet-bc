## Load in NHS data from both nhs1 and nhs2 cohorts
## Impute half min, probit transform, remove samples and features
## With too much missing data



# Input file names
met_diet_assoc_file <- 
  paste0(manuscript_folder, 
         "/Processed_Data_and_Results/known_mets_association_results.xlsx")

# Output file names
filtered_betas_file <- 
  paste0(manuscript_folder, "/Processed_Data_and_Results/filtered_betas.xlsx")

###### MAPLET PREPROCESSING ######
# Load raw data (no transformation, no imputation)

all.raw.NHS2 = merge_metab_data(
  merge_type = "union", combine_cohorts = T, collection_to_use = "first",
  cohorts = c("nhs2"), endpoints = c("breast"), keep_failed_pm_metabolites = F,
  impute_cutoff = 0, archive_date = "2024-10-09", 
  transformation = "transform_none")

#Extract pData,fData and mData from raw data

data_nhs2<-all.raw.NHS2$expr_set$all_cohorts

nhs2_pData<-pData(data_nhs2)
nhs2_fData<-fData(data_nhs2)
nhs2_mData<-exprs(data_nhs2)

## Make sure the fData matches the rows and pData matches the columns
## Of mData (both should return TRUE)

all(rownames(nhs2_mData) == rownames(nhs2_fData))
all(colnames(nhs2_mData) == rownames(nhs2_pData))


#' Preprocess broad data
#'
#' Perform classic metabolomics preprocessing steps including
#' removing missing data, filtering out QC with high CoV
#' quotient normalize, log transform and scale, imput missing values
#' 
#' @param D SummarizedExperiment containing raw Broad metabolomics data

#'
#' @return Preprocessed SummarizedExperiment
#'
#' @noRd
preprocess_broad_data <- function(D){
  
  processed_dat <- D %>% 
    
    # Filter out metabolites with over 25% missingness
    mt_pre_filter_missingness(feat_max=0.25) %>%
    
    #Filter out samples with over 25% missingness
    mt_pre_filter_missingness(samp_max=0.25) %>%
    
    #Quotient Normalize
    mt_pre_norm_quot() %>%
  
    #Log transform
    mt_pre_trans_log() %>%
    
    # metabolic outlier detection followed by imputation
    mt_pre_outlier_to_na() %>% 
    
    # Filter out metabolites with over 25% missingness
    mt_pre_filter_missingness(feat_max=0.25) %>%
    
    #Impute NAs
    mt_pre_impute_knn() %>% 
    
    # Scale data
    mt_pre_trans_scale()
  
  
  processed_dat
}


NHS2_SE <- SummarizedExperiment(assays = nhs2_mData,
                                colData = nhs2_pData,
                                rowData = nhs2_fData) %>% 
  preprocess_broad_data()

# Save NHS2 processed data as ExpressionSet

NHS2_ExprSet_processed <- 
  ExpressionSet(assayData = assay(NHS2_SE), 
                phenoData = AnnotatedDataFrame(data.frame(colData(NHS2_SE))),
                featureData = AnnotatedDataFrame(data.frame(rowData(NHS2_SE))))
save(NHS2_ExprSet_processed,
     file = "Processed_Data_and_Results/NHS2_ExprSet_processed_maplet.RData")

###### Get metabolites and betas from diet study results#####


# Known NHS metabolites
NHS2_knowns <- rownames(fData(NHS2_ExprSet_processed))

sig_diff_mets <- read.xlsx(met_diet_assoc_file)

NHS2_betas_all_keto <- sig_diff_mets %>% 
  filter(KD_change_p_adj < 0.05) %>% 
  select(HMDB_ID, estimate.x) %>% 
  dplyr::rename(beta = estimate.x) %>% 
  filter(HMDB_ID %in% NHS2_knowns)


NHS2_betas_all_lfd <- sig_diff_mets %>% 
  filter(LFD_change_p_adj < 0.05) %>% 
  select(HMDB_ID, estimate.y) %>% 
  dplyr::rename(beta = estimate.y) %>% 
  filter(HMDB_ID %in% NHS2_knowns)


wb <- openxlsx::createWorkbook()

# create worksheet
openxlsx::addWorksheet(wb,"NHS2_betas_all_keto")
# write data
openxlsx::writeData(wb, "NHS2_betas_all_keto", NHS2_betas_all_keto, 
                    rowNames = F, colNames = T)


# create worksheet
openxlsx::addWorksheet(wb,"NHS2_betas_all_lfd")
# write data
openxlsx::writeData(wb, "NHS2_betas_all_lfd", NHS2_betas_all_lfd, 
                    rowNames = F, colNames = T)

# write out
openxlsx::saveWorkbook (wb, file=filtered_betas_file, overwrite=TRUE)


