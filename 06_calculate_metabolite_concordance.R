# Calculation of concordance values for each metabolite in diet scores
# And comparison of concordance values with fold changes
load(paste0(manuscript_folder, 
            "/Processed_Data_and_Results/NHS_metadata_scores_df_maplet.Rdata"))
load(paste0(manuscript_folder, 
            "/Processed_Data_and_Results/NHS2_ExprSet_processed_maplet.RData"))

# inputs
betas_file <- paste0(manuscript_folder, 
  "/Processed_Data_and_Results/filtered_betas.xlsx")


## KD

# Get fold changes and metabolite names from NHS data
kd_annotation_col =  read_excel(betas_file, 
                             sheet="NHS2_betas_all_keto") %>% 
  arrange(desc(beta)) %>% 
  data.frame() %>% 
  dplyr::rename(Metabolite_FC= beta) %>% 
  merge(fData(NHS2_ExprSet_processed), by.x = "HMDB_ID", by.y = "row.names") %>% 
  select(HMDB_ID, Metabolite_FC, biochemical_name)

# Get processed metabolite data from NHS2
kd_met_dat <- t(exprs(NHS2_ExprSet_processed)) %>% as.data.frame() %>% 
  select(kd_annotation_col$HMDB_ID)

rownames(kd_annotation_col) <- colnames(kd_met_dat)
kd_met_dat$id <- substr(rownames(kd_met_dat),1,6)

kd_scores = NHS2_all_keto_dat %>% select(id, beta.score, case)

# Merge metabolomite abundances with diet scores by subject
kd_merged = merge(kd_met_dat, kd_scores, by = "id") %>% 
  arrange(desc(beta.score)) %>% 
  dplyr::rename(KD_Score = "beta.score")

annotation_row = data.frame(KD_Score = kd_merged$KD_Score)

rownames(annotation_row) = rownames(kd_merged) = kd_merged$id


kd_full_data <- as.matrix(kd_merged %>% select(-c(id, case, KD_Score)))

# Perform pearson correlation between metabolite abundance and 
# Diet score to calculate concordance
kd_order <- data.frame(met = colnames(kd_full_data), 
                       order = lapply(1:ncol(kd_full_data), function(i){
                         cor.test(kd_full_data[,i], annotation_row$KD_Score)$estimate
                       }) %>% unlist()) %>% arrange(desc(order))


kd_merged <- merge(kd_annotation_col, kd_order, by.x = "HMDB_ID", by.y = "met") %>% 
  dplyr::rename(NHSII_concordance = "order")


## LFD

lfd_annotation_col =  read_excel(betas_file, 
                                sheet="NHS2_betas_all_lfd") %>% 
  arrange(desc(beta)) %>% 
  data.frame() %>% 
  dplyr::rename(Metabolite_FC= beta) %>% 
  merge(fData(NHS2_ExprSet_processed), by.x = "HMDB_ID", by.y = "row.names") %>% 
  select(HMDB_ID, Metabolite_FC, biochemical_name)

# Get processed metabolite data from NHS2
lfd_met_dat <- t(exprs(NHS2_ExprSet_processed)) %>% as.data.frame() %>% 
  select(lfd_annotation_col$HMDB_ID)

rownames(lfd_annotation_col) <- colnames(lfd_met_dat)
lfd_met_dat$id <- substr(rownames(lfd_met_dat),1,6)

lfd_scores = NHS2_all_lfd_dat %>% select(id, beta.score, case)

lfd_merged = merge(lfd_met_dat, lfd_scores, by = "id") %>% 
  arrange(desc(beta.score)) %>% 
  dplyr::rename(LFD_Score = "beta.score")

annotation_row = data.frame(LFD_Score = lfd_merged$LFD_Score)

rownames(annotation_row) = rownames(lfd_merged) = lfd_merged$id

lfd_full_data <- as.matrix(lfd_merged %>% select(-c(id, case, LFD_Score)))

# Perform pearson correlation between metabolite abundance and 
# Diet score to calculate concordance
lfd_order <- data.frame(met = colnames(lfd_full_data), 
                        order = lapply(1:ncol(lfd_full_data), function(i){
                          cor.test(lfd_full_data[,i], annotation_row$LFD_Score)$estimate
                        }) %>% unlist()) %>% arrange(desc(order))


lfd_merged <- merge(lfd_annotation_col, lfd_order, by.x = "HMDB_ID", by.y = "met") %>% 
  dplyr::rename(NHSII_concordance= "order") 


save(kd_merged, lfd_merged, 
     file = paste0(manuscript_folder, 
                   "/Processed_Data_and_Results/Concordance_information.RData"))
