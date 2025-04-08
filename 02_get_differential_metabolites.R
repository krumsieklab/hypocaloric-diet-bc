library(glmnet)
library(caret)
library(nlme)
library(parallel)
library(ggplot2)
library(ggpubr)
library(openxlsx)

##### Data input and output file names #####
manuscript_folder <- "/udd/nhast/keto-metabolomics/SchweickartAnnalise_HypocaloricDiets_022025"

source(paste0(manuscript_folder, "/Scripts/internal_functions.R"))


KD_knowns_file <- 
  paste0(manuscript_folder, 
         "/Processed_Data_and_Results/KD_known_metabolites_processed_batch_corrected.xlsx")
LFD_knowns_file <- 
  paste0(manuscript_folder, 
  "/Processed_Data_and_Results/LFD_known_metabolites_processed.xlsx")

# outputs
known_mets_analyzed <- 
  paste0(manuscript_folder, 
         '/Processed_Data_and_Results/known_mets_association_results.xlsx')
sig_diff_mets <- 
  paste0(manuscript_folder, 
         "/Processed_Data_and_Results/significant_diet_associated_mets.xlsx")


##### DIFFERENTIAL TESTS #####


#' Single diet paired t-test wrapper
#' wrapper function for getting adjusted p-values of all metabolites 
#' when testing abundance change between first and last week of a diet
#'
#'
#' @param D SummarizedExperiment of metabolites for samples of a signle diet
#' @param correction_method which p correction method to use, default "fdr"
#' @param cutoff whihc p-value cutoff to use, default 0.05
#' @param cores number of cores to run this function over, default 4
#'
#' @return dataframe with significantly different metabolites
#'
#'
paired_test_wrapper <- function(
    D, 
    correction_method = "fdr", 
    cutoff = 0.05,
    cores = 4 
){
  
  # Extract subject data out of Summarized Experiment
  data=t(assay(D))
  time=colData(D)$Week
  sampleID=colData(D)$Subject
  group=colData(D)$Diet
  Universal_ID = colnames(data)
  HMDB_ID = rowData(D)$HMDB_ID
  
  # Get raw p-value comparing first and last week in specified diet per metabolite
  raw_ps <- mclapply(1:ncol(data), function(i){
    paired_test_ind(i, data, sampleID, time, group)
  },mc.cores= cores) %>% bind_rows()
  
  raw_ps$Metabolite = rowData(D)$Metabolite[raw_ps$index]
  raw_ps$p_adj = p.adjust(raw_ps$pval, method = correction_method)
  raw_ps$wilcox_p_adj <- p.adjust(raw_ps$wilcox_p, method = correction_method)
  raw_ps$Universal_ID = Universal_ID
  raw_ps$HMDB_ID = HMDB_ID
  raw_ps %>% 
    dplyr::select(index, Metabolite, p_adj, wilcox_p_adj, estimate, Universal_ID, HMDB_ID) 
}




#' Diet paired test for single metabolite
#' Gets raw p-values of a single metabolite  
#' when testing abundance change between first and last week of a diet
#'
#'
#' @param i index of single metabolite
#' @param data expression data of metabolites
#' @param sampleID sample IDs corresponding to metabolite abundance values
#' @param time time corresponding to metabolite abundance values
#' @param group diet corresponding to metabolite abundance values
#'
#' @return data.frame with index of metabolite, its raw p-value, and the direction of change
#'
#'
paired_test_ind<-function(i, data, sampleID, time, group){
  
  # Create dataframe of subject data
  dataM <- data.frame(Rep = sampleID, Group = as.factor(group),
                      time = as.numeric(time), Expr = as.numeric(data[,i]))
  
  # For each subject
  test <- lapply(unique(dataM$Rep), function(rep){
    # Get the starting abundance
    start_val <- dataM %>% 
      dplyr::filter(Rep == rep) %>% 
      dplyr::filter(time == min(time)) %>% 
      dplyr::select(Rep, Group, Expr)
    
    # Get the ending abundace
    end_val <- dataM %>% 
      dplyr::filter(Rep == rep) %>% 
      dplyr::filter(time == max(time)) %>% 
      dplyr::select(Rep, Group, Expr)
    df <- bind_rows(start_val, end_val)
    df$timepoint = c("before","after")
    df
  }) %>% bind_rows()
  
  # Perform paired t-test
  t_test  <- t.test(test$Expr[test$timepoint=="after"], 
                    test$Expr[test$timepoint=="before"], paired=T)
  
  wilcox_test  <- wilcox.test(test$Expr[test$timepoint=="after"], 
                              test$Expr[test$timepoint=="before"], paired=T)
  
  
  # Return dataframe with metabolite index, raw p-value and direction of change
  data.frame(index = i, pval = t_test$p.value, estimate = t_test$estimate, wilcox_p = wilcox_test$p.value)
}




##### Load Data #####
KD_SE <- mt_load_se_xls(file = KD_knowns_file)

LFD_SE <- mt_load_se_xls(file = LFD_knowns_file)

rownames(KD_SE) <- rowData(KD_SE)$HMDB_ID

rownames(LFD_SE) <- rowData(LFD_SE)$HMDB_ID

# Annotate control platforms into first week (pre diet) and last week (post diet)
colData(KD_SE)$pre_post <- assign_pre_post(KD_SE)
colData(LFD_SE)$pre_post <- assign_pre_post(LFD_SE)

KD_fc_data <- get_diet_fc(KD_SE)
LFD_fc_data <- get_diet_fc(LFD_SE)

#### Get p-values for metabolite changes pre- to post- diet for both groups ####

# Get metabolite adjusted p-values of differential analysis 
# between the first and last time points in KD subjects
paired_ps_KD <- paired_test_wrapper(KD_SE)

# Get metabolite adjusted p-values of differential analysis 
# between the first and last time points in SD subjects
paired_ps_lfd <- paired_test_wrapper(LFD_SE)

# Change names
paired_ps_KD %<>% dplyr::rename(KD_change_p_adj = p_adj) %>% 
  dplyr::select(-c(index, Universal_ID))

paired_ps_lfd %<>% dplyr::rename(LFD_change_p_adj = p_adj) %>% 
  dplyr::select(-c(index, Universal_ID))

# Merge data
all_changes <- merge(paired_ps_KD, paired_ps_lfd,by="HMDB_ID") %>% 
  dplyr::rename(Metabolite = Metabolite.x) %>% 
  select(-Metabolite.y)

# Get fold change information with confidence intervals for both diet groups

KD_ci_data <- data.frame(t(sapply(KD_fc_data, Rmisc::CI))) %>% 
  dplyr::rename(KD_lower = lower) %>% 
  dplyr::rename(KD_upper = upper) %>% 
  dplyr::rename(KD_mean = mean)

KD_ci_data$HMDB_ID <- rownames(assay(KD_SE))


lfd_ci_data <- data.frame(t(sapply(LFD_fc_data, Rmisc::CI)))%>% 
  dplyr::rename(LFD_lower = lower) %>% 
  dplyr::rename(LFD_upper = upper) %>% 
  dplyr::rename(LFD_mean = mean)

lfd_ci_data$HMDB_ID <- rownames(assay(LFD_SE))


# Combine p-value and fold change information
all_changes %<>% merge(KD_ci_data, by = "HMDB_ID") %>% 
  merge(lfd_ci_data, by = "HMDB_ID")

# Annotate metabolites into groups:
# --> Caloric Restriction Metabolite significantly changes in the same direction for both diets
# --> LFD Metabolite significantly changes in LFD but not KD
# --> Keto Metabolite significantly changes in KD but not LFD
# --> Diet Divergent Metabolite significantly changes in both diets but in different directions

all_changes$annotation <- unlist(lapply(1:nrow(all_changes), function(i) {
  if (all_changes$KD_change_p_adj[i] > 0.05 & all_changes$LFD_change_p_adj[i] < 0.05){
    "LFD Metabolite"
  }
  
  # If KD is significant but LFD isn't, this is a keto metabolite
  else if (all_changes$KD_change_p_adj[i] < 0.05 & all_changes$LFD_change_p_adj[i] > 0.05){
    "Keto Metabolite"
  }
  
  # If both LFD and KD are significant in the same direction, 
  # Labelled weightloss metabolite, otherwise keto metabolite (more important annotation for NHS)
  else if (all_changes$KD_change_p_adj[i] < 0.05 & all_changes$LFD_change_p_adj[i] < 0.05){
    if (sign(all_changes$estimate.x[i]) == sign(all_changes$estimate.y[i])){
      'Caloric Restriction Metabolite'
    }
    else{'Diet Divergent Metabolite'}
  }
  else{NA}
}))


# Load in manual annotations for all metabolites
all_changes %<>% merge(data.frame(rowData(KD_SE)), by ="HMDB_ID", all.x=T) %>% 
  distinct()


# Parse out lipid species to get chain length and saturation
fa_parsed <- lapply(all_changes$Fat_type, function(fa_name){
  fa_type = strsplit(fa_name, "\\(")[[1]][1]
  length = gsub("(?<=\\()[^()]*(?=\\:)(*SKIP)(*F)|.", "", fa_name, perl=T) %>% as.numeric()
  saturation = gsub("(?<=\\:)[^()]*(?=\\))(*SKIP)(*F)|.", "", fa_name, perl=T) %>% as.numeric()
  data.frame(Fat_type = fa_name, fa_type, length, saturation)
}) %>% bind_rows() %>% na.omit()

fa_types <- unique(fa_parsed$fa_type)

all_changes <- merge(fa_parsed, all_changes, by = "Fat_type", all.y=T) %>% distinct()


# write out all metabolites analyzed
out <- 'known_mets_analyzed'
wb <- openxlsx::createWorkbook()
# creat worksheet
openxlsx::addWorksheet(wb,sprintf('%s', out))
# write data
openxlsx::writeData(wb, sprintf('%s', out), all_changes, rowNames = F, colNames = T)
# create and add a style to the column headers
headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
# style for body
bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
# apply style
addStyle(wb, sheet = sprintf('%s', out), bodyStyle, rows = 1:nrow(all_changes), cols = 1:ncol(all_changes), gridExpand = TRUE)
addStyle(wb, sheet = sprintf('%s', out), headerStyle, rows = 1, cols = 1:ncol(all_changes), gridExpand = TRUE)
# write out
openxlsx::saveWorkbook (wb, file=known_mets_analyzed, overwrite=TRUE)


