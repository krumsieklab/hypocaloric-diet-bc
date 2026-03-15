# Load and preprocess controlled trials data measured on
# Broad metabolomics platform 


#######################################
##### Data input and output files #####
#######################################

# Raw Broad Metabolomics files of the Control Ketogenic Diet Studies
# MSK is the endometrial cancer cohort
# OSU is the weight loss cohort

MSK_raw <- paste0(manuscript_folder, "/Input_Data/MSK.xlsx")
OSU_raw <- paste0(manuscript_folder, "/Input_Data/OSU.xlsx")

# OSU metadata
OSU_metadata_file <- paste0(manuscript_folder, "/Input_Data/OSU_Metadata.xlsx")
Broad_Annotations_file <- paste0(manuscript_folder, "/Input_Data/Broad_Metabolomics_Annotations.xlsx")

# NHS-Control Alignment files for HILIC-pos and C8-pos metabolites
HP_align_file <- paste0(manuscript_folder, "/Input_Data/HP_align.csv")
CP_align_file <- paste0(manuscript_folder, "/Input_Data/CP_align.csv")

# Standardized lipid names mapping file
LM_standardized_names_file <-  
  paste0(manuscript_folder, "/Input_Data/LIPIDMAPS_standardised_names.csv")

# Write out processed datasets

known_mets_control_file <- 
  paste0(manuscript_folder, 
         "/Processed_Data_and_Results/processed_control_known_metabolites.xlsx")


KD_knowns_corrected_data <- 
  paste0(manuscript_folder, 
         "/Processed_Data_and_Results/KD_known_metabolites_processed_batch_corrected.xlsx")
LFD_knowns_data <- 
  paste0(manuscript_folder, 
         "/Processed_Data_and_Results/LFD_known_metabolites_processed.xlsx")

############################
##### HELPER FUNCTIONS #####
############################


#' Filter out repeat metabolites based on HMDB
#'
#' For multiple measurements of the same metabolite (same HMDB)
#' keeps the measurement with the least QC Coefficient of Variance (CoV)
#'
#' @param SE Summarized Experiment object 
#'
#' @return SE without repeated HMDBs
#'
#' @noRd
filter_repeat_hmdbs <- function(SE){
  
  # Get which metabolites don't have an HMDB ID (no repeats)
  no_hmdb_id_indices <- which(is.na(rowData(SE)$HMDB_ID))
  
  # Get all unique hmdb ids
  unique_hmdbs <- rowData(SE) %>% 
    data.frame() %>% 
    filter(!is.na(HMDB_ID)) %>% 
    .$HMDB_ID %>% 
    unique()
  
  # Filter out repeat metabolites based on QC coefficient of variance
  # saving the metabolite with the lowest
  indices_to_keep <- lapply(unique_hmdbs, function(i){
    indices <- which(rowData(SE)$HMDB_ID == i)
    if (length(indices) > 1){
      indices[which(rowData(SE)$QC_cv[indices]==min(rowData(SE)$QC_cv[indices]))]
    }
    else{
      indices
    }
  }) %>% unlist()
  
  # Return the summarized experiment with no repeats
  SE[c(indices_to_keep,no_hmdb_id_indices),]
}



#' Load the Broad metabolomics dataset from excel
#'
#' Loads the metabolomics data from the Broad platform
#' in excel sheet file format
#'
#' @param D SummarizedExperiment object into which to write (optional)
#' @param file File name of the excel sheet to read in
#' @param sheet Sheet name to read in (these are generally separated based on LC platform)
#' @param study Which study this file comes from, either MSK or OSU
#'
#' @return SummarizedExperiment containing raw data, sample and molecular annotations
#'
#' @noRd
load_broad <-function (D, file, sheet, study = c("MSK", "OSU")) 
{
  # Read in the raw data
  result = list()
  result$info$file <- file
  result$info$sheet <- sheet
  raw = readxl::read_excel(path = file, sheet = sheet, col_names = F)
  
  # Get the metabolite and sample headers
  
  imetheader = min(which(!is.na(raw[, 1]))) # Where metabolite info begins
  imetlast = max(which(apply(is.na(raw), 1, sum) < dim(raw)[2])) # Where metabolite info ends
  isampheader = min(which(!is.na(raw[1, ]))) # Where sample info begins
  isamplast = max(which(apply(is.na(raw), 2, sum) < dim(raw)[1])) # Where sample info ends
  
  # Separate metabolite and sample name headers (found in the same cell)
  overl = stringr::str_replace(gsub("\\s+", " ", stringr::str_trim(raw[imetheader, 
                                                                       isampheader])), "B", "b")
  overl = strsplit(overl, " ")[[1]]
  overlmet = overl[1]
  overlsamp = overl[2]
  
  # Read in specifically metabolite information using indices found above
  result$metinfo <- readxl::read_excel(path = file, sheet = sheet, 
                                       col_names = T, range = readxl::cell_limits(ul = c(imetheader, 
                                                                                         1), lr = c(imetlast, isampheader)))
  result$metinfo <- as.data.frame(result$metinfo)
  colnames(result$metinfo) <- gsub(" ", "_", colnames(result$metinfo))
  colnames(result$metinfo)[ncol(result$metinfo)] = overlmet
  rownames(result$metinfo) <- result$metinfo$Name
  
  # Subselect raw table for sample information
  result$sampleinfo = data.frame(t(raw[1:imetheader, (isampheader + 
                                                        1):isamplast]), stringsAsFactors = F)
  colnames(result$sampleinfo) = c(as.vector(as.matrix(raw[1:imetheader - 
                                                            1, isampheader])), overlsamp)
  rownames(result$sampleinfo) = result$sampleinfo$Name
  result$sampleinfo %<>% dplyr::mutate_all(readr::parse_guess)
  
  # Save raw data using indices above with samples as rows and metabolites as columns
  result$data <- t(raw[(imetheader + 1):imetlast, (isampheader + 
                                                     1):isamplast])
  
  result$data <- as.data.frame(apply(result$data, 2, as.numeric))
  result$data[result$data==0] = NA
  
  # Annotate samples with their diet group
  result$sampleinfo$Diet<- get_sample_metadata(study, result$sampleinfo$Subject, type = "diet")
  result$sampleinfo$Sex<- get_sample_metadata(study, result$sampleinfo$Subject, type = "sex")
  
  
  # Annotate samples with their week (needs to be parsed out in OSU)
  if (study=="OSU"){
    result$sampleinfo$Week <- sapply(result$sampleinfo$Sample, function(i) strsplit(i,"K")[[1]][2])
  }
  
  # Convert to SummarizedExperiment
  colnames(result$data) = result$metinfo$Name
  rownames(result$data) = result$sampleinfo$Group
  D <- SummarizedExperiment(assay = t(result$data), colData = result$sampleinfo, 
                            rowData = result$metinfo)
  D
}


#' Assign diet to study subjects
#'
#' Check following file for OSU assignments:
#' Metagenics Data Update for Dr. Volek.xlsx
#' Following information is from Marcus (all subjects are Female):
#' Diet       Study ID
#'  KD          1
#'  SD          2
#'  SD          3
#'  KD          4
#'  KD          5
#'  KD          6
#'  SD          7
#'  SD          8
#'  KD          9
#'  KD          10
#'  KD          11
#'  KD          12
#' @param study Name of study (MSK or OSU)
#' @param name_vec Vector of subject ids
#' @param type Either "diet" or sex"
#'
#' @return Vector of diets (either KD or LFD) or sexes (F or M) corresponding to order of name_vec
#'
#' @noRd
get_sample_metadata <- function(
    study, 
    name_vec,
    type = c("diet", "sex")
){
  
  # If MSK, manually assign diet, all sex is F
  if (study=="MSK"){
    if (type == "diet"){
      MSK_SD<- c("\"02\"","\"03\"","\"07\"","\"08\"")
      return(unlist(lapply(name_vec,function(i){ if (i %in% MSK_SD){return ("LFD")}else{return("KD")}})))
    }
    else{
      return(rep("F", length(name_vec)))
    }
  }
  
  # If OSU, pull metadata from original datafiles and map
  else{
    OSU_metadata <- readxl::read_xlsx(OSU_metadata_file,
                                      sheet="Raw Chart - No Estimated Means") %>% 
      select(ID, GROUP, SEX) %>% 
      mutate(DIET = case_when(GROUP == "SUPP" ~ "KD", GROUP == "PLAC" ~ "KD",
                              GROUP == "LF" ~ "LFD"))
    OSU_metadata$Subject <- lapply(OSU_metadata$ID, function(i) as.numeric(str_split(i, " ")[[1]][2]))
    
    if (type == "diet"){
      return(OSU_metadata$DIET[match(name_vec, OSU_metadata$Subject)])
    }
    else{
      return(OSU_metadata$SEX[match(name_vec, OSU_metadata$Subject)])
    }
  }
  
  
}



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

    # Filter out metabolites with cv greater than 25%
    mt_modify_filter_features(QC_cv <0.25) %>%

    #Filter out samples with over 25% missingness
    mt_pre_filter_missingness(samp_max=0.25) %>%

    #Quotient Normalize
    mt_pre_norm_quot() %>%

    mt_plots_dilution_factor(in_col = "Subject") %>%

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



#' Data quality check
#'
#' Removes metabolites that show significant differences in the first week
#' 
#' @param combined_SE SummarizedExperiment containing processed OSU and MSK data

#'
#' @return SummarizedExperiment with metabolites that have first week differences removed
#'
#' @noRd
quality_check <- function(
    combined_SE
){
  
  OSU_dat<- combined_SE[,combined_SE$Sample_type=="Sample-STAK"]
  MSK_dat<- combined_SE[,combined_SE$Sample_type=="Sample-Makker"]
  
  # Don't want metabolites that significantly differ in first week (unreliable measurement)
  kd_0<-intersect(which(combined_SE$Diet =="KD"),which(combined_SE$Week ==0))
  sd_0<-intersect(which(combined_SE$Diet =="LFD"),which(combined_SE$Week ==0))
  init_diffs<-sapply(1:nrow(assay(combined_SE)), function(i) t.test(assay(combined_SE)[i, kd_0],assay(combined_SE)[i, sd_0])$p.value)
  same_week0<- which(p.adjust(init_diffs, method="bonferroni")>0.05)
  print(paste((nrow(combined_SE)-length(same_week0)), " metabolites significantly differ in the first week, removing"))
  
  
  combined_SE[same_week0,]
  
}



#' Loads, processes, and combines the control platforms (MSK and OSU)
#'
#' Filters for identified metabolites (that have HMDB IDs) leaving only known metabolites
#' Combines those metabolites with NHS identifiers
#'
#' @return SummarizedExperiment with preprocessed combined control platforms
#'
only_known_controls <-  function(){
  
  ## HILIC-pos
  # Read in HMDB information from NHS data for HP platform
  # This maps known metabolites from the NHS data to the control data
  HP_align <- read.csv(HP_align_file)
  
  MSK_hilic_pos<-load_broad(file=MSK_raw, 
                            sheet="HILIC-pos",
                            study= "MSK")
  
  rowData(MSK_hilic_pos)<- rowData(MSK_hilic_pos)[,c("HMDB_ID",
                                                     "Metabolite",
                                                     "Compound_ID")]
  
  names(rowData(MSK_hilic_pos))<- c("HMDB_ID","Metabolite", "MSK_ID")
  
  colData(MSK_hilic_pos)<- colData(MSK_hilic_pos)[,c("Sample_type","Sex",
                                                     "BMI","Subject",
                                                     "Name","Diet",'Week')]
  
  OSU_hilic_pos<-load_broad(file=OSU_raw, 
                            sheet="HILIC-pos",
                            study= "OSU")
  
  rowData(OSU_hilic_pos)<- rowData(OSU_hilic_pos)[,c("HMDB_ID",
                                                     "Metabolite", 
                                                     "Compound_ID")]
  
  names(rowData(OSU_hilic_pos))<- c("HMDB_ID","Metabolite", "OSU_ID")
  
  colData(OSU_hilic_pos)<- colData(OSU_hilic_pos)[,c("Sample_type","Sex",
                                                     "BMI","Subject",
                                                     "Name","Diet",'Week')]
  
  hilic_pos<- cbind(MSK_hilic_pos,OSU_hilic_pos) %>%  
    
    mt_pre_cv(Sample_type=="QC-pooled_plasma", out_col = "QC_cv")%>%
    
    maplet::mt_modify_filter_samples(filter=!is.na(Week))
  
  rowData(hilic_pos)$row_order <- 1:nrow(hilic_pos)
  merged_row_data <- merge(data.frame(rowData(hilic_pos)), HP_align, 
                           by.x = "OSU_ID", by.y = "Compound_ID_OSU",
                           all.x = T)
  rowData(hilic_pos) <- merged_row_data[order(merged_row_data$row_order),]
  
  # Coalesce HMDB IDs and Metabolites to keep knowns
  rowData(hilic_pos)$HMDB_ID <- dplyr::coalesce(rowData(hilic_pos)$HMDB_ID.x,
                                                rowData(hilic_pos)$HMDB_ID.y)
  
  rowData(hilic_pos)$Metabolite <- 
    dplyr::coalesce(rowData(hilic_pos)$Metabolite.x,
                    rowData(hilic_pos)$Metabolite.y)
  
  rowData(hilic_pos) %<>% data.frame() %>% 
    dplyr::select(OSU_ID, MSK_ID, QC_cv, HMDB_ID, Metabolite)
  
  hilic_pos <- hilic_pos[!is.na(rowData(hilic_pos)$HMDB_ID),]
  hilic_pos <- hilic_pos[!(rowData(hilic_pos)$HMDB_ID %in% 
                             c("NA","redundant ion", "internal standard")),]
  
  hilic_pos %<>% filter_repeat_hmdbs()
  rownames(hilic_pos) <- paste("hilic_pos", rowData(hilic_pos)$HMDB_ID, sep="_")
  rownames(rowData(hilic_pos)) <- rownames(hilic_pos)
  
  ##C8-pos
  
  # Read in HMDB information from NHS data for CP platform
  CP_align <- read.csv(CP_align_file)
  
  MSK_c8_pos<-load_broad(file=MSK_raw, 
                         sheet="C8-pos",
                         study= "MSK")
  
  rowData(MSK_c8_pos)<- rowData(MSK_c8_pos)[,c("HMDB_ID",
                                               "Metabolite",
                                               "Compound_ID")]
  
  names(rowData(MSK_c8_pos))<- c("HMDB_ID","Metabolite", "MSK_ID")
  
  colData(MSK_c8_pos)<- colData(MSK_c8_pos)[,c("Sample_type","Sex",
                                               "BMI", "Subject","Name",
                                               "Diet",'Week')]
  
  OSU_c8_pos<-load_broad(file=OSU_raw, 
                         sheet="C8-pos",
                         study= "OSU")
  
  rowData(OSU_c8_pos)<- rowData(OSU_c8_pos)[,c("Metabolite", "Compound_ID")]
  
  names(rowData(OSU_c8_pos))<- c("Metabolite", "OSU_ID")
  
  colData(OSU_c8_pos)<- colData(OSU_c8_pos)[,c("Sample_type","Sex",
                                               "BMI","Subject",
                                               "Name","Diet",'Week')]
  
  c8_pos<- cbind(MSK_c8_pos,OSU_c8_pos) %>%  
    
    mt_pre_cv(Sample_type=="QC-pooled_plasma", out_col = "QC_cv")%>%
    
    maplet::mt_modify_filter_samples(filter=!is.na(Week))
  
  rowData(c8_pos)$row_order <- 1:nrow(c8_pos)
  merged_row_data <- merge(data.frame(rowData(c8_pos)), HP_align, 
                           by.x = "OSU_ID", by.y = "Compound_ID_OSU", 
                           all.x = T)
  rowData(c8_pos) <- merged_row_data[order(merged_row_data$row_order),]
  
  # Coalesce HMDB IDs and Metabolites to keep knowns
  rowData(c8_pos)$HMDB_ID <- dplyr::coalesce(rowData(c8_pos)$HMDB_ID.x, 
                                             rowData(c8_pos)$HMDB_ID.y)
  
  rowData(c8_pos)$Metabolite <- dplyr::coalesce(rowData(c8_pos)$Metabolite.x, 
                                                rowData(c8_pos)$Metabolite.y)
  
  rowData(c8_pos) %<>% data.frame() %>% 
    dplyr::select(OSU_ID, MSK_ID, QC_cv, HMDB_ID, Metabolite)
  
  c8_pos <- c8_pos[!is.na(rowData(c8_pos)$HMDB_ID),]
  c8_pos <- c8_pos[!(rowData(c8_pos)$HMDB_ID %in% c("NA","redundant ion", 
                                                    "internal standard")),]
  
  c8_pos %<>% filter_repeat_hmdbs()
  rownames(c8_pos) <- paste("c8_pos", rowData(c8_pos)$HMDB_ID, sep="_")
  rownames(rowData(c8_pos)) <- rownames(c8_pos)
  
  ##HILIC-neg
  MSK_hilic_neg<-load_broad(file=MSK_raw, 
                            sheet="HILIC-neg",
                            study= "MSK")
  
  rowData(MSK_hilic_neg)<- rowData(MSK_hilic_neg)[,c("HMDB_ID",
                                                     "Metabolite",
                                                     "Compound_ID")]
  
  names(rowData(MSK_hilic_neg))<- c("HMDB_ID","Metabolite", "MSK_ID")
  
  colData(MSK_hilic_neg)<- colData(MSK_hilic_neg)[,c("Sample_type","Sex",
                                                     "BMI","Subject",
                                                     "Name","Diet",'Week')]
  
  OSU_hilic_neg<-load_broad(file=OSU_raw, 
                            sheet="HILIC-neg",
                            study= "OSU")
  
  rowData(OSU_hilic_neg)<- rowData(OSU_hilic_neg)[,c("Metabolite", "Compound_ID")]
  
  names(rowData(OSU_hilic_neg))<- c("Metabolite", "OSU_ID")
  
  colData(OSU_hilic_neg)<- colData(OSU_hilic_neg)[,c("Sample_type","Sex",
                                                     "BMI","Subject",
                                                     "Name","Diet",'Week')]
  identical(rownames(MSK_hilic_neg), rownames(OSU_hilic_neg))
  
  hilic_neg<- cbind(MSK_hilic_neg,OSU_hilic_neg) %>%  
    
    mt_pre_cv(Sample_type=="QC-pooled_plasma", out_col = "QC_cv")%>%
    
    maplet::mt_modify_filter_samples(filter=!is.na(Week))
  
  
  hilic_neg <- hilic_neg[!is.na(rowData(hilic_neg)$HMDB_ID),]
  hilic_neg <- hilic_neg[!(rowData(hilic_neg)$HMDB_ID %in% 
                             c("NA","redundant ion", "internal standard")),]
  
  hilic_neg %<>% filter_repeat_hmdbs()
  
  rownames(hilic_neg) <- paste("hilic_neg", rowData(hilic_neg)$HMDB_ID, sep="_")
  rownames(rowData(hilic_neg)) <- rownames(hilic_neg)
  
  ##C18-neg
  MSK_c18_neg<-load_broad(file=MSK_raw, 
                          sheet="C18-neg",
                          study= "MSK")
  
  rowData(MSK_c18_neg)<- rowData(MSK_c18_neg)[,c("HMDB_ID","Metabolite", 
                                                 "Compound_ID")]
  
  names(rowData(MSK_c18_neg))<- c("HMDB_ID","Metabolite", "MSK_ID")
  
  colData(MSK_c18_neg)<- colData(MSK_c18_neg)[,c("Sample_type","Sex",
                                                 "BMI","Subject",
                                                 "Name","Diet",'Week')]
  
  OSU_c18_neg<-load_broad(file=OSU_raw, 
                          sheet="C18-neg",
                          study= "OSU")
  
  rowData(OSU_c18_neg)<- rowData(OSU_c18_neg)[,c("Metabolite", "Compound_ID")]
  
  names(rowData(OSU_c18_neg))<- c("Metabolite", "OSU_ID")
  
  colData(OSU_c18_neg)<- colData(OSU_c18_neg)[,c("Sample_type","Sex",
                                                 "BMI","Subject","Name",
                                                 "Diet",'Week')]
  
  c18_neg<- cbind(MSK_c18_neg,OSU_c18_neg) %>%  
    
    mt_pre_cv(Sample_type=="QC-pooled_plasma", out_col = "QC_cv")%>%
    
    maplet::mt_modify_filter_samples(filter=!is.na(Week))
  
  c18_neg <- c18_neg[!is.na(rowData(c18_neg)$HMDB_ID),]
  c18_neg <- c18_neg[!(rowData(c18_neg)$HMDB_ID %in% c("NA","redundant ion", 
                                                       "internal standard")),]
  
  c18_neg %<>% filter_repeat_hmdbs()
  
  rownames(c18_neg) <- paste("c18_neg", rowData(c18_neg)$HMDB_ID, sep="_")
  rownames(rowData(c18_neg)) <- rownames(c18_neg)
  
  combined_assay<- do.call("rbind", list(assay(hilic_pos),
                                         assay(c8_pos),
                                         assay(hilic_neg),
                                         assay(c18_neg)))
  
  combined_rows <- do.call("rbind",list(data.frame(rowData(hilic_pos)),
                                        data.frame(rowData(c8_pos)),
                                        data.frame(rowData(hilic_neg)),
                                        data.frame(rowData(c18_neg)))) %>% 
    left_join(read.csv(LM_standardized_names_file) %>% 
                distinct(), by = "Metabolite")
  
  combined_SE<-SummarizedExperiment(assays=combined_assay, 
                                    colData= colData(hilic_pos), 
                                    rowData=combined_rows) %>% 
    # Filter out subject 11 from MSK, only has first week
    maplet::mt_modify_filter_samples(filter=  Subject!= "\"11\"") %>% 
    preprocess_broad_data() %>% 
    mt_pre_confounding_correction(formula = ~BMI+Sex) %>% 
    quality_check()
  
  # Unify week  and sex numbering between cohorts
  combined_SE %<>% 
    mt_anno_mutate(anno_type = "samples", col_name = "Week",
                   term = ifelse(Sample_type=="Sample-Makker", 
                                 as.character(as.integer(Week)-1), Week)) %>% 
    mt_anno_mutate(anno_type = "samples", col_name = "Sex",
                   term = ifelse(Sex == "FALSE", "F", Sex))
  
  combined_SE %>% 
    filter_repeat_hmdbs() %>% 
    mt_anno_xls(file = Broad_Annotations_file, sheet = 1, 
                anno_type = 'features', anno_id_col = "HMDB_ID") %>% 
    mt_anno_mutate(anno_type = "features", col_name = 'Fat_type',
                   term =  ifelse(Fat_type=="NA", NA, Fat_type)) %>% 
    mt_anno_mutate(anno_type = "features", 
                   col_name = 'Molecular.Annotation', 
                   term =  ifelse(Molecular.Annotation=="NA", 
                                  NA, Molecular.Annotation))
  
  
}


##### Batch Correction for KD #####
combined_control_knowns_D <- only_known_controls() 
KD_knowns_D <- combined_control_knowns_D %>% 
  mt_modify_filter_samples(Diet == "KD")

annotation_rows <- KD_knowns_D %>% 
  colData() %>% 
  data.frame() %>% 
  dplyr::select(Sample_type,Sex,BMI, Subject, Week) %>% 
  mutate(Cohort = as.factor(Sample_type)) %>% 
  mutate(Sex = as.factor(Sex))


# Identifying and correcting batch effects
batch <- annotation_rows$Cohort 

# Apply ComBat
adjusted_data <- ComBat(dat=assay(KD_knowns_D), batch=batch, 
                        par.prior=TRUE, prior.plots=FALSE)


KD_batch_corrected <- SummarizedExperiment(assays = adjusted_data,
                                           colData = colData(KD_knowns_D),
                                           rowData = rowData(KD_knowns_D))

##### Data Saving #####

# Only known metabolites without batch correction
mt_write_se_xls(combined_control_knowns_D, file=known_mets_control_file)



# Only LFD samples
LFD_D <- combined_control_knowns_D %>% 
  mt_modify_filter_samples(Diet=="LFD") %>% 
  mt_modify_filter_samples(Sample_type == "Sample-STAK")

mt_write_se_xls(LFD_D, file=LFD_knowns_data)

# KD samples, batch corrected

mt_write_se_xls(KD_batch_corrected, file=KD_knowns_corrected_data)
