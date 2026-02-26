rm(list = ls())

library(dplyr)
library(maplet)
library(reshape)
library(readxl)
library(sva)
library(stringr)
library(ggplot2)
library(glmnet)
library(caret)
library(nlme)
library(parallel)
library(ggpubr)
library(openxlsx)
library(chanmetab)
library(sva)
library(sas7bdat)
library(tidyverse)
library(survival) 
library(mgcv)
library(data.table)
library(Biobase)
library(glue)
library(gt)
library(gtsummary)
library(viridis)
library(ggrepel)
library(cowplot)
library(grid)

setwd("~/keto-metabolomics/SchweickartAnnalise_HypocaloricDiets_022025")

source("Scripts/internal_functions.R")


manuscript_folder <- 
  "/udd/nhast/keto-metabolomics/SchweickartAnnalise_HypocaloricDiets_022025"


# Script Run order:

# 01: Preprocess and save controlled trials data, filter metabolites
#     Outputs ExpressionSets in excel sheets to 
#     "Processed_Data_and_Results" folder
source("Scripts/01_Preprocess_Controlled_Trials_data.R")

# 01: Find differential metabolites by diet
#     Outputs betas in excel sheets to 
#     "Processed_Data_and_Results" folder
source("Scripts/02_get_differential_metabolites.R")

# 03: Preprocess and save NHS data, filter metabolites
#     Outputs ExpressionSets in excel sheets to 
#     "Processed_Data_and_Results" folder
source("Scripts/03_Preprocess_NHS_data.R")

# 04: Calculate diet score for each subject, combine with metadata
#     Outputs combined diet scores and metadata for both NHS1 and NHS2 and
#     KD and LFD to "Processed_Data_and_Results" folder
source("Scripts/04_load_save_data.R")

# 05: Runs MWAS of metabolites with breast cancer outcome
#     Outputs results to excel file in 
#     "Processed_Data_and_Results"
source("Scripts/05_run_mwas.R")

# 06: Calculation of concordance values for each metabolite in diet scores
#     And comparison of concordance values with fold changes
#     Outputs concordance information in .RData file to
#     "Processed_Data_and_Results"

source("Scripts/06_calculate_metabolite_concordance.R")

# Create Manuscript Tables
rmarkdown::render("Schweickart_2025_Tables.Rmd")

# Create Manuscript Figures
rmarkdown::render("Schweickart_2025_Figures.Rmd")

# Get Manuscript Data in Text
rmarkdown::render("Schweickart_2025_Data_In_Text.Rmd")
