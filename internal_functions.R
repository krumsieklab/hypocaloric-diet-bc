# Functions for loading data and running association analysis
#remotes::install_version("glue", version = "1.7.0")
library(glue)
library(sas7bdat)
library(tidyverse)
library(readxl)
library(survival) 
library(mgcv)
library(openxlsx)
library(data.table)
library(Biobase)
library(dplyr)
library(gt)
library(metafor)
library(magrittr)
library(officer)
library(htmltools)
library(maplet)


#' Annotate control platforms into first week (pre diet) and last week (post diet)
#'
#'
#' @param D SummarizedExperiment of metabolites with "Week" column in colData
#' 
#' @return Vector with 0 for first week, 1 for last week, and NA for all other weeks aligned with SE Week column
#'
#'
assign_pre_post <- function(D){
  pre_post <- c()
  for (i in 1:nrow(data.frame(colData(D)))){
    subj <- colData(D)$Subject[i]
    weeks <- colData(D)[which(colData(D)$Subject==subj), "Week"]
    if (colData(D)$Week[i] == min(weeks)){
      pre_post <- c(pre_post, 0)
    }
    else if (colData(D)$Week[i] == max(weeks)){
      pre_post <- c(pre_post, 1)
    }
    else{
      pre_post <- c(pre_post, NA)
    }
  }
  pre_post
}


#' Get metabolite fold change from pre- and post- data points
#'
#'
#' @param D SummarizedExperiment of metabolites with "Subject" and "pre_post" already in the colData
#' 
#' @return Dataframe with the fold changes for each metabolite for each subject
#'
#'
get_diet_fc <-function(D){
  orig_df <- cbind(data.frame(colData(D)[,c("Subject", "pre_post")]),t(assay(D)))
  
  diff_df <- lapply(unique(orig_df$Subject), function(i){
    sub_df <- orig_df %>% filter(Subject==i)
    post_row = which(sub_df$pre_post == 1)
    pre_row = which(sub_df$pre_post == 0)
    diffs <- sub_df[post_row, rownames(assay(D))] - sub_df[pre_row, rownames(assay(D))]
    data.frame(diffs)
  }) %>% bind_rows()
  rownames(diff_df)<-unique(orig_df$Subject)
  diff_df
}

#' Get BMI change from pre- and post- data points
#'
#'
#' @param D SummarizedExperiment of metabolites with "Subject" and "pre_post" already in the colData
#' 
#' @return Dataframe with the BMI changes for each subject
#'
#'
get_BMI_change <-function(D){
  orig_df <- data.frame(colData(D)[,c("Subject", "pre_post", "BMI")])
  
  diff_df <- lapply(unique(orig_df$Subject), function(i){
    sub_df <- orig_df %>% filter(Subject==i)
    post_row = which(sub_df$pre_post == 1)
    pre_row = which(sub_df$pre_post == 0)
    BMI_change <- sub_df[post_row, "BMI"] - sub_df[pre_row, "BMI"]
    data.frame(BMI_change)
  }) %>% bind_rows()
  rownames(diff_df)<-unique(orig_df$Subject)
  diff_df
}

# Get the quartile and median for each score in the NHSII cohort
get_quartiles_and_medians <- function(data_for_assoc, controls_only = T){
 
   # Set quartiles and medians for p-trend
  if(controls_only){
    control=data_for_assoc[data_for_assoc$case==0,]
    first_quart <- summary(control$beta.score)["1st Qu."]
    med <- summary(control$beta.score)["Median"]
    third_quart <- summary(control$beta.score)["3rd Qu."]
  }
  else{
    first_quart <- summary(data_for_assoc$beta.score)["1st Qu."]
    med <- summary(data_for_assoc$beta.score)["Median"]
    third_quart <- summary(data_for_assoc$beta.score)["3rd Qu."]
  }

  
  # Assign quartiles ("beta.score.cat") for all samples
  data_for_assoc$beta.score.cat <- 
    case_when(data_for_assoc$beta.score < first_quart ~ 1,
              data_for_assoc$beta.score < med ~ 2,
              data_for_assoc$beta.score < third_quart ~ 3,
              data_for_assoc$beta.score > third_quart ~ 4)
  
  # Assign each sample the median of their quartile
  median_first_quart <- data_for_assoc %>% filter(beta.score.cat == 1) %>% 
    .$beta.score %>% median()
  median_second_quart <- data_for_assoc %>% filter(beta.score.cat == 2) %>% 
    .$beta.score %>% median()
  median_third_quart <- data_for_assoc %>% filter(beta.score.cat == 3) %>% 
    .$beta.score %>% median()
  median_fourth_quart <- data_for_assoc %>% filter(beta.score.cat == 4) %>% 
    .$beta.score %>% median()
  
  # Add to new variable beta.score.cat.median 
  data_for_assoc$beta.score.cat.median <- 
    case_when(data_for_assoc$beta.score < first_quart ~ median_first_quart,
              data_for_assoc$beta.score.cat==2 ~ median_second_quart,
              data_for_assoc$beta.score.cat==3 ~ median_third_quart,
              data_for_assoc$beta.score.cat==4 ~ median_fourth_quart)
  data_for_assoc
}


# Performs conditional logistic regression analysis 
log_reg_association <- function(
    data_for_assoc,
    model_type = c("linear", "quartile", "trend"),
    bmi, # Include or exclude BMI covariates
    conditional, # Conditional or unconditional
    meno = T
    ){
  
  # Conditional vs unconditional covars
  conditional_covars = "strata(matchid)"
  if(!conditional){
    # Matching factors that are removed in conditional logistic regression
    if(meno){
      conditional_covars = 
        "ageyr + as.factor(fast) + as.factor(menopmh) + blddate"
    }
    else{
      conditional_covars = 
        "ageyr + as.factor(fast) + blddate"
    }
  }

  # Set covariates based on inclusion of BMI variables
  if(bmi){
    covars = 
    'act + as.factor(current.smk) + ahei_noal + alco + menage + 
    as.factor(afbpar) + familyhx + bbddx + bmi18 + as.factor(changecat)'
  }
  else{
    covars = 'act + as.factor(current.smk) + ahei_noal + alco +  
    menage + as.factor(afbpar) + familyhx + bbddx'
  }

  # Set score type (continuous, categorical, medians)
  if(model_type == "linear"){ x = "scale(beta.score)" }
  else if(model_type == "quartile"){ x = "as.factor(beta.score.cat)" }
  else{ x = "beta.score.cat.median" }
  
  form <- as.formula(glue("case ~ {x} + {conditional_covars} + {covars}"))
  
  # If conditional, use clogit with matchid
  if (conditional){  res <- clogit(formula = form, data = data_for_assoc) }
  
  # If not conditional, use generalized linear model
  else{ res <- glm(formula = form, family = binomial, data = data_for_assoc) }
  
  summary(res)
}


# Extract stats from regression model
get_stats <- function(res, 
                      model_type = "linear", 
                      conditional = T){
  
  # If conditional, clogit was used
  if (conditional){
    if (model_type == "linear"){
      p_value <- res %>% .$coefficients %>% .[1, "Pr(>|z|)"]%>% round(3)
      coef <- res %>% .$conf.int %>% .[1, "exp(coef)"]%>% round(2)
      lower <- res %>% .$conf.int %>% .[1, "lower .95"]%>% round(2)
      upper <- res %>% .$conf.int %>% .[1, "upper .95"]%>% round(2)
      res_list <- list(glue("{coef} ({lower}, {upper})"), p_value)
      names(res_list) <- c("Per SD increase", "P-value")
    }
    else{
      Q2_lower <- res %>% .$conf.int %>% .[1,"lower .95"] %>% round(2)
      Q2_upper <- res %>% .$conf.int %>% .[1, "upper .95"]%>% round(2)
      Q2_est <-  res %>% .$conf.int %>% .[1, "exp(coef)"]%>% round(2)
      Q3_lower <- res %>% .$conf.int %>% .[2, "lower .95"]%>% round(2)
      Q3_upper <- res %>% .$conf.int %>% .[2, "upper .95"]%>% round(2)
      Q3_est <- res %>% .$conf.int %>% .[2, "exp(coef)"]%>% round(2)
      Q4_lower <- res %>% .$conf.int %>% .[3, "lower .95"]%>% round(2)
      Q4_upper <- res %>% .$conf.int %>% .[3, "upper .95"]%>% round(2)
      Q4_est <- res %>% .$conf.int %>% .[3, "exp(coef)"]%>% round(2)
      res_list <- list(glue("{Q2_est} ({Q2_lower}, {Q2_upper})"),
                       glue("{Q3_est} ({Q3_lower}, {Q3_upper})"),
                       glue("{Q4_est} ({Q4_lower}, {Q4_upper})"))
      names(res_list) <- c("Q2", "Q3", "Q4")
    }
  }
  
  # Otherwise, glm was used
  else{
    if (model_type == "linear"){
      p_value <- res %>% .$coefficients %>% .[2, "Pr(>|z|)"] %>% round(3)
      est <-  res %>% .$coefficients %>% .[2, "Estimate"] %>% round(2)
      std_err <- res %>% .$coefficients %>% .[2, "Std. Error"] %>% round(2)
      
      coef <- exp(est) %>%  round(2)
      lower <- (est - 1.96 * std_err) %>% exp() %>%  round(2) 
      upper <-(est + 1.96 * std_err) %>% exp() %>%  round(2) 
      
      res_list <- list(glue("{coef} ({lower}, {upper})"), p_value)
      names(res_list) <- c("Per SD increase", "P-value")
    }
    else{
      Q2_est <-  res %>% .$coefficients %>% .[2, "Estimate"]
      Q2_se <-  res %>% .$coefficients %>% .[2, "Std. Error"]
      
      Q2_lower <- (Q2_est - 1.96 * Q2_se) %>% exp() %>%  round(2) 
      Q2_upper <- (Q2_est + 1.96 * Q2_se) %>% exp() %>%  round(2) 
      Q2_est_exp <-  Q2_est %>% exp() %>% round(2)
      
      Q3_est <-  res %>% .$coefficients %>% .[3, "Estimate"]
      Q3_se <-  res %>% .$coefficients %>% .[3, "Std. Error"]
      
      Q3_lower <- (Q3_est - 1.96 * Q3_se) %>% exp() %>%  round(2) 
      Q3_upper <- (Q3_est + 1.96 * Q3_se) %>% exp() %>%  round(2) 
      Q3_est_exp <-  Q3_est %>% exp() %>% round(2)
      
      Q4_est <-  res %>% .$coefficients %>% .[4, "Estimate"]
      Q4_se <-  res %>% .$coefficients %>% .[4, "Std. Error"]
      
      Q4_lower <- (Q4_est - 1.96 * Q4_se) %>% exp() %>%  round(2) 
      Q4_upper <- (Q4_est + 1.96 * Q4_se) %>% exp() %>%  round(2) 
      Q4_est_exp <-  Q4_est %>% exp() %>% round(2)
      
      res_list <- list(glue("{Q2_est_exp} ({Q2_lower}, {Q2_upper})"),
                       glue("{Q3_est_exp} ({Q3_lower}, {Q3_upper})"),
                       glue("{Q4_est_exp} ({Q4_lower}, {Q4_upper})"))
      names(res_list) <- c("Q2", "Q3", "Q4")
    }
  }
  
  res_list
}

table_output <- function(data_for_assoc, title, 
                         outfile,
                         controls_only = T,
                         conditional = T, # Conditional or unconditional
                         stratification = c("ER", "meno_at_bd", 
                                            "meno_at_dx","invasive"),
                         meno=T,
                         save = T){
  if (!conditional){
    # Stratify based on stratification variable
    if(stratification == "ER"){
      # ER positive cases
      data_for_assoc_1 <- rbind(data_for_assoc %>% filter(case == 0),
                        data_for_assoc %>% filter(erpos == 1)) %>% 
        get_quartiles_and_medians(controls_only)
      
      # ER negative cases
      data_for_assoc_2<- rbind(data_for_assoc %>% filter(case == 0),
                               data_for_assoc %>% filter(erpos == 0)) %>% 
        get_quartiles_and_medians(controls_only)
      
      lab1 = "ER + "
      lab2 = "ER - "
    }
    else if(stratification == "meno_blood_draw"){
      # Premenopausal at blood draw
      data_for_assoc_1 <- 
        rbind(data_for_assoc %>% filter(case == 0),
              data_for_assoc %>% filter(case == 1 & menopmh == "pre")) %>% 
        get_quartiles_and_medians(controls_only)
      
      # Postmenopausal at blood draw
      data_for_assoc_2 <- 
        rbind(data_for_assoc %>% filter(case == 0),
              data_for_assoc %>% filter(case == 1 & menopmh != "pre")) %>% 
        get_quartiles_and_medians(controls_only)

      lab1 = "Pre-menopausal at blood draw "
      lab2 = "Post-menopausal at blood draw "
    }
    else if(stratification == "meno_dx"){
      # Premenopausal at diagnosis in cases
      data_for_assoc_1 <- 
        rbind(data_for_assoc %>% filter(case == 0),
              data_for_assoc %>% filter(case == 1 & menopmh_dx == 1)) %>% 
        get_quartiles_and_medians(controls_only)
      
      # Postmenopausal at diagnosis in cases
      data_for_assoc_2 <- 
        rbind(data_for_assoc %>% filter(case == 0),
              data_for_assoc %>% filter(case == 1 & menopmh_dx != 1)) %>% 
        get_quartiles_and_medians(controls_only)
      
      lab1 = "Pre-menopausal cancer "
      lab2 = "Post-menopausal cancer "
    }
    #invasive
    else {
      # In situ breast cancer in cases
      data_for_assoc_1 <- 
        rbind(data_for_assoc %>% filter(case == 0),
              data_for_assoc %>% filter(case == 1 & invas == 0)) %>% 
        get_quartiles_and_medians(controls_only)
      
      # Invasive breast cancer in cases
      data_for_assoc_2<- 
        rbind(data_for_assoc %>% filter(case == 0),
              data_for_assoc %>% filter(case == 1 & invas == 1)) %>% 
        get_quartiles_and_medians(controls_only)
      
      lab1 = "In situ breast cancer "
      lab2 = "Invasive breast cancer "
    }
    
    case_num_1 = table(data_for_assoc_1$case)[2]
    control_num_1 =  table(data_for_assoc_1$case)[1]
    
    table_df_1 <- make_base_stats_table(data_for_assoc_1, 
                                        conditional = F, meno=meno)
    
    case_num_2 = table(data_for_assoc_2$case)[2]
    control_num_2 =  table(data_for_assoc_2$case)[1]
    
    table_df_2 <- make_base_stats_table(data_for_assoc_2, 
                                        conditional = F, meno=meno)

    title1 = glue("{lab1} ({case_num_1} cases/{control_num_1} controls)")
    title2 = glue("{lab2} ({case_num_2} cases/{control_num_2} controls)")
    
    # Continuous analysis no bmi
    continuous_res_no_bmi_1 <- data_for_assoc_1 %>% 
      log_reg_association(model_type = "linear", bmi = F,
                          conditional = conditional, meno=meno) 
    
    continuous_res_w_bmi_1 <- data_for_assoc_1 %>% 
      log_reg_association(model_type = "linear", bmi = T,
                          conditional = conditional, meno=meno) 
    
    # Continuous analysis with bmi
    continuous_res_no_bmi_2 <- data_for_assoc_2 %>% 
      log_reg_association(model_type = "linear", bmi = F,
                          conditional = conditional, meno=meno)
    
    continuous_res_w_bmi_2 <- data_for_assoc_2 %>% 
      log_reg_association(model_type = "linear", bmi = T,
                          conditional = conditional, meno=meno)
    
    het_p_no_bmi <- get_heterogeneity_p_value(continuous_res_no_bmi_1,
                                              continuous_res_no_bmi_2)
    het_p_w_bmi <- get_heterogeneity_p_value(continuous_res_w_bmi_1,
                                              continuous_res_w_bmi_2)
    table_df_1$group = rep(title1,2)
    table_df_1$row = rownames(table_df_1)
    
    table_df_2$group = rep(title2,2)
    table_df_2$row = rownames(table_df_2)
    
    combined <- rbind(table_df_1, table_df_2)
    combined$het_p <- c(het_p_no_bmi, het_p_w_bmi, NA, NA)
    
    combined %<>% gt(rowname_col = "row", groupname_col = "group") %>% 
      tab_spanner(label = "Categorical OR (95% CI)", 
                  columns = names(combined)[1:5]) %>% 
      tab_spanner(label = "Continuous OR (95% CI)", 
                  columns =names(combined)[6:8]) %>% 
      tab_header(title = md(title),
                 subtitle = paste0("N = ", nrow(data_for_assoc))) 
    if(save){
      combined %>% 
        gtsave(filename = outfile)
    }
    else{combined}
    
  }
  
  else{
    table_df <- make_base_stats_table(
      data_for_assoc %>% 
       get_quartiles_and_medians(controls_only))
    
    table_df %<>% gt(rownames_to_stub = T) %>% 
      tab_spanner(label = "Categorical OR (95% CI)", 
                  columns = names(table_df)[1:5]) %>% 
      tab_spanner(label = "Continuous OR (95% CI)", 
                  columns =names(table_df)[6:7]) %>% 
      tab_stubhead(label = "Quartiles (cases/controls)") %>% 
      tab_header(title = md(title),
                 subtitle = paste0("N = ", nrow(data_for_assoc))) 
    if(save){
      table_df %>% 
        gtsave(filename = outfile)
    }
    else{table_df}
  }
  
}

# Makes a table with stats from quartile, trend, and linear regressions
make_base_stats_table <- function(data_for_assoc, conditional = T, meno = T){
  
  # Continuous analysis
  continuous_res_no_bmi <- data_for_assoc %>% 
    log_reg_association(model_type = "linear", bmi = F,
                        conditional = conditional, meno = meno) %>% 
    get_stats(conditional = conditional)
  
  continuous_res_w_bmi <- data_for_assoc %>% 
    log_reg_association(model_type = "linear", bmi = T,
                        conditional = conditional, meno = meno) %>% 
    get_stats(conditional= conditional)
  
  # P-trend analysis
  trend_p_no_bmi <- data_for_assoc %>% 
    log_reg_association(model_type = "trend", bmi = F,
                        conditional = conditional, meno = meno) %>% 
    
    .$coefficients %>% .["beta.score.cat.median", "Pr(>|z|)"] %>% round(2)
  
  trend_p_w_bmi <- data_for_assoc %>% 
    log_reg_association(model_type = "trend", bmi = T,
                        conditional = conditional, meno = meno) %>%
    
    .$coefficients %>% .["beta.score.cat.median", "Pr(>|z|)"] %>% round(2)
  
  # Quartile analysis
  quart_res_no_bmi <- data_for_assoc %>% 
    log_reg_association(model_type = "quartile", bmi = F,
                        conditional = conditional, meno = meno) %>% 
    
    get_stats(model_type = "quartile", conditional = conditional)
  
  quart_res_w_bmi <- data_for_assoc %>% 
    log_reg_association(model_type = "quartile", bmi = T,
                        conditional = conditional, meno = meno) %>% 
    
    get_stats(model_type = "quartile", conditional = conditional)
  
  # Get cases and controls for each quartile
  quart_counts <- table(data_for_assoc %>% select(beta.score.cat, case)) %>% 
    data.frame() %>% group_by(beta.score.cat) %>% arrange(desc(case)) %>% 
    summarise("N cases/controls" = paste(Freq, collapse="/"))
  
  table_df <- cbind(
    data.frame(Q1 = c("1.0 (ref)", "1.0 (ref)")),
    bind_rows(quart_res_no_bmi, quart_res_w_bmi),
    data.frame("P value of trend" = c(trend_p_no_bmi, trend_p_w_bmi)),
    bind_rows(continuous_res_no_bmi, continuous_res_w_bmi))
  
  # names(table_df)<- c(paste0(names(table_df)[1:4], " (",
  #                            quart_counts$`N cases/controls`, ")"), 
  #                     "P-value of trend",
  #                     names(table_df)[6:7])
  
  names(table_df)<- c("Q1", "Q2","Q3","Q4", 
                      "P-value of trend",
                      names(table_df)[6:7])
  rownames(table_df)<- c("Model 1", "Model 2")
  
  table_df
}

# Perform Cochran's Q test on two results to test for heterogeneity
get_heterogeneity_p_value <- function(res_1,res_2){
  
  est_1 <-  res_1 %>% .$coefficients %>% .[2, "Estimate"] %>% round(2)
  std_err_1 <- res_1 %>% .$coefficients %>% .[2, "Std. Error"] %>% round(2)
  
  est_2 <-  res_2 %>% .$coefficients %>% .[2, "Estimate"] %>% round(2)
  std_err_2 <- res_2 %>% .$coefficients %>% .[2, "Std. Error"] %>% round(2)
  
  rma(yi = c(est_1, est_2), sei = c(std_err_1, std_err_2), method = "FE") %>% 
    .$QEp %>% round(2)
  
}

subgroup_log_reg_association <- function(
    data_for_assoc, bmi = T){
  # Set covariates based on inclusion of BMI variables
  if(bmi){
    covars = 'ageyr + as.factor(fast) + as.factor(menopmh) + act + 
      as.factor(current.smk) + ahei_noal + alco + blddate +  menage + 
      as.factor(afbpar) + familyhx + bbddx + bmi18 +as.factor(changecat)'
  }
  else{
    covars = 'ageyr + as.factor(fast) + as.factor(menopmh) + act + 
      as.factor(current.smk) + ahei_noal + alco + blddate +  menage + 
      as.factor(afbpar) + familyhx + bbddx '
  }
  
  # Set score type 
  x = "scale(beta.score)" 
  
  form <- as.formula(glue("case ~ {x} + {covars}"))
  res <- glm(form, data_for_assoc, family = "binomial")
  res
}

get_subgroup_stats <- function(res){
  p_val <- summary(res)$coefficients["scale(beta.score)", "Pr(>|z|)"] %>% 
    round(2)
  
  estimate <- summary(res)$coefficients["scale(beta.score)", "Estimate"]
  std_err <- summary(res)$coefficients["scale(beta.score)", "Std. Error"]
  
  OR <- exp(coef(res))["scale(beta.score)"] %>%  round(2)
  OR_lower <- 
    exp(summary(res)$coefficients[,1] - 1.96*summary(res)$coefficients[,2]) %>% 
    .["scale(beta.score)"] %>% 
  round(2)
  OR_upper <- 
    exp(summary(res)$coefficients[,1] + 1.96*summary(res)$coefficients[,2]) %>% 
    .["scale(beta.score)"] %>% 
    round(2)
  
  data.frame("P-value" = p_val,
             Estimate = estimate,
             std_error = std_err,
             "OR (95% CI)" = glue("{OR} ({OR_lower}, {OR_upper})"))
  
}


subgroup_table_output <- function(data_for_assoc, title, results_folder,
                                  outfile, meno_or_er = "ER", plot = F,
                                  save = T){
  
  # ER status
  if(meno_or_er == "ER"){
    dat_erpos<- rbind(data_for_assoc %>% filter(case == 0),
                      data_for_assoc %>% filter(erpos == 1))
    
    erpos_cases = table(dat_erpos$case)[2]
    erpos_controls =  table(dat_erpos$case)[1]
    
    no_bmi_res_1 <- subgroup_log_reg_association(dat_erpos,
                                        bmi=F) %>% 
      get_subgroup_stats()
    bmi_res_1 <- subgroup_log_reg_association(dat_erpos, bmi=T) %>% 
      get_subgroup_stats()
    
    dat_erneg<- rbind(data_for_assoc %>% filter(case == 0),
                      data_for_assoc %>% filter(erpos == 0))
    
    erneg_cases = table(dat_erneg$case)[2]
    erneg_controls =  table(dat_erneg$case)[1]
    
    no_bmi_res_2 <- subgroup_log_reg_association(dat_erneg, bmi=F) %>% 
      get_subgroup_stats()
    bmi_res_2 <- subgroup_log_reg_association(dat_erneg, bmi=T) %>% 
      get_subgroup_stats()
    
    lab1 = glue("ER + ({erpos_cases} cases/{erpos_controls} controls)")
    lab2 = glue("ER - ({erneg_cases} cases/{erneg_controls} controls)")
    
  }
  else if(meno_or_er == "meno_blood_draw"){
    dat_pre<- rbind(data_for_assoc %>% filter(case == 0),
                    data_for_assoc %>% filter(case == 1 & menopmh == "pre"))
    
    pre_cases = table(dat_pre$case)[2]
    pre_controls =  table(dat_pre$case)[1]
    
    no_bmi_res_1 <- subgroup_log_reg_association(dat_pre, bmi=F) %>% 
      get_subgroup_stats()
    bmi_res_1 <- subgroup_log_reg_association(dat_pre, bmi=T) %>% 
      get_subgroup_stats()
    
    dat_post<- rbind(data_for_assoc %>% filter(case == 0),
                     data_for_assoc %>% filter(case == 1 & menopmh != "pre"))
    
    post_cases = table(dat_post$case)[2]
    post_controls =  table(dat_post$case)[1]
    
    no_bmi_res_2 <- subgroup_log_reg_association(dat_post, bmi=F) %>% 
      get_subgroup_stats()
    bmi_res_2 <- subgroup_log_reg_association(dat_post, bmi=T) %>% 
      get_subgroup_stats()
    
    lab1 = 
      glue("Pre-menopausal at blood draw \n({pre_cases} cases/{pre_controls} controls)")
    lab2 = 
      glue("Post-menopausal at blood draw \n({post_cases} cases/{post_controls} controls)")
  }
  else if(meno_or_er == "meno_dx"){
    dat_pre<- rbind(data_for_assoc %>% filter(case == 0),
                    data_for_assoc %>% filter(case == 1 & menopmh_at_dx == 1))
    
    pre_cases = table(dat_pre$case)[2]
    pre_controls =  table(dat_pre$case)[1]
    
    no_bmi_res_1 <- subgroup_log_reg_association(dat_pre, bmi=F) %>% 
      get_subgroup_stats()
    bmi_res_1 <- subgroup_log_reg_association(dat_pre, bmi=T) %>% 
      get_subgroup_stats()
    
    dat_post<- rbind(data_for_assoc %>% filter(case == 0),
                     data_for_assoc %>% filter(case == 1 & menopmh_at_dx != 1))
    
    post_cases = table(dat_post$case)[2]
    post_controls =  table(dat_post$case)[1]
    
    no_bmi_res_2 <- subgroup_log_reg_association(dat_post, bmi=F) %>% 
      get_subgroup_stats()
    bmi_res_2 <- subgroup_log_reg_association(dat_post, bmi=T) %>% 
      get_subgroup_stats()
    
    lab1 = 
      glue("Pre-menopausal cancer \n({pre_cases} cases/{pre_controls} controls)")
    lab2 = 
      glue("Post-menopausal cancer \n({post_cases} cases/{post_controls} controls)")
  }
  #invasive
  else {
    dat_non<- rbind(data_for_assoc %>% filter(case == 0),
                    data_for_assoc %>% filter(case == 1 & invas == 0))
    
    non_cases = table(dat_non$case)[2]
    non_controls =  table(dat_non$case)[1]
    
    no_bmi_res_1 <- subgroup_log_reg_association(dat_non, bmi=F) %>% 
      get_subgroup_stats()
    bmi_res_1 <- subgroup_log_reg_association(dat_non, bmi=T) %>% 
      get_subgroup_stats()
    
    dat_invas<- rbind(data_for_assoc %>% filter(case == 0),
                     data_for_assoc %>% filter(case == 1 & invas == 1))
    
    invas_cases = table(dat_invas$case)[2]
    invas_controls =  table(dat_invas$case)[1]
    
    no_bmi_res_2 <- subgroup_log_reg_association(dat_invas, bmi=F) %>% 
      get_subgroup_stats()
    bmi_res_2 <- subgroup_log_reg_association(dat_invas, bmi=T) %>% 
      get_subgroup_stats()
    
    lab1 = 
      glue("Non-invasive cancer \n({non_cases} cases/{non_controls} controls)")
    lab2 = 
      glue("Invasive cancer \n({invas_cases} cases/{invas_controls} controls)")
  }
  
  if(plot){
    cairo_pdf(paste0( results_folder, title, ".pdf"), onefile=T)
    forest(rma(yi = c(no_bmi_res_1$Estimate,
                      no_bmi_res_2$Estimate),
               sei = c(no_bmi_res_1$std_error,
                       no_bmi_res_2$std_error),
               method = "FE"),
           slab = c(lab1, lab2),
           main = "Model 1 meta-analysis")
    forest(rma(yi = c(bmi_res_1$Estimate, 
                      bmi_res_2$Estimate), 
               sei = c(bmi_res_1$std_error, 
                       bmi_res_2$std_error),
               method = "FE"),
           slab = c(lab1, lab2),
           main = "Model 2 meta-analysis")
    dev.off()
  }
  
  het_p_no_bmi <- rma(yi = c(no_bmi_res_1$Estimate, 
                             no_bmi_res_2$Estimate), 
                      sei = c(no_bmi_res_1$std_error, 
                              no_bmi_res_2$std_error),
                      method = "FE")$QEp %>% round(2)
  
  
  het_p_with_bmi <- rma(yi = c(bmi_res_1$Estimate, 
                               bmi_res_2$Estimate), 
                        sei = c(bmi_res_1$std_error, 
                                bmi_res_2$std_error),
                        method = "FE")$QEp %>% round(2)
  
  table_df <- cbind(
    bind_rows(no_bmi_res_1 %>% select(-c(Estimate, std_error)), 
              bmi_res_1 %>% select(-c(Estimate, std_error))),
    bind_rows(no_bmi_res_2 %>% select(-c(Estimate, std_error)), 
              bmi_res_2%>% select(-c(Estimate, std_error))),
    c(het_p_no_bmi, het_p_with_bmi))
  
  rownames(table_df) <- c("Model 1", "Model 2")
  names(table_df) <- c("P-value ", "OR (95% CI) ", 
                       "P-value", "OR (95% CI)",
                       "Heterogeneity p-value")
  
  table_df %<>% 
    gt(rownames_to_stub = T) %>% 
    tab_spanner(label = lab1, columns = names(table_df[1:2])) %>% 
    tab_spanner(label = lab2,columns = names(table_df[3:4])) %>%
    tab_header(title)
  if(save){
    table_df %>% 
    gtsave(filename = outfile)
  }
  else{
    table_df
  }
}

