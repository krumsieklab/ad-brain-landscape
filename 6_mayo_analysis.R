# Script to calculate metabolites associated with diagnosis in the Mayo cohort

# Generates 1 file for supplementary table 11


# libraries
library(maplet) # omics analysis 
library(bigutilsr)
library(tidyverse) # tidy
library(magrittr) # %<>%
library(glue)
library(openxlsx) # excel
library(readxl)
source("internal_functions.R") # customized functions

# input 
mayo_metabolomics_input <- 'input/metabolomics_mayo_brain_ms.xlsx'
mayo_metabolomics_input_checksum <- "e3041fcf33dd2815e85e324a1f334a0a"
mayo_metadata <- 'input/metadata_mayo_brain_ms.xlsx'
# this is result of script 4_*.R
rosmap_results <- 'results/supplementary_table_6_metabolomics_associations_with_ad.xlsx'
# output
mayo_metabolomics_output <- 'results/supplementary_table_11_mayo_ad_tcx_metabolomics_associations.xlsx'
# adjusted p-value cutoff
pcut <- 0.05 

# loading and preprocessing metabolomics data using our default protocol ----
D  <- 
  # load data and annotations
  mt_load_metabolon_v1(file=mayo_metabolomics_input, sheet= "OrigScale") %>%
  # make sure file is unchanged
  mt_load_checksum(file = mayo_metabolomics_input, checksum = mayo_metabolomics_input_checksum) %>%
  # add annotations
  mt_anno_xls(file = mayo_metadata , sheet = 1, anno_type = "samples", data_id_col = "CLIENT IDENTIFIER", anno_id_col = "SampleID_TUBE") %>%
  # create standard variable name SubjectI
  mt_anno_mutate(anno_type ='samples', col_name = 'SubjectID', term=CLIENT.IDENTIFIER) %>% # copy from projid
  # set zeros to NA
  mt_pre_zero_to_na() %>%
  # select only temporal cortex samples
  mt_modify_filter_samples(filter = BrainRegion=='TCX') %>% 
  # select only AD and control samples
  mt_modify_filter_samples(filter = Diagnosis %in%c('AD', 'Control')) %>% 
  ### preprocess data
  # filter samples with 99.9..% because the function uses a <= operator
  mt_pre_filter_missingness(samp_max = 0.9999) %>% 
  # filter metabolites with >25% missing values
  mt_pre_filter_missingness(feat_max = 0.25) %>%
  # quotient normalization
  mt_pre_norm_quot() %>%
  # log2 transform data
  mt_pre_trans_log() %>%
  # kNN imputation
  mt_pre_impute_knn() %>%
  # sample outlier removal
  mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>%
  # metabolic outlier detection followed by imputation
  mt_pre_outlier_to_na(use_quant = T, quant_thresh = 0.025) %>% 
  # kNN imputation
  mt_pre_impute_knn() 

# clean data
D %<>%
  # remove samples without sex, apoe and age info
  mt_modify_filter_samples(!is.na(Sex)) %>% 
  mt_modify_filter_samples(!is.na(AgeAtDeath)) %>% 
  # count the number of apoe4s and make an ordinal variable
  mt_anno_mutate(anno_type = 'samples', 
                 col_name = 'apoe_genotype',
                 term = case_when(APOE == '44' ~ 2, APOE == '24' ~ 1,
                                  APOE == '34' ~ 1, TRUE ~ 0))%>%
  mt_anno_mutate(anno_type = 'samples', 
                 col_name = 'any4',
                 term = case_when(APOE == '44' ~ 1, APOE == '24' ~ 1,
                                  APOE == '34' ~ 1, TRUE ~ 0))%>%
  # numeric sex
  mt_anno_mutate(anno_type = 'samples', 
                 col_name = 'Sex',
                 term = case_when(Sex == 'F' ~ 0, TRUE ~ 1))%>%
  # convert 90+ to 90
  mt_anno_mutate(anno_type = "samples", col_name = 'AgeAtDeath', 
                 term = case_when(AgeAtDeath=='90+' ~ '90', TRUE ~ AgeAtDeath)) %>%
  # convert sex to factors and age to numeric
  mt_anno_mutate(anno_type = "samples", col_name = 'Sex', 
                 term = as.factor(as.factor(Sex))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'AgeAtDeath', 
                 term = as.numeric(as.matrix(AgeAtDeath))) %>% 
  mt_anno_mutate(anno_type = "samples", col_name = 'apoe_genotype', 
                 term = as.factor(as.matrix(apoe_genotype)))

mt_write_se_xls(D, file=mayo_metabolomics_output)

# compute and save metabolic associations ----
# logistic regression of AD
D1 <- D %>% mt_modify_filter_samples(filter=Diagnosis%in%c('AD', 'Control')) %>% 
  # remove samples without apoe info
  mt_modify_filter_samples(!is.na(APOE)) %>% 
  # convert diagnosis to numeric 
  mt_anno_mutate(anno_type = "samples", col_name = 'Diagnosis', 
                 term = case_when(Diagnosis=='AD' ~ 1, Diagnosis=='Control' ~ 0, TRUE~NA_real_)) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'Diagnosis', 
                 term = as.factor(as.factor(Diagnosis)))
# untargeted association analysis
res_m1 <- association_analysis(D1, outcome = 'Diagnosis', 
                                  outcome_type =  'twofactor')
# rosmap associations
# reading rosmap results
res <- rosmap_results %>%
  excel_sheets() %>%
  purrr::set_names() %>%
  map(read_excel, path = rosmap_results) %>%
purrr::imap(~mutate(.x, group = .y)) %>% bind_rows()
# signficant metabolites in rosmap
rosmap_mets <- res %>% filter(adj_p<=0.05) %>% pull(COMP_ID) %>% unique()
# targeted mayo associations w/o confounders
res_m1_1 <- res_m1 %>% filter(COMP_ID%in%rosmap_mets) %>% 
  mutate(t_adj_p = p.adjust(p_value, method='BH'))
# targeted mayo associations w confounders
res_m2 <- association_analysis(D1, outcome = 'Diagnosis', 
                                  outcome_type =  'twofactor', 
                                  conf_formula = 'Sex + AgeAtDeath + apoe_genotype')
# finally significant metabolites in mayo
reps <- intersect(res_m1_1 %>% filter(t_adj_p<=pcut) %>% pull(name),
          res_m2 %>% filter(p_value<=pcut) %>% pull(name))
# mayo output formatting
mayo_res_ad <- left_join(res_m1_1 %>% 
                        select(name, estimate, std_error, statistic, p_value, t_adj_p) %>% 
                        dplyr::rename(estimate1=estimate, std_error1=std_error,
                         statistic1=statistic, p_value1=p_value, t_adj_p1=t_adj_p), 
                      res_m2 %>% select(name, estimate, std_error, statistic, p_value,
                                        SUPER_PATHWAY, SUB_PATHWAY, PUBCHEM, CAS, KEGG, HMDb,
                                        COMP_ID) %>%
                        dplyr::rename(estimate2=estimate, std_error2=std_error,
                                      statistic2=statistic, p_value2=p_value, 
                                      HMDB=HMDb), 
                      by='name') %>% 
  mutate(replicates= case_when(name%in%reps ~ 'Yes', TRUE ~ 'No'), 
          outcome='AD') %>%
  select(name, outcome, estimate2, std_error2, statistic2, p_value2, t_adj_p1, replicates, p_value1, estimate1, std_error1, statistic1, SUPER_PATHWAY, SUB_PATHWAY, PUBCHEM, CAS, KEGG, HMDB, COMP_ID)

# write out results
wb <- openxlsx::createWorkbook()
# creat worksheet
openxlsx::addWorksheet(wb,'AD')
# write data
openxlsx::writeData(wb, 'AD', mayo_res_ad, rowNames = F, colNames = T)
# create and add a style to the column headers
headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
# style for body
bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
# apply style
addStyle(wb, sheet = 'AD', bodyStyle, rows = 1:nrow(mayo_res_ad), cols = 1:ncol(mayo_res_ad), gridExpand = TRUE)
addStyle(wb, sheet = 'AD', headerStyle, rows = 1, cols = 1:ncol(mayo_res_ad), gridExpand = TRUE)
# write workbook
openxlsx::saveWorkbook (wb, file=mayo_metabolomics_output, overwrite=TRUE)

## finished
print("Done! Mayo metabolomics association analysis completed.") 

print("Generated excel file with supplementary table 11 in results folder!") 
