# Script to preprocess ROS/MAP metabolomics data. Generates a file with
# preprocesed data and 2 files: Supplementary table 2 and 2 sheets of
# supplementary table 3, last sheet of supplementary table was manually assembled from ROS/MAP codebooks.

# Note: Medication correction at line 98 is time intensive (~2.5h on a single core)


# input files
metabolomics_input <- 'input/metabolomics_rosmap514_brain_ms.xlsx'
metabolomics_input_checksum <- "2de2f322785c32888eacf3967b79b65e"
metabolomics_metadata <- 'input/metadata_rosmap514_brain_ms.xlsx'
medication_columns <- 'input/medication_columns.xlsx'
# output files
met_pdata <- 'results/tmp_rosmap_metabolomics_processed_medcor_data.xlsx'
metabolites_analyzed <- 'results/supplementary_table_2_metabolites_analyzed.xlsx'
medication_output <- 'results/supplementary_table_3_medication_effect_metabolites.xlsx'

# libraries
library(maplet) # omics analysis 
library(tidyverse) # tidy
library(openxlsx) # excel things

# data loading ----
D  <- 
  # load data and annotations
  mt_load_metabolon_v1(file = metabolomics_input, sheet = "OrigScale") %>%
  # make sure file is unchanged
  mt_load_checksum(file = metabolomics_input, checksum = metabolomics_input_checksum) %>%
  # create standard variable name SubjectI
  mt_anno_mutate(anno_type = 'samples', col_name = 'SubjectID', term = CLIENT.IDENTIFIER) %>% # copy from projid
  # set zeros to NA
  mt_pre_zero_to_na()

# identify and remove the qc samples in the dataset
D %<>% mt_modify_filter_samples(filter = !SubjectID %in% colData(D)$SubjectID[grep("NIST", colData(D)$SubjectID)]) %>%
  mt_modify_filter_samples(filter = !SubjectID %in% colData(D)$SubjectID[grep("MTRX", colData(D)$SubjectID)]) 

# preprocessing ----
D %<>% 
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

# attach and clean metadata ---- 
D %<>% # attach metadata
  mt_anno_xls(file=metabolomics_metadata, sheet=1, 
              anno_type = 'samples', anno_id_col = 'projid', data_id_col = 'SubjectID') %>%
  # dichotomize reagan score
  mt_anno_mutate(anno_type = 'samples',
                 col_name = 'niareagansc',
                 term = case_when(niareagansc >=3 ~ 0, niareagansc < 3 ~ 1, TRUE ~ NA_real_))%>%
  # count the number of apoe4s and make an ordinal variable
  mt_anno_mutate(anno_type = 'samples',
                 col_name = 'apoe_genotype',
                 term = case_when(apoe_genotype == 44 ~ 2, apoe_genotype == 24 ~ 1,
                                  apoe_genotype == 34 ~ 1, is.na(apoe_genotype) ~ NA_real_, TRUE ~ 0))%>%
  # turn fasting nas and 9s into 3s i.e don't know status
  mt_anno_mutate(anno_type = 'samples',
                 col_name = 'fasting',
                 term = case_when(fasting == 1 ~ 1, fasting == 2 ~ 2,TRUE ~ 3)) %>%
  # turn the category 6 into NAs; 6 represents non-AD dementia
  mt_anno_mutate(anno_type = 'samples',
                 col_name = 'cogdx',
                 term = case_when(cogdx == 6 ~ NA_real_, TRUE ~ cogdx)) %>%
  # turn the category 6 into NAs; 6 represents non-AD dementia
  mt_anno_mutate(anno_type = 'samples',
                 col_name = 'dcfdx',
                 term = case_when(dcfdx == 6 ~ NA_real_, TRUE ~ dcfdx)) %>%
  # create a diagnosis based on Mayo recommendations
  # AD cases: Braak score ≥ 4, CERAD score ≤ 2.
  # Controls: Braak score ≤ 3, CERAD score ≥ 3.
  mt_anno_mutate(anno_type = 'samples',
                 col_name = 'diagnosis',
                 term = case_when((braaksc >= 4 & ceradsc <= 2) ~ 1,
                                  (braaksc <= 3 & ceradsc >= 3) ~ 0, TRUE ~ NA_real_)) %>%
  # converting amyloid load to its square root value
  mt_anno_mutate(anno_type = 'samples',
                 col_name = 'sqrt_amyloid',
                 term = case_when(!is.na(amyloid) ~ sqrt(amyloid), is.na(amyloid) ~ NA_real_)) %>%
  # converting cogng_random_slope - cognitive decline to positively correlate with AD pathology
  mt_anno_mutate(anno_type = 'samples',
                 col_name = 'cogng_random_slope',
                 term = -1*cogng_random_slope) %>%
  mt_reporting_data()

# medication correction ----
Dmc <- D %>% # converting all medication columns to factors
  mt_anno_class(anno_type = 'samples', file=medication_columns, sheet=1) %>%
  # medication correction
  mt_pre_confounding_correction_stepaic(cols_to_correct = names(colData(D))[grep("_rx", names(colData(D)))], cols_to_exclude = c("ad_rx", "neurologic_rx"))

# quantifying effects of medications left out in previous step
Dmc2 <- Dmc %>% 
  mt_pre_confounding_correction_stepaic(cols_to_correct = c("ad_rx", "neurologic_rx"))

# delete samples with missing confounders to be used in association analysis
Dmc %<>% mt_modify_filter_samples(!is.na(bmi)) %>% 
  mt_modify_filter_samples(!is.na(msex)) %>% 
  mt_modify_filter_samples(!is.na(educ)) %>%
  mt_modify_filter_samples(!is.na(apoe_genotype)) %>% 
  mt_modify_filter_samples(!is.na(age_death)) %>%
  mt_modify_filter_samples(!is.na(pmi)) %>% 
  mt_reporting_data()

# save outputs ---- 
# saving medication corrected data for following analysis
mt_write_se_xls(Dmc, file=met_pdata)

# metabolites analyzed
mets <- Dmc %>% rowData() %>% data.frame() %>%
  select(-PATHWAY_SORTORDER, - BIOCHEMICAL) %>% dplyr::rename(HMDB=ID) %>%
  select(name, everything())

# write out metabolites i.e supplementary table 2
out <- 'metabolites_analyzed'; this_res <- mets
wb <- openxlsx::createWorkbook()
# creat worksheet
openxlsx::addWorksheet(wb,sprintf('%s', out))
# write data
openxlsx::writeData(wb, sprintf('%s', out), this_res, rowNames = F, colNames = T)
# create and add a style to the column headers
headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
# style for body
bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
# apply style
addStyle(wb, sheet = sprintf('%s', out), bodyStyle, rows = 1:nrow(this_res), cols = 1:ncol(this_res), gridExpand = TRUE)
addStyle(wb, sheet = sprintf('%s', out), headerStyle, rows = 1, cols = 1:ncol(this_res), gridExpand = TRUE)
# write out
openxlsx::saveWorkbook (wb, file=metabolites_analyzed, overwrite=TRUE)

# medication effect
df_medcor <- metadata(Dmc)$results[[grep('pre_confounding_correction_stepaic', names(metadata(Dmc)$results))]]$output %>% 
  data.frame() %>% bind_cols(rowData(Dmc) %>% data.frame() %>% select(name)) %>%
  select(-feature) %>% select(name, covariates, model.rsq, model.pvalue)%>%
  dplyr::rename(medications=covariates,
                model_rsq = model.rsq,
                model_pval = model.pvalue)

ad_medcor <- metadata(Dmc2)$results[[grep('pre_confounding_correction_stepaic', names(metadata(Dmc2)$results))[2]]]$output %>% 
  data.frame() %>% bind_cols(rowData(Dmc2) %>% data.frame() %>% select(name)) %>% 
  select(-feature) %>% select(name, covariates, model.rsq, model.pvalue)%>%
  dplyr::rename(medications=covariates,
                model_rsq = model.rsq,
                model_pval = model.pvalue)

# write out medication effects i.e supplementary table 3
out <- 'medication_correction'; this_res <- df_medcor
wb <- openxlsx::createWorkbook()
# creat worksheet
openxlsx::addWorksheet(wb,sprintf('%s', out))
# write data
openxlsx::writeData(wb, sprintf('%s', out), this_res, rowNames = F, colNames = T)
# create and add a style to the column headers
headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
# style for body
bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
# apply style
addStyle(wb, sheet = sprintf('%s', out), bodyStyle, rows = 1:nrow(this_res), cols = 1:ncol(this_res), gridExpand = TRUE)
addStyle(wb, sheet = sprintf('%s', out), headerStyle, rows = 1, cols = 1:ncol(this_res), gridExpand = TRUE)

out <- 'ad_correction_effect'; this_res <- ad_medcor
# creat worksheet
openxlsx::addWorksheet(wb,sprintf('%s', out))
# write data

openxlsx::writeData(wb, sprintf('%s', out), this_res, rowNames = F, colNames = T)
# create and add a style to the column headers
headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
# style for body
bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
# apply style
addStyle(wb, sheet = sprintf('%s', out), bodyStyle, rows = 1:nrow(this_res), cols = 1:ncol(this_res), gridExpand = TRUE)
addStyle(wb, sheet = sprintf('%s', out), headerStyle, rows = 1, cols = 1:ncol(this_res), gridExpand = TRUE)

# write out
openxlsx::saveWorkbook (wb, file=medication_output, overwrite=TRUE)

## finished
print("Done! preprocessing ROS/MAP metabolomics data completed.") 
print("Generated excel files with supplementary tables 2 and 3 in results folder!") 