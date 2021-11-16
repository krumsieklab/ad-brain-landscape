# Customized functions used in other scripts.

# This script includes the following functions which can
# be used for any association analysis 
# (a) association_analysis
# for heterogeneity analysis of sex and apoe4 allele
# (b) beta_analysis, internal_beta_analysis
# for formating the output of ordinal regression
# (c) get_model_stats

# libraries
library(maplet) # omics analysis
library(rms) # ordinal regression

# statistical association analysis
association_analysis <- function (
  D,
  # SE
  outcome,
  # outcome to be used for response
  outcome_type,
  # outcome type - numeric/binary/ordinal
  conf_formula = NULL,
  # confounders to be corrected for
  int_w_analyte = NULL,
  # name of covariate that interacts with metabolite
  all_vals = F,
  # return all estimate and p-values from the model
  padj_method = 'BH'){# p-value adjustment method
  # merge data with sample info
  Ds <-
    D %>% maplet:::mti_format_se_samplewise () # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK
  # metabolites in data
  mets <- D %>% assay () %>% rownames ()
  # univariate analysis of numeric outcomes
  if (outcome_type == 'numeric') {
    # turn this outcome variable into numeric
    Ds %<>% mutate (!!sym(outcome) := as.numeric(as.matrix(!!sym(outcome))))
    # loop over metabolites
    univ_stats <- lapply (mets, function(x) {
      # formula for this metabolite
      # with interaction?
      if (is.null (int_w_analyte) == F) {
        this_formula <-
          as.formula(glue('{outcome} ~ {x} + {conf_formula} + {int_w_analyte}*{x}'))
        # with covariates?
      } else if (is.null(conf_formula) == F) {
        this_formula <-
          as.formula(glue('{outcome} ~ {x} + {conf_formula}'))
        # only metabolite? 
      } else{
        this_formula <- as.formula(glue('{outcome} ~ {x}'))
      }
      
      # linear regression
      this_fit <- lm(this_formula, data = Ds)
      # results summary
      this_res <-
        this_fit %>% summary() %>% coefficients() %>% data.frame()
      # format results
      names(this_res) <-
        c('estimate', 'std_error', 'statistic', 'p_value')
      # return estimate values for all covariates? 
      if (all_vals) {
        # formating output accordingly 
        tmp <- lapply(
          2:nrow(this_res),
          FUN = function(i)
            this_res[i, ]
        )
        tmp_names <-
          lapply(
            rownames(this_res)[2:nrow(this_res)],
            FUN = function(x)
              paste0(
                x,
                sep = '_',
                c('estimate', 'std_error', 'statistic', 'p_value')
              )
          )
        tmp_names[[1]] <-
          c('estimate', 'std_error', 'statistic', 'p_value')
        tmp_names[[length(tmp)]] <-
          paste0(
            paste0('analyte:', int_w_analyte),
            sep = '_',
            c('estimate', 'std_error', 'statistic', 'p_value')
          )
        this_res <-  do.call(cbind, tmp)
        names(this_res) <- unlist(tmp_names)
        # output only for the outcome and metabolite?
      } else {
        this_res <- this_res [2, ]
      }
      
      # format results
      this_res %<>% mutate(
        analyte = rownames(this_res),
        outcome = outcome,
        covariates = conf_formula
      )
      # order output columns
      this_res %<>% select(analyte, outcome, everything())
      return(this_res)
    }) %>% # create data from of results
      do.call(rbind, .) %>% data.frame() %>%
      # bind rowData
      bind_cols(D %>% rowData() %>% data.frame() %>% select(name, everything()))
    # univariate tests of binary outcomes
  } else if (outcome_type == 'twofactor') {
    # turn the outcome variable into factor
    Ds %<>% mutate(!!sym(outcome) := as.factor(as.matrix(!!sym(outcome))))
    # loop over metabolites
    univ_stats <- lapply(mets, function(x) {
      # formula for this metabolite
      # with interaction term?
      if (is.null(int_w_analyte) == F) {
        this_formula <-
          as.formula(glue('{outcome} ~ {x} + {conf_formula} + {int_w_analyte}*{x}'))
        # with covariates?
      } else if (is.null(conf_formula) == F) {
        this_formula <-
          as.formula(glue('{outcome} ~ {x} + {conf_formula}'))
        # only metabolite?
      } else{
        this_formula <- as.formula(glue('{outcome} ~ {x}'))
      }
      # logistic regression
      this_fit <-
        glm(this_formula,
            data = Ds,
            family = binomial(link = 'logit'))
      # results summary
      this_res <-
        this_fit %>% summary() %>% coefficients() %>% data.frame()
      # format results
      names(this_res) <-
        c('estimate', 'std_error', 'statistic', 'p_value')
      # return estimate values for all covariates? 
      if (all_vals) {
        # format output accordingly
        tmp <- lapply(
          2:nrow(this_res),
          FUN = function(i)
            this_res[i, ]
        )
        tmp_names <-
          lapply(
            rownames(this_res)[2:nrow(this_res)],
            FUN = function(x)
              paste0(
                x,
                sep = '_',
                c('estimate', 'std_error', 'statistic', 'p_value')
              )
          )
        tmp_names[[1]] <-
          c('estimate', 'std_error', 'statistic', 'p_value')
        tmp_names[[length(tmp)]] <-
          paste0(
            paste0('analyte:', int_w_analyte),
            sep = '_',
            c('estimate', 'std_error', 'statistic', 'p_value')
          )
        this_res <-  do.call(cbind, tmp)
        names(this_res) <- unlist(tmp_names)
        # output only outcome and metabolite
      } else {
        this_res <- this_res [2, ]
      }
      this_res %<>% mutate(
        analyte = rownames(this_res),
        outcome = outcome,
        covariates = conf_formula
      )
      # order output columns
      this_res %<>% select(analyte, outcome, everything())
      return(this_res)
    }) %>% # create data from of results
      do.call(rbind, .) %>% data.frame() %>%
      # bind rowData
      bind_cols(D %>% rowData() %>% data.frame() %>% select(name, everything()))
    # univariate tests of ordinal outcomes
  } else if (outcome_type == 'ordinal') {
    # turn the outcome variable into factor
    Ds %<>% mutate(!!sym(outcome) := as.factor(as.matrix(!!sym(outcome))))
    # loop over metabolites
    univ_stats <- lapply(mets, function(x) {
      # formula for this metabolite
      # with interaction term?
      if (is.null(int_w_analyte) == F) {
        this_formula <-
          as.formula(glue('{outcome} ~ {x} + {conf_formula} + {int_w_analyte}*{x}'))
        # with covariates?
      } else if (is.null(conf_formula) == F) {
        this_formula <-
          as.formula(glue('{outcome} ~ {x} + {conf_formula}'))
        # only metabolite?
      } else{
        this_formula <- as.formula(glue('{outcome} ~ {x}'))
      }
      # ordinal regression
      this_fit <- rms::orm(this_formula, data = Ds, family = probit)
      # results summary
      this_res <-
        get_model_stats(this_fit) %>% .$coefs %>% data.frame()
      # format results
      names(this_res) <-
        c('analyte',
          'std_error',
          'estimate',
          'statistic',
          'p_value',
          'other')
      # return estimate values for all covariates? 
      if (all_vals) {
        # format output accordingly
        starting_row <- grep(x, this_res$analyte)[1]
        tmp <-
          lapply(
            starting_row:nrow(this_res),
            FUN = function(i)
              this_res[i, ]
          )
        tmp_names <-
          lapply(
            this_res$analyte[starting_row:nrow(this_res)],
            FUN = function(x)
              paste0(
                x,
                sep = '_',
                c(
                  'analyte',
                  'std_error',
                  'estimate',
                  'statistic',
                  'p_value',
                  'other'
                )
              )
          )
        tmp_names[[1]] <-
          c('analyte',
            'std_error',
            'estimate',
            'statistic',
            'p_value',
            'other')
        tmp_names[[length(tmp)]] <-
          paste0(
            paste0('analyte:', int_w_analyte),
            sep = '_',
            c(
              'analyte',
              'std_error',
              'estimate',
              'statistic',
              'p_value',
              'other'
            )
          )
        this_res <-  do.call(cbind, tmp)
        names(this_res) <- unlist(tmp_names)
        # output only metabolite and outcome
      } else {
        this_res <- this_res %>% filter(analyte %in% x)
      }
      
      # format output
      this_res %<>% mutate(#analyte=rownames(this_res),
        outcome = outcome,
        covariates = conf_formula)
      # order outout columns
      this_res %<>% select(analyte, outcome, everything()) %>% select(-other)
      return(this_res)
      
    }) %>% # create data from of results
      do.call(rbind, .) %>% data.frame() %>%
      # bind rowdatas
      bind_cols(D %>% rowData() %>% data.frame() %>% select(name, everything()))
    
  }
  # adjust pvalues
  univ_stats %<>% mutate(adj_p = p.adjust(p_value, method = padj_method))
  # order by adjusted p-values
  univ_stats <- univ_stats[order(univ_stats$adj_p), ]
  # return results
  return(univ_stats)
}

# function called by the function beta_analysis
internal_beta_analysis <- function(
  df_grp1,
  # one of the two groups to be compared
  df_grp2,
  # one of the two groups to be compared
  group_names,
  # names of the groups (eg. male, female)
  pcut = 0.05,
  # adjusted p-value threshold for within group analysis
  pcor_method = 'BH'
  # p-value adjustment for within group analysis
) {
  # merge input tables
  this_res <-
    left_join(df_grp1,
              df_grp2,
              by = 'name',
              suffix = sprintf('_%s', group_names))
  # compute z score out of beta estimates and correspnding p-value
  tmp <- this_res %>% mutate(
    z_num = !!sym(sprintf('estimate_%s', group_names[1])) -!!sym(sprintf('estimate_%s', group_names[2])),
    z_den = sqrt((!!sym(
      sprintf('std_error_%s', group_names[1])
    )) ^ 2 + (!!sym(
      sprintf('std_error_%s', group_names[2])
    )) ^ 2),
    z_grp = z_num / z_den
  ) %>% select(name, z_grp)
  # format out data
  this_res %<>% left_join (tmp, by = 'name') %>%
    mutate(
      p_grp = 2 * pnorm(-abs(z_grp)),
      grp_adj_p = p.adjust(p_grp, method = pcor_method),
      grp_diff = case_when((
        !!sym(sprintf('adj_p_%s', group_names[1])) <= pcut |
          !!sym(sprintf('adj_p_%s', group_names[2])) <= pcut
      )  & p_grp <= pcut ~ 'Yes',
      TRUE ~ 'No')
    ) %>%
    select(name, z_grp, p_grp, grp_diff, everything()) %>%
    .[order(.$grp_diff, decreasing = T), ]
  return(this_res)
}

# comparing the beta estimates of two groups (eg. males and females)
beta_analysis <- function(
  D,
  # SE
  outcomes,
  # outcome to be used for response
  groups_variable,
  # variable to form the groups (eg. sex)
  group_names = c('grp1', 'grp2'),
  # this variable should match the sequence of levels(as.factor(colData(D)[[groups_variable]]))
  param_conf_formula = NULL
  # if any covariates to be used in the model
) {
  # create SE subsets based on groups_variable
  groups_list <-
    lapply(
      levels(as.factor(colData(D)[[groups_variable]])),
      FUN = function(x) {
        D %>% mt_modify_filter_samples(!!sym(groups_variable) == !!as.vector(x))
      }
    )
  names(groups_list) <- group_names
  # loop over groups
  res_grps <- lapply(
    group_names,
    FUN = function(grp) {
      # loop over outcomes
      tmp <- lapply(
        1:nrow(outcomes),
        FUN = function(i) {
          this_res <- association_analysis(
            D = groups_list[[grp]],
            outcome = outcomes$outcome[i],
            outcome_type = outcomes$type[i],
            conf_formula = param_conf_formula
          )
          return(this_res)
        }
      ) %>% # end outcome loop
        do.call(rbind, .)
    }
  )# end group loop
  
  # loop over outcomes for beta estimate comparison
  
  beta_list <- lapply(
    outcomes$outcome,
    FUN = function(out) {
      internal_beta_analysis(
        df_grp1 = res_grps[[1]] %>% filter(outcome == out),
        df_grp2 = res_grps[[2]] %>% filter(outcome ==
                                             out),
        group_names = group_names
      )
    }
  )
  return(beta_list)
}
# parser from stackoverflow by CoderGuy123
get_model_stats <- function(x, precision = 60) {
  # remember old number formatting function
  # (which would round and transforms p-values to formats like "<0.01")
  old_format_np <- rms::formatNP
  # substitute it with a function which will print out as many digits as we want
  assignInNamespace("formatNP", function(x, ...)
    formatC(x, format = "f", digits = precision), "rms")
  
  # remember old width setting
  old_width <- options('width')$width
  # substitute it with a setting making sure the table will not wrap
  options(width = old_width + 4 * precision)
  
  # actually print the data and capture it
  cap <- capture.output(print(x))
  
  # restore original settings
  options(width = old_width)
  assignInNamespace("formatNP", old_format_np, "rms")
  
  # model stats
  stats <- c()
  stats$R2.adj <-
    str_match(cap, "R2 adj\\s+ (\\d\\.\\d+)") %>% na.omit() %>% .[, 2] %>% as.numeric()
  
  # coef stats lines
  coef_lines <-
    cap[which(str_detect(cap, "Coef\\s+S\\.E\\.")):(length(cap) - 1)]
  
  # parse
  coef_lines_table <-
    
    suppressWarnings(readr::read_table(coef_lines %>% stringr::str_c(collapse = "\n")))
  colnames(coef_lines_table)[1] <- "Predictor"
  
  list(stats = stats,
       coefs = coef_lines_table)
}