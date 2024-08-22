# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#                  Location-Scale Models in Psychotherapy Research:         #
#                       Predictors of Heterogeneity                         #
#                         Sensitivity Analyses                              #
#                                                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Note: The terms "predictor" and "scale moderator" are used interchangeably 

# for a detailed documentation of the dataset see: 
# https://docs.metapsy.org/databases/depression-psyctr/

# load dependencies
pacman::p_load(
  metapsyData,
  dplyr,
  metafor,
  skimr,
  sjPlot,
  metapsyTools,
  purrr,
  tidyr,
  forcats,
  readxl)


# load prepared data with aggregation of effect sizes on an arm-level
load("data/data_agg.rda")
load("results/res_mod.rda")
data <- data_agg

# load study names for that the format got aggregated
load("results/comparisons_aggregate_format.rda")

# source function to fit univariate moderations
source("utils/fit.models.R")



# 1 RoB sensitivity analysis -----------------------------------------
# without adjusting for RoB

# preparation
data <- within(data, {
  condition_arm1 = fct_relevel(condition_arm1, "cbt", "3rd","bat", "dyn", "ipt", 
                               "lrt", "other psy", "pst", "sup")
  format= fct_relevel(format, "ind", "grp", "gsh", "tel", "oth")
  target_group = fct_relevel(target_group, "adul", "child", "adol", "stud", 
                             "old", "med", "ppd", "oth")
  country= fct_relevel(country, "eu", "au", "can", "eas", "uk", "us", "oth")
  age_group = fct_relevel(age_group, "adul", "child", "adol", "yadul", "old")
  diagnosis = fct_relevel(diagnosis, "mdd", "sub", "mood", "cut", "chr")
})

## 1.1 analysis ------------------------------------------------------------

data$year100 <- data$year/100
data$totaln_bl10 <- data$totaln_bl/10
data$instrument_red <- as.factor(data$instrument_red)

# create a reduced dataset for the moderation analysis of "format" without 
# studies with aggregated values
data_format <- subset(data, !(comparison_id %in% comparisons_aggregate_format))

# drop studies with parent-reported outcomes as n= 2
data_rating <- subset(data, data$rating !='parent_report')
data_rating$rating <- data_rating$rating %>% droplevels()


# run univariate moderation models
res_mod_sens_rob <- tibble(
  model_name= c(
    # Intervention and control characteristics
    "condition_arm1", "condition_arm2", "n_sessions_arm1", "format",
    
    # study characteristics
    "country", "rating",  "rob",  "year100", "totaln_bl10", "target_group",
    "instrument_red",
    
    # participant characteristics
    "age_group", "comorbid_mental", "diagnosis", "percent_women", "recruitment"
  ),
  moderators = c(
    # Intervention characteristics
    ~condition_arm1,
    ~condition_arm2,              # control group condition
    ~n_sessions_arm1,             # IG: Average number of sessions received
    ~format,                      # therapy format (ind= individual; grp= group; gsh= guided self-help; tel= telephone; cpl= couple therapy; oth= other (mixed formats); ush= unguided self-help
    
    
    # Study characteristics
    ~country,                     # us= USA; uk= United Kingdom; eu= Europe; can= Canada; au= australia; eas= east asia; oth= other
    ~rating,                      # self-reported (\"self-report\") or clinician-rated (\"clinician\")
    ~rob,                         # overall risk of bias: 0 (high risk)-4 (low risk))
    ~year100,                     # year of publication
    ~totaln_bl10,                 # total N across both groups at baseline
    ~target_group,                # adul= adults, old= older adults, stud= student population, ppd= women with perinatal depression; med= comorbid medical disorder; oth= other
    ~instrument_red,
    
    # Participant characteristics (study-level)
    ~age_group,                   # adul= adults, 5= older adults (≥55 years); old= older old adults (≥75 years)
    ~comorbid_mental,             # comrbid mental disorder at baseline (yes/no))
    ~diagnosis,                   # mdd= major depression; mood= mood disorder; cut= cut-off score; sub= subclinical depression;chr= chronic depression
    ~percent_women,               # % of women at baseline
    ~recruitment                  # com= community; clin= clinical; oth= other
  ),
  
  database= list(data, data, data, data_format, data, data_rating, data, data, 
                 data, data, data, data, data, data, data, data),
  model= pmap(list(moderators, model_name, database), fitModels),
  
  n_missing = data %>% 
    select(c(
      "condition_arm1", "condition_arm2", "n_sessions_arm1", "format",
      "country", "rating",  "rob",  "year100", "totaln_bl10", "target_group", 
      "instrument_red", "age_group", "comorbid_mental", "diagnosis", "percent_women", 
      "recruitment"
    )) %>% 
    map_dbl(.,~sum(is.na(.x)))
)


# inspect results
res_mod_sens_rob$model 

# name models
names(res_mod_sens_rob$model) <- res_mod_sens_rob$model_name


## 1.2 LR tests ------------------------------------------------------------
# conduct a likelihood ratio test comparing a model that allows τ2 to differ 
# across subgroups with a model assuming homoscedastic τ2

reduced_model <- rma(yi=.g, sei=.g_se, scale= ~1, data=data, method= "REML",
                     test="knha", control= list(optimizer="BFGS"))

reduced_model_format <- rma(yi=.g, sei=.g_se, scale= ~1, data=data_format, 
                            method= "REML", test="knha", control= list(optimizer="BFGS"))

reduced_model_rating <- rma(yi=.g, sei=.g_se, scale= ~1, data=data_rating, 
                            method= "REML", test="knha", control= list(optimizer="BFGS"))

reduced_model_session <- data %>% 
  select(.g, .g_se, n_sessions_arm1) %>% drop_na() %>% 
  rma(yi=.g, sei=.g_se, scale= ~1, data=., method= "REML", test="knha", 
      control= list(optimizer="BFGS"))

reduced_model_baselineN <- data %>% 
  select(.g, .g_se, totaln_bl10) %>% drop_na() %>% 
  rma(yi=.g, sei=.g_se, scale= ~1, data=., method= "REML", test="knha", 
      control= list(optimizer="BFGS"))

reduced_model_percentWomen <- data %>% 
  select(.g, .g_se, percent_women) %>% drop_na() %>% 
  rma(yi=.g, sei=.g_se, scale= ~1, data=., method= "REML", test="knha", 
      control= list(optimizer="BFGS"))


# likelihood ratio tests (separately for format and predictors with missing values)
res_likelihoodRatio <- list(
  n_sessions_arm1 = anova(res_mod_sens_rob$model[c("n_sessions_arm1")][[1]], reduced_model_session),
  format = anova(res_mod_sens_rob$model[c("format")][[1]], reduced_model_format),
  rating = anova(res_mod$model[c("rating")][[1]], reduced_model_rating),
  totaln_bl10 = anova(res_mod_sens_rob$model[c("totaln_bl10")][[1]], reduced_model_baselineN),
  percent_women = anova(res_mod_sens_rob$model[c("percent_women")][[1]], reduced_model_percentWomen),
  map(res_mod_sens_rob$model[c( "condition_arm1", "condition_arm2", "country", 
                       "year100", "target_group","instrument_red", "age_group", "rob",
                       "comorbid_mental", "diagnosis",  "recruitment")], 
      ~anova(.x,reduced_model))
) 

# order the results as in res_mod_sens_rob and save in results list
res_mod_sens_rob$res_likelihoodRatio <- c(res_likelihoodRatio[-length(res_likelihoodRatio)], 
                                 res_likelihoodRatio[[6]])%>% 
  .[match(res_mod_sens_rob$model_name, names(.))]


res_mod_sens_rob$res_likelihoodRatio


## 1.3 profile likelihood CIs ----------------------------------------------

# for all non-fixed scale parameters
res_mod_sens_rob$profileLikelihoodCI <- map(res_mod_sens_rob$model_name, 
                                            ~confint(res_mod_sens_rob$model[c(.x)][[1]], 
                                            digits= 3))

profileLikelihoodCI= res_mod_sens_rob$profileLikelihoodCI

#save(profileLikelihoodCI, file="results/res_mod_sens_rob_profile_likelihood.rda")
#res_mod_sens_rob$profileLikelihoodCI  = profileLikelihoodCI

## 1.4 predicted tau^2 -----------------------------------------------------


# function to iterate over predictor levels for the calculation 
newScales <- function(moderator) {
  # determine number of levels
  levels_count <- length(levels(moderator))-1
  # create Diagonal Matrix
  newscales <- diag(levels_count)
  newscales_list <- split(newscales, row(newscales))
  newscales_list_reordered <- c(list(rep(0, levels_count)), newscales_list)
  
  return(newscales_list_reordered)
}

# Iterate over newscales_list and compute predictions (all for studies with low RoB)
res_mod_sens_rob$pred_tau2 <- list(
  condition_arm1= map(newScales(data$condition_arm1),
                      ~predict(res_mod_sens_rob$model[c("condition_arm1")][[1]], 
                               newscale = .x, digits = 3, transf = exp)),
  condition_arm2= map(newScales(data$condition_arm2), 
                      ~predict(res_mod_sens_rob$model[c("condition_arm2")][[1]], 
                               newscale = .x, digits = 3, transf = exp)),
  n_sessions_arm1= list(predict(res_mod_sens_rob$model[c("n_sessions_arm1")][[1]], 
                                newscale= c(0), 
                                digits = 3, transf = exp),
                        predict(res_mod_sens_rob$model[c("n_sessions_arm1")][[1]], 
                           newscale= c(mean(data$n_sessions_arm1, na.rm=TRUE)), 
                           digits = 3, transf = exp)),
  format= map(newScales(data$format), 
              ~predict(res_mod_sens_rob$model[c("format")][[1]],
                       newscale = .x, digits = 3, transf = exp)),
  country= map(newScales(data$country), 
               ~predict(res_mod_sens_rob$model[c("country")][[1]], 
                        newscale = .x, digits = 3, transf = exp)),
  rating= map(newScales(data_rating$rating), 
              ~predict(res_mod_sens_rob$model[c("rating")][[1]], 
                       newscale = .x, digits = 3, transf = exp)),
  rob= list(predict(res_mod_sens_rob$model[c("rob")][[1]], 
               newscale= c(0), digits = 3, transf = exp),
            predict(res_mod_sens_rob$model[c("rob")][[1]], 
                    newscale= c(mean(data$rob,  na.rm=TRUE)), digits = 3, transf = exp)),
  year100= list(predict(res_mod_sens_rob$model[c("year100")][[1]], 
                   newscale= c(0), digits = 3, transf = exp),
                predict(res_mod_sens_rob$model[c("year100")][[1]], 
                        newscale= c(mean(data$year100, na.rm=TRUE)), digits = 3, transf = exp)),
  totaln_bl10= list(predict(res_mod_sens_rob$model[c("totaln_bl10")][[1]], 
                       newscale= c(0), 
                       digits = 3, transf = exp),
                    predict(res_mod_sens_rob$model[c("totaln_bl10")][[1]], 
                            newscale= c(mean(data$totaln_bl10,  na.rm=TRUE)), 
                            digits = 3, transf = exp)),
  target_group= map(newScales(data$target_group), 
                    ~predict(res_mod_sens_rob$model[c("target_group")][[1]],
                             newscale = .x, digits = 3, transf = exp)),
 instrument_red= map(newScales(data$instrument_red), 
                      ~predict(res_mod_sens_rob$model[c("instrument_red")][[1]], 
                               newscale = .x, digits = 3, transf = exp)),
  age_group= map(newScales(data$age_group), 
                 ~predict(res_mod_sens_rob$model[c("age_group")][[1]],  
                          newscale = .x, digits = 3, transf = exp)),
  comorbid_mental= map(newScales(data$comorbid_mental), 
                       ~predict(res_mod_sens_rob$model[c("comorbid_mental")][[1]], 
                                newscale = .x, digits = 3, transf = exp)),
  diagnosis= map(newScales(data$diagnosis), 
                 ~predict(res_mod_sens_rob$model[c("diagnosis")][[1]],   
                          newscale = .x, digits = 3, transf = exp)),
  percent_women= list(predict(res_mod_sens_rob$model[c("percent_women")][[1]], 
                         newscale= c(0), 
                         digits = 3, transf = exp),
                      predict(res_mod_sens_rob$model[c("percent_women")][[1]], 
                              newscale= c(mean(data$percent_women,  na.rm=TRUE)), 
                              digits = 3, transf = exp)),
  recruitment= map(newScales(data$recruitment), 
                   ~predict(res_mod_sens_rob$model[c("recruitment")][[1]], 
                            newscale = .x, digits = 3, transf = exp))
)


## 1.5 save ----------------------------------------------------------------

# location estimates + omnibus test
tibble(
  moderator = res_mod_sens_rob$model_name,
  N = map_dbl(res_mod_sens_rob$model, "k"),
  loc_coef = map_dbl(res_mod_sens_rob$model, "b")%>% round(3),
  loc_coef_CI_lower = map_dbl(res_mod_sens_rob$model, "ci.lb")%>% round(3),
  loc_coef_CI_upper = map_dbl(res_mod_sens_rob$model, "ci.ub")%>% round(3),
  pval_omnibus_test_scale_coefs = map_dbl(res_mod_sens_rob$model, "QSp")%>% round(3),
  pval_likelihoodRatio = unlist(map(res_mod_sens_rob$res_likelihoodRatio, ~round(.x$pval, digits = 3)))
)-> res_sensitivity_rob_location

res_sensitivity_rob_location


# scale estimates
tibble(
  moderator = rep(res_mod_sens_rob$model_name,
                  times= map(res_mod_sens_rob$model, "q") %>% unlist),
  coefficients = unlist(map(res_mod_sens_rob$model, ~dimnames(.x$Z)[[2]])),
  scale_coef = unlist(map(res_mod_sens_rob$model, "alpha"))%>% round(3),
  scale_coef_CI_lower= unlist(map(res_mod_sens_rob$model, "ci.lb.alpha"))%>% round(3),
  scale_coef_CI_lower_profileLikelihod = map(res_mod_sens_rob$profileLikelihoodCI,
                                           ~map_chr(.x, ~unlist(.x)["random2"]) %>% 
                                             .[-length(.)] %>% as.numeric()) %>% 
  unlist() %>% round(3),
  scale_coef_CI_upper= unlist(map(res_mod_sens_rob$model, "ci.ub.alpha"))%>% round(3),
  scale_coef_CI_upper_profileLikelihod = map(res_mod_sens_rob$profileLikelihoodCI,
                                             ~map_chr(.x, ~unlist(.x)["random3"]) %>% 
                                               .[-length(.)] %>% as.numeric()) %>% 
    unlist() %>% round(3),
  scale_coef_se = unlist(map(res_mod_sens_rob$model, "se.alpha"))%>% round(3),
  tau2 = res_mod_sens_rob$pred_tau2 %>% unlist %>% .[grep("pred", names(.))] %>% as.numeric() %>% na.omit %>% round(3),
  tau2_CI_lower =res_mod_sens_rob$pred_tau2 %>% unlist %>% .[grep("ci.lb", names(.))] %>% 
    as.numeric() %>% na.omit %>% round(3),
  tau2_CI_upper = res_mod_sens_rob$pred_tau2 %>% unlist %>% .[grep("ci.ub", names(.))] %>% 
    as.numeric() %>% na.omit %>% round(3),
  F_t_z_val = unlist(map(res_mod_sens_rob$model, "zval.alpha"))%>% round(3),
  pval = unlist(map(res_mod_sens_rob$model, "pval.alpha")) %>% round(3)
) -> res_sensitivity_rob_scale



# save in excel sheet
write.xlsx(res_sensitivity_rob_location, file = "results/results.xlsx",
           sheetName = "RoB sensitivity analysis (location +omnibus)", append = TRUE)
# add a second data set in a new worksheet
write.xlsx(res_sensitivity_rob_scale, file = "results/results.xlsx", 
           sheetName="RoB sensitivity analysis (scale)", append=TRUE)

# save as rda file
save(res_mod_sens_rob, file="results/res_mod_sens_rob.rda")




# 2 Sample size sensitivity analysis ---------------------------------------------
# with adjustment for baseline sample size 


# preparation
data <- within(data, {
  condition_arm1 = fct_relevel(condition_arm1, "cbt", "3rd","bat", "dyn", "ipt", "lrt", "other psy", "pst", "sup")
  format= fct_relevel(format, "ind", "grp", "gsh", "tel", "oth")
  target_group = fct_relevel(target_group, "adul", "child", "adol", "stud", "old", "med", "ppd", "oth")
  country= fct_relevel(country, "eu", "au", "can", "eas", "uk", "us", "oth")
  age_group = fct_relevel(age_group, "adul", "child", "adol", "yadul", "old")
  diagnosis = fct_relevel(diagnosis, "mdd", "sub", "mood", "cut", "chr")
})

## 2.1 analysis ------------------------------------------------------------

data$year100 <- data$year/100
data$totaln_bl10 <- data$totaln_bl/10
data$instrument_red <- as.factor(data$instrument_red)

# create a reduced dataset for the moderation analysis of "format" without 
# studies with aggregated values
data_format <- subset(data, !(comparison_id %in% comparisons_aggregate_format))

# drop studies with parent-reported outcomes as n= 2
data_rating <- subset(data, data$rating !='parent_report')
data_rating$rating <- data_rating$rating %>% droplevels()

# run univariate moderation models
res_mod_sens_n <- tibble(
  model_name= c(
    # Intervention and control characteristics
    "condition_arm1", "condition_arm2", "n_sessions_arm1", "format",
    
    # study characteristics
    "country", "rating",  "rob",  "year100", "totaln_bl10", "target_group",
    "instrument_red",
    
    # participant characteristics
    "age_group", "comorbid_mental", "diagnosis", "percent_women", "recruitment"
  ),
  moderators = c(
    # Intervention characteristics
    ~condition_arm1+ totaln_bl10,
    ~condition_arm2+totaln_bl10,              # control group condition
    ~n_sessions_arm1+totaln_bl10,             # IG: Average number of sessions received
    ~format+totaln_bl10,                      # therapy format (ind= individual; grp= group; gsh= guided self-help; tel= telephone; cpl= couple therapy; oth= other (mixed formats); ush= unguided self-help
    
    
    # Study characteristics
    ~country+totaln_bl10,                     # us= USA; uk= United Kingdom; eu= Europe; can= Canada; au= australia; eas= east asia; oth= other
    ~rating+totaln_bl10,                      # self-reported (\"self-report\") or clinician-rated (\"clinician\")
    ~rob+totaln_bl10,                         # overall risk of bias: 0 (high risk)-4 (low risk))
    ~year100+totaln_bl10,                     # year of publication
    ~totaln_bl10,                             # total N across both groups at baseline
    ~target_group+totaln_bl10,                # adul= adults, old= older adults, stud= student population, ppd= women with perinatal depression; med= comorbid medical disorder; oth= other
    ~instrument_red+totaln_bl10,
    # Participant characteristics (study-level)
    ~age_group+totaln_bl10,                   # adul= adults, 5= older adults (≥55 years); old= older old adults (≥75 years)
    ~comorbid_mental+totaln_bl10,             # comrbid mental disorder at baseline (yes/no))
    ~diagnosis+totaln_bl10,               # mdd= major depression; mood= mood disorder; cut= cut-off score; sub= subclinical depression;chr= chronic depression
    ~percent_women+totaln_bl10,               # % of women at baseline
    ~recruitment+totaln_bl10              # com= community; clin= clinical; oth= other
  ),
  
  database= list(data, data, data, data_format, data, data_rating, data, data, 
                 data, data, data, data, data, data, data, data),
  model= pmap(list(moderators, model_name, database), fitModels),
  
  n_missing = data %>% 
    select(c(
      "condition_arm1", "condition_arm2", "n_sessions_arm1", "format",
      "country", "rating",  "rob",  "year100", "totaln_bl10", "target_group",
      "instrument_red", "age_group", "comorbid_mental", "diagnosis", "percent_women", 
      "recruitment"
    )) %>% 
    map_dbl(.,~sum(is.na(.x)))
)


# inspect results
res_mod_sens_n$model 

# name models
names(res_mod_sens_n$model) <- res_mod_sens_n$model_name


## 2.2 LR tests ------------------------------------------------------------
# conduct a likelihood ratio test comparing a model that allows τ2 to differ 
# across subgroups with a model assuming homoscedastic τ2

reduced_model <- rma(yi=.g, sei=.g_se, scale= ~totaln_bl10, data=data, method= "REML",
                     test="knha", control= list(optimizer="BFGS"))

reduced_model_format <- rma(yi=.g, sei=.g_se, scale= ~totaln_bl10, data=data_format, 
                            method= "REML", test="knha", control= list(optimizer="BFGS"))

reduced_model_rating <- rma(yi=.g, sei=.g_se, scale= ~totaln_bl10, data=data_rating, 
                            method= "REML", test="knha", control= list(optimizer="BFGS"))

reduced_model_session <- data %>% 
  select(.g, .g_se,totaln_bl10, n_sessions_arm1) %>% drop_na() %>% 
  rma(yi=.g, sei=.g_se, scale= ~totaln_bl10, data=., method= "REML", test="knha", 
      control= list(optimizer="BFGS"))

reduced_model_baselineN <- data %>% 
  select(.g, .g_se, totaln_bl10) %>% drop_na() %>% 
  rma(yi=.g, sei=.g_se, scale= ~1, data=., method= "REML", test="knha", 
      control= list(optimizer="BFGS"))

reduced_model_percentWomen <- data %>% 
  select(.g, .g_se, percent_women, totaln_bl10) %>% drop_na() %>% 
  rma(yi=.g, sei=.g_se, scale= ~totaln_bl10, data=., method= "REML", test="knha", 
      control= list(optimizer="BFGS"))


# likelihood ratio tests (separately for format and predictors with missing values)
res_likelihoodRatio <- list(
  n_sessions_arm1 = anova(res_mod_sens_n$model[c("n_sessions_arm1")][[1]], reduced_model_session),
  format = anova(res_mod_sens_n$model[c("format")][[1]], reduced_model_format),
  rating = anova(res_mod_sens_n$model[c("rating")][[1]], reduced_model_rating),
  totaln_bl10 = anova(res_mod_sens_n$model[c("totaln_bl10")][[1]], reduced_model_baselineN),
  percent_women = anova(res_mod_sens_n$model[c("percent_women")][[1]], reduced_model_percentWomen),
  map(res_mod_sens_n$model[c( "condition_arm1", "condition_arm2", "country",   
                               "year100", "target_group","instrument_red", 
                              "age_group", "rob", "comorbid_mental", "diagnosis",  
                              "recruitment")], 
      ~anova(.x,reduced_model))
) 

# order the results as in res_mod_sens_se and save in results list
res_mod_sens_n$res_likelihoodRatio <- c(res_likelihoodRatio[-length(res_likelihoodRatio)], 
                                         res_likelihoodRatio[[6]])%>% 
  .[match(res_mod_sens_n$model_name, names(.))]


res_mod_sens_n$res_likelihoodRatio


## 1.3 profile likelihood CIs ----------------------------------------------

# for all non-fixed scale parameters
res_mod_sens_n$profileLikelihoodCI <- map(res_mod_sens_n$model_name, 
                                           ~confint(res_mod_sens_n$model[c(.x)][[1]], 
                                                    digits= 3))

profileLikelihoodCI = res_mod_sens_n$profileLikelihoodCI 
save(profileLikelihoodCI, file= "results/res_sens_n_profile_likelihood.rda")

#res_mod_sens_n$profileLikelihoodCI = profileLikelihoodCI

## 2.4 predicted tau^2 -----------------------------------------------------


# function to iterate over predictor levels for the calculation 
newScales <- function(moderator) {
  # determine number of levels
  levels_count <- length(levels(moderator))
  # create Diagonal Matrix
  newscales <- diag(levels_count)
  # set last cell to mean N across studies 
  newscales[, levels_count] <- mean(data$totaln_bl10, na.rm=T)
  # convert the matrix to a list of vectors
  newscales_list <- split(newscales, row(newscales))
  
  # reorder the list of vectors, moving the last row to the first position 
  newscales_list_reordered <- c(newscales_list[length(newscales_list)], newscales_list[-length(newscales_list)])
  
  # add a last row to calulculate the PI for sample size
  newscales_list_reordered <-c(newscales_list_reordered, newscales_list_reordered[1])
  
  return(newscales_list_reordered)
}


# Iterate over newscales_list and compute predictions (all for studies with low RoB)
res_mod_sens_n$pred_tau2 <- list(
  condition_arm1= map(newScales(data$condition_arm1),
                      ~predict(res_mod_sens_n$model[c("condition_arm1")][[1]], 
                               newscale = .x, digits = 3, transf = exp)),
  condition_arm2= map(newScales(data$condition_arm2), 
                      ~predict(res_mod_sens_n$model[c("condition_arm2")][[1]], 
                               newscale = .x, digits = 3, transf = exp)),
  n_sessions_arm1= list(
    predict(res_mod_sens_n$model[c("n_sessions_arm1")][[1]], 
            newscale= c(0,0), digits = 3, transf = exp),
    predict(res_mod_sens_n$model[c("n_sessions_arm1")][[1]], 
                           newscale= c(mean(data$n_sessions_arm1, na.rm=TRUE),
                                       mean(data$totaln_bl10, na.rm=T)), digits = 3, transf = exp),
    predict(res_mod_sens_n$model[c("n_sessions_arm1")][[1]], 
            newscale= c(0,
                        mean(data$totaln_bl10, na.rm=T)), digits = 3, transf = exp)),
  format= map(newScales(data$format), 
              ~predict(res_mod_sens_n$model[c("format")][[1]], 
                       newscale = .x, digits = 3, transf = exp)),
  country= map(newScales(data$country), 
               ~predict(res_mod_sens_n$model[c("country")][[1]], 
                        newscale = .x, digits = 3, transf = exp)),
  rating= map(newScales(data_rating$rating), 
              ~predict(res_mod_sens_n$model[c("rating")][[1]], 
                       newscale = .x, digits = 3, transf = exp)),
  rob= list(
    predict(res_mod_sens_n$model[c("rob")][[1]], 
            newscale= c(0,0), digits = 3, transf = exp),
    predict(res_mod_sens_n$model[c("rob")][[1]], 
               newscale= c(mean(data$rob,  na.rm=T),
                           mean(data$totaln_bl10, na.rm=T)), digits = 3, transf = exp),
    predict(res_mod_sens_n$model[c("rob")][[1]], 
            newscale= c(0,mean(data$totaln_bl10, na.rm=T)), digits = 3, transf = exp)),
  year100= list(
    predict(res_mod_sens_n$model[c("year100")][[1]], 
            newscale= c(0,0), digits = 3, transf = exp),
    predict(res_mod_sens_n$model[c("year100")][[1]], 
                   newscale= c(mean(data$year100, na.rm=T), 
                               mean(data$totaln_bl10, na.rm=T)), digits = 3, transf = exp),
    predict(res_mod_sens_n$model[c("year100")][[1]], 
            newscale= c(0, mean(data$totaln_bl10, na.rm=T)), digits = 3, transf = exp)),
  totaln_bl10= list(
    predict(res_mod_sens_n$model[c("totaln_bl10")][[1]], 
                       newscale= c(0), digits = 3, transf = exp),
    predict(res_mod_sens_n$model[c("totaln_bl10")][[1]], 
            newscale= c(mean(data$totaln_bl10,  na.rm=T)), digits = 3, transf = exp)),
  target_group= map(newScales(data$target_group),
                    ~predict(res_mod_sens_n$model[c("target_group")][[1]], 
                             newscale = .x, digits = 3, transf = exp)),
  instrument_red= map(newScales(data$instrument_red),
                      ~predict(res_mod$model[c("instrument_red")][[1]], 
                               newscale = .x, digits = 3, transf = exp)),
  age_group= map(newScales(data$age_group), 
                 ~predict(res_mod_sens_n$model[c("age_group")][[1]], 
                          newscale = .x, digits = 3, transf = exp)),
  comorbid_mental= map(newScales(data$comorbid_mental), 
                       ~predict(res_mod_sens_n$model[c("comorbid_mental")][[1]], 
                                newscale = .x, digits = 3, transf = exp)),
  diagnosis= map(newScales(data$diagnosis), 
                 ~predict(res_mod_sens_n$model[c("diagnosis")][[1]], 
                          newscale = .x, digits = 3, transf = exp)),
  percent_women= list(
    predict(res_mod_sens_n$model[c("percent_women")][[1]], 
            newscale= c(0,0), digits = 3, transf = exp),
    predict(res_mod_sens_n$model[c("percent_women")][[1]], 
                         newscale= c(mean(data$percent_women, na.rm=TRUE),
                                     mean(data$totaln_bl10, na.rm=T)), digits = 3, transf = exp),
    predict(res_mod_sens_n$model[c("percent_women")][[1]], 
            newscale= c(0,
                        mean(data$totaln_bl10, na.rm=T)), digits = 3, transf = exp)),
  recruitment= map(newScales(data$recruitment), 
                   ~predict(res_mod_sens_n$model[c("recruitment")][[1]],
                            newscale = .x, digits = 3, transf = exp))
  
)

res_mod_sens_n$pred_tau2$condition_arm2

## 2.5 save ----------------------------------------------------------------

# location estimates + omnibus test
tibble(
  moderator = res_mod_sens_n$model_name,
  N = map_dbl(res_mod_sens_n$model, "k"),
  loc_coef = map_dbl(res_mod_sens_n$model, "b")%>% round(3),
  loc_coef_CI_lower = map_dbl(res_mod_sens_n$model, "ci.lb")%>% round(3),
  loc_coef_CI_upper = map_dbl(res_mod_sens_n$model, "ci.ub")%>% round(3),
  pval_omnibus_test_scale_coefs = map_dbl(res_mod_sens_n$model, "QSp")%>% round(3),
  pval_likelihoodRatio = unlist(map(res_mod_sens_n$res_likelihoodRatio, ~round(.x$pval, digits = 3)))
)-> res_sensitivity_N_location

res_sensitivity_N_location



# scale estimates
tibble(
  moderator = rep(res_mod_sens_n$model_name,
                  times= map(res_mod_sens_n$model, "q") %>% unlist),
  coefficients = unlist(map(res_mod_sens_n$model, ~dimnames(.x$Z)[[2]])),
  scale_coef = unlist(map(res_mod_sens_n$model, "alpha"))%>% round(3),
  scale_coef_CI_lower= unlist(map(res_mod_sens_n$model, "ci.lb.alpha"))%>% round(3),
  scale_coef_CI_lower_profileLikelihod = map(res_mod_sens_n$profileLikelihoodCI,
                                             ~map_chr(.x, ~unlist(.x)["random2"]) %>% 
                                               .[-length(.)] %>% as.numeric()) %>% 
    unlist() %>% round(3),
  scale_coef_CI_upper= unlist(map(res_mod_sens_n$model, "ci.ub.alpha"))%>% round(3),
  scale_coef_CI_upper_profileLikelihod = map(res_mod_sens_n$profileLikelihoodCI,
                                             ~map_chr(.x, ~unlist(.x)["random3"]) %>% 
                                               .[-length(.)] %>% as.numeric()) %>% 
    unlist() %>% round(3),
  scale_coef_se = unlist(map(res_mod_sens_n$model, "se.alpha"))%>% round(3),
  tau2 = res_mod_sens_n$pred_tau2 %>% unlist %>% .[grep("pred", names(.))] %>% as.numeric() %>% na.omit %>% round(3),
  tau2_CI_lower =  res_mod_sens_n$pred_tau2 %>% unlist %>% .[grep("ci.lb", names(.))] %>% 
    as.numeric() %>% na.omit %>% round(3),
  tau2_CI_upper = res_mod_sens_n$pred_tau2 %>% unlist %>% .[grep("ci.ub", names(.))] %>% 
    as.numeric() %>% na.omit %>% round(3),
  F_t_z_val = unlist(map(res_mod_sens_n$model, "zval.alpha"))%>% round(3),
  pval = unlist(map(res_mod_sens_n$model, "pval.alpha")) %>% round(3)
) -> res_sensitivity_N_scale

res_sensitivity_N_scale

# save in excel sheet
write.xlsx(res_sensitivity_N_location, file = "results/results.xlsx",
           sheetName = "N sensitivity analysis (location +omnibus)", append = TRUE)
# add a second data set in a new worksheet
write.xlsx(res_sensitivity_N_scale, file = "results/results.xlsx", 
           sheetName="N sensitivity analysis (scale)", append=TRUE)

# save as rda file
save(res_mod_sens_n, file="results/res_mod_sens_n.rda")


