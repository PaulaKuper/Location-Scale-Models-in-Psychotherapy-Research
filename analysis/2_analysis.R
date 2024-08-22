# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#                  Location-Scale Models in Psychotherapy Research:         #
#                       Predictors of Heterogeneity                         #
#                          2. Analysis - univariate                         #
#                                                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# the code is based on the analysis code of Viechtbauer & Lopéz-Lopéz (2022):
# https://osf.io/qkxcu

# Note: The terms "predictor" and "scale moderator" are used interchangeably 

# load packages
pacman::p_load(
  metafor,
  dplyr,
  glmulti,
  purrr,
  tidyr,
  xlsx,
  skimr,
  forcats
)

# load prepared data with aggregation of effect sizes on an arm-level
load("data/data_agg.rda")
data <- data_agg

# load study names for that the format got aggregated
load("results/comparisons_aggregate_format.rda")

# source function to fit univariate moderations 
source("utils/fit.models.R")

# 0 preparation -----------------------------------------------------------
# reorder factor levels
# R automatically dummy codes the factor variable, using the first level as the 
# reference level, and creates dummy variables for the remaining levels;
# results of the omnibus test of a factor will be identical regardless of which 
# level is chosen as the reference level 
# see: https://www.metafor-project.org/doku.php/tips:models_with_or_without_intercept)

data <- within(data, {
  condition_arm1 = fct_relevel(condition_arm1, "cbt", "3rd","bat", "dyn", "ipt", "lrt", "other psy", "pst", "sup")
  format= fct_relevel(format, "ind", "grp", "gsh", "tel", "oth")
  target_group = fct_relevel(target_group, "adul", "child", "adol", "stud", "old", "med", "ppd", "oth")
  country= fct_relevel(country, "eu", "au", "can", "eas", "uk", "us", "oth")
  age_group = fct_relevel(age_group, "adul", "child", "adol", "yadul", "old")
  diagnosis = fct_relevel(diagnosis, "mdd", "sub", "mood", "cut", "chr")
  instrument_red= as.factor(instrument_red)
})


# 1 main effect ---------------------------------------------------------

res <- rma(yi=.g, sei=.g_se, data= data, method= "REML", test="knha")
summary(res)
confint(res, type= "PL", digits= 2)

# obtain prediction interval
predict(res,  digits= 3) 



# 2 univariate analyses -------------------------------------------------

# use contrast coding for factors to test whether there are differences
# between the different factor levels; results of the omnibus test of a factor 
# will be identical regardless of which level is chosen as the reference level
# see also: https://www.metafor-project.org/doku.php/tips:models_with_or_without_intercept

# transform variables with very different scales compared to the standardized 
# effect size measure to a comparable scale
data$year100 <- data$year/100
data$totaln_bl10 <- data$totaln_bl/10 


# create a reduced dataset for the moderation analysis of "format" without 
# studies with aggregated values
data_format <- subset(data, !(comparison_id %in% comparisons_aggregate_format))

# drop studies with parent-reported outcomes as n= 2
data_rating <- subset(data, data$rating !='parent_report')
data_rating$rating <- data_rating$rating %>% droplevels()

# run univariate moderation models
res_mod <- tibble(
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
    ~condition_arm1+rob,
    ~condition_arm2+rob,          # control group condition
    ~n_sessions_arm1+rob,         # IG: Average number of sessions received
    ~format +rob,                 # therapy format (ind= individual; grp= group; gsh= guided self-help; tel= telephone; cpl= couple therapy; oth= other (mixed formats); ush= unguided self-help
    
    
    # Study characteristics
    ~country+rob,                 # us= USA; uk= United Kingdom; eu= Europe; can= Canada; au= australia; eas= east asia; oth= other
    ~rating+rob,                  # self-reported (\"self-report\") or clinician-rated (\"clinician\")
    ~rob,                         # overall risk of bias: 0 (high risk)-4 (low risk))
    ~year100+rob,                 # year of publication
    ~totaln_bl10+rob,             # total N across both groups at baseline
    ~target_group+rob,            # adul= adults, old= older adults, stud= student population, ppd= women with perinatal depression; med= comorbid medical disorder; oth= other
    ~instrument_red+ rob,
    
    # Participant characteristics (study-level)
    ~age_group+rob,               # adul= adults, 5= older adults (≥55 years); old= older old adults (≥75 years)
    ~comorbid_mental+rob,         # comrbid mental disorder at baseline (yes/no))
    ~diagnosis+rob,               # mdd= major depression; mood= mood disorder; cut= cut-off score; sub= subclinical depression;chr= chronic depression
    ~percent_women+rob,           # % of women at baseline
    ~recruitment+rob              # com= community; clin= clinical; oth= other
  ),
  
  database = list(data, data, data, data_format, data, data_rating, data, data, 
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
res_mod$model

# inspect missingness for each predictor
res_mod$n_missing  
map(res_mod$n_missing, ~./nrow(data))# percent missingness

# name models
names(res_mod$model) <- res_mod$model_name

# save
save(res_mod, file= "results/res_mod.rda")


# 3 likelihood ratio tests ----------------------------------------------
# conduct a likelihood ratio test comparing a model that allows τ2 to differ 
# across subgroups with a model containing only RoB as scale predictor

reduced_model <- rma(yi=.g, sei=.g_se, scale= ~rob, data=data, method= "REML",
                     test="knha", control= list(optimizer="BFGS"))

reduced_model_format <- rma(yi=.g, sei=.g_se, scale= ~rob, data=data_format, 
                            method= "REML", test="knha", control= list(optimizer="BFGS"))

reduced_model_rating <- rma(yi=.g, sei=.g_se, scale= ~rob, data=data_rating, 
                            method= "REML", test="knha", control= list(optimizer="BFGS"))

reduced_model_session <- data %>% 
  select(.g, .g_se, n_sessions_arm1, rob) %>% drop_na() %>% 
  rma(yi=.g, sei=.g_se, scale= ~rob, data=., method= "REML", test="knha", 
      control= list(optimizer="BFGS"))

reduced_model_baselineN <- data %>% 
  select(.g, .g_se, totaln_bl10, rob) %>% drop_na() %>% 
  rma(yi=.g, sei=.g_se, scale= ~rob, data=., method= "REML", test="knha", 
      control= list(optimizer="BFGS"))

reduced_model_percentWomen <- data %>% 
  select(.g, .g_se, percent_women, rob) %>% drop_na() %>% 
  rma(yi=.g, sei=.g_se, scale= ~rob, data=., method= "REML", test="knha", 
      control= list(optimizer="BFGS"))

reduced_model_rob <- rma(yi=.g, sei=.g_se, scale= ~1, data=data, method= "REML", 
                         test="knha", control= list(optimizer="BFGS"))

# likelihood ratio tests (separately for format, rating and predictors with missing values)
res_likelihoodRatio <- list(
 n_sessions_arm1 = anova(res_mod$model[c("n_sessions_arm1")][[1]], reduced_model_session),
 format = anova(res_mod$model[c("format")][[1]], reduced_model_format),
 rating = anova(res_mod$model[c("rating")][[1]], reduced_model_rating),
 totaln_bl10 = anova(res_mod$model[c("totaln_bl10")][[1]], reduced_model_baselineN),
 percent_women = anova(res_mod$model[c("percent_women")][[1]], reduced_model_percentWomen),
 rob = anova(res_mod$model[c("rob")][[1]], reduced_model_rob),
 map(res_mod$model[c( "condition_arm1", "condition_arm2", "country",  
                      "year100", "target_group", "instrument_red",  "age_group", 
                     "comorbid_mental", "diagnosis",  "recruitment")], 
    ~anova(.x,reduced_model))
) 

# order the results as in res_mod and save in results list
res_mod$res_likelihoodRatio <- c(res_likelihoodRatio[-length(res_likelihoodRatio)], 
                                 res_likelihoodRatio[[7]])%>% 
  .[match(res_mod$model_name, names(.))]


res_mod$res_likelihoodRatio



# 4 profile likelihood CIs ------------------------------------------------
# Profile likelihood CIs for the scale coefficients as the sampling distribution 
# of the scale coefficients may not be not normal here. For the calculation, 
# the scale parameter is varied over a range of values while keeping the other 
# parameters fixed.  The interval is then constructed by finding the range of 
# values for which the (restricted) log-likelihood falls within a threshold of 
# the REML value that is determined by the confidence level and the degrees of 
# freedom (Viechtbauer & López-López, 2022).


# for all non-fixed scale parameters
res_mod$profileLikelihoodCI <- map(res_mod$model_name, 
                                   ~confint(res_mod$model[c(.x)][[1]], digits= 3))

profileLikelihoodCI = res_mod$profileLikelihoodCI

save(profileLikelihoodCI, file= "results/res_profile_likelihood.rda")
save(res_mod, file= "results/res_mod.rda")


# 5 profile likelihood plots for model diagnostics ---------------------------
# for the scale coefficients in the model to make sure that the bounds obtained 
# for the profile likelihood CIs are sensible

load("results/res_mod.rda")

# profile likelihood plots for the scale coefficients in the model
res_profile <- map(res_mod$model_name, 
                   ~profile(res_mod$model[c(.x)][[1]]))
names(res_profile) = res_mod$model_name

save(res_profile, file= "results/res_profile.rda")


# plot 1/2
jpeg(filename = "results/plots/profile_plots1.jpg",  width = 5000, height = 2000, res= 300)
par(mfrow= c(5, 5), mar = c(2, 2, 2, 2))
# study characteristics
plot(res_profile$condition_arm2[[2]], main= "CG: Other CG")
plot(res_profile$condition_arm2[[3]], main= "CG: Waitlist")

plot(res_profile$country[[2]], main= "Country: Australia")
plot(res_profile$country[[3]], main= "Country: Canada")
plot(res_profile$country[[4]], main= "Country: East Asia")
plot(res_profile$country[[5]], main= "Country: UK")
plot(res_profile$country[[6]], main= "Country: USA")
plot(res_profile$country[[7]], main= "Country: Other")

plot(res_profile$rating[[2]], main= "Rating: Self-report")

plot(res_profile$rob[[2]], main= "RoB")
plot(res_profile$year100[[2]], main = "Year (/100)")
plot(res_profile$totaln_bl10[[2]], main= "Baseline N (/10)")

# intervention characteristics
plot(res_profile$target_group[[2]], main= "Target group: Children")
plot(res_profile$target_group[[3]], main= "Target group: Adolescents")
plot(res_profile$target_group[[4]], main= "Target group: Students")
plot(res_profile$target_group[[5]], main= "Target group: Older adults")
plot(res_profile$target_group[[6]], main= "Target group: Med")
plot(res_profile$target_group[[7]], main= "Target group: PPD")
plot(res_profile$target_group[[8]], main= "Target group: Other")

plot(res_profile$target_group[[2]], main= "Instrument: BDI-2")
plot(res_profile$target_group[[3]], main= "Instrument: CES-D")
plot(res_profile$target_group[[4]], main= "Instrument: EPDS")
plot(res_profile$target_group[[5]], main= "Instrument: HDRS")
plot(res_profile$target_group[[6]], main= "Instrument: Other")
plot(res_profile$target_group[[7]], main= "Instrument: PHQ-9")

dev.off()

# plot 2/2
jpeg(filename = "results/plots/profile_plots2.jpg",  width = 5000, height = 2000, res= 300)
par(mfrow= c(5, 5), mar = c(2, 2, 2, 2))

plot(res_profile$condition_arm1[[2]], main= "N sessions")

plot(res_profile$condition_arm1[[2]], main= "IG: 3rd")
plot(res_profile$condition_arm1[[3]], main= "IG: BAT")
plot(res_profile$condition_arm1[[4]], main= "IG: Dyn")
plot(res_profile$condition_arm1[[5]], main= "IG: IPT")
plot(res_profile$condition_arm1[[6]], main= "IG: LRT")
plot(res_profile$condition_arm1[[7]], main= "IG: Other psy")
plot(res_profile$condition_arm1[[8]], main= "IG: PST")
plot(res_profile$condition_arm1[[9]], main= "IG: Sup")

plot(res_profile$format[[2]], main= "Format: Group")
plot(res_profile$format[[3]], main= "Format: Guided self-help")
plot(res_profile$format[[4]], main= "Format: Telephone")
plot(res_profile$format[[5]], main= "Format: Other")

# aggregate participant characteristics
plot(res_profile$age_group[[2]], main= "Age group: Children")
plot(res_profile$age_group[[3]], main= "Age group: Adolescents")
plot(res_profile$age_group[[4]], main= "Age group: Young adults")
plot(res_profile$age_group[[5]], main= "Age group: Older adults")

plot(res_profile$diagnosis[[2]], main= "Diagnosis: Subclinical depression")
plot(res_profile$diagnosis[[3]], main= "Diagnosis: Mood disorder")
plot(res_profile$diagnosis[[4]], main= "Diagnosis: Cut-off score")
plot(res_profile$diagnosis[[5]], main= "Diagnosis: Chronic depression")

plot(res_profile$recruitment[[2]], main= "Recruitment: Community")
plot(res_profile$recruitment[[3]], main= "Recruitment: Other")

dev.off()


# 6 predicted tau^2 ------------------------------------------------
# computes predicted tau^2 for the specified predictors values (scale part)

# function to iterate over predictor levels for the calculation 
newScales <- function(moderator) {
  # determine number of levels
  levels_count <- length(levels(moderator))
  
  # create Diagonal Matrix
  newscales <- diag(levels_count)
  
  # set last cells to 4 (RoB =4= low risk studies)
  newscales[, levels_count] <- 4
  
  # convert the matrix to a list of vectors
  newscales_list <- split(newscales, row(newscales))
  
  # reorder the list of vectors, moving the last row to the first position 
  newscales_list_reordered <- c(newscales_list[length(newscales_list)], newscales_list[-length(newscales_list)])
  
  # add a last row to calulculate the PI for RoB
  newscales_list_reordered <-c(newscales_list_reordered, newscales_list_reordered[1])
  
  return(newscales_list_reordered)
}


# iterate over newscales_list and compute predictions (all for studies with low RoB)
res_mod$pred_tau2 <- list(
  condition_arm1= map(newScales(data$condition_arm1), 
                      ~predict(res_mod$model[c("condition_arm1")][[1]], 
                               newscale = .x, digits = 3, transf = exp)),
  condition_arm2= map(newScales(data$condition_arm2), 
                      ~predict(res_mod$model[c("condition_arm2")][[1]], 
                               newscale = .x, digits = 3, transf = exp)),
  n_sessions_arm1= list(predict(res_mod$model[c("n_sessions_arm1")][[1]], 
                                newscale= c(0,0), digits = 3, transf = exp),
                        predict(res_mod$model[c("n_sessions_arm1")][[1]], 
                                newscale= c(mean(data$n_sessions_arm1,  
                                                 na.rm=TRUE),4), digits = 3, 
                                transf = exp),
                        predict(res_mod$model[c("n_sessions_arm1")][[1]], 
                                newscale= c(0,4), digits = 3, transf = exp)), 
  format= map(newScales(data_format$format), 
              ~predict(res_mod$model[c("format")][[1]], 
                       newscale = .x, digits = 3, transf = exp)),
  country= map(newScales(data$country), 
               ~predict(res_mod$model[c("country")][[1]], 
                        newscale = .x, digits = 3, transf = exp)),
  rating= map(newScales(data_rating$rating), 
              ~predict(res_mod$model[c("rating")][[1]], 
                       newscale = .x, digits = 3, transf = exp)),
  rob=list(predict(res_mod$model[c("rob")][[1]], 
                   newscale= c(0), digits = 3, transf = exp),
           predict(res_mod$model[c("rob")][[1]], 
                   newscale= c(mean(data$rob, na.rm=TRUE)), digits = 3, 
                   transf = exp)),
  year100= list(predict(res_mod$model[c("year100")][[1]], 
                        newscale= c(0,0), digits = 3, transf = exp),
                predict(res_mod$model[c("year100")][[1]], 
                        newscale= c(mean(data$year100, na.rm=TRUE),4), digits = 3, 
                        transf = exp),
                predict(res_mod$model[c("year100")][[1]], 
                        newscale= c(0,4), digits = 3, transf = exp)),
  totaln_bl10= list(predict(res_mod$model[c("totaln_bl10")][[1]], 
                            newscale= c(0,0), digits = 3, transf = exp),
                    predict(res_mod$model[c("totaln_bl10")][[1]], 
                            newscale= c(mean(data$totaln_bl10,  na.rm=TRUE),4), 
                            digits = 3, transf = exp),
                    predict(res_mod$model[c("totaln_bl10")][[1]], 
                            newscale= c(0,4), digits = 3, transf = exp)),
  target_group= map(newScales(data$target_group), 
                    ~predict(res_mod$model[c("target_group")][[1]],
                             newscale = .x, digits = 3, transf = exp)),
  instrument_red= map(newScales(data$instrument_red),
                      ~predict(res_mod$model[c("instrument_red")][[1]], 
                               newscale = .x, digits = 3, transf = exp)),
  age_group= map(newScales(data$age_group), 
                 ~predict(res_mod$model[c("age_group")][[1]], 
                          newscale = .x, digits = 3, transf = exp)),
  comorbid_mental= map(newScales(data$comorbid_mental),
                       ~predict(res_mod$model[c("comorbid_mental")][[1]],
                                newscale = .x, digits = 3, transf = exp)),
  diagnosis= map(newScales(data$diagnosis), 
                 ~predict(res_mod$model[c("diagnosis")][[1]], 
                          newscale = .x, digits = 3, transf = exp)),
  percent_women= list(predict(res_mod$model[c("percent_women")][[1]], 
                              newscale= c(0,0), digits = 3, transf = exp),
                      predict(res_mod$model[c("percent_women")][[1]], 
                              newscale= c(mean(data$percent_women, na.rm=TRUE),4), 
                              digits = 3, transf = exp),
                      predict(res_mod$model[c("percent_women")][[1]], 
                              newscale= c(0,4), digits = 3, transf = exp)),
  recruitment= map(newScales(data$recruitment), 
                   ~predict(res_mod$model[c("recruitment")][[1]],
                            newscale = .x, digits = 3, transf = exp))
  )

# rating
res_mod$model$rating
map(newScales(data_rating$rating), ~predict(res_mod$model[c("rating")][[1]], 
                                            newscale = .x, digits = 3, transf = exp))


# N
predict(res_mod$model[c("totaln_bl10")][[1]], newscale= c(5,4), digits = 3, transf = exp) # N= 50, RoB low
predict(res_mod$model[c("totaln_bl10")][[1]], newscale= c(20,4), digits = 3, transf = exp) # N= 200, RoB low


# RoB 
predict(res_mod$model[c("rob")][[1]], newscale= c(4), digits = 3, transf = exp) # low RoB
predict(res_mod$model[c("rob")][[1]], newscale= c(0), digits = 3, transf = exp) # high RoB

# Country 
predict(res_mod$model[c("country")][[1]], newscale= c(0,0,1,0,0,0,4), digits = 3, transf = exp) # East Asia, low RoB
predict(res_mod$model[c("country")][[1]], newscale= c(0,0,0,0,0,0,4), digits = 3, transf = exp) # Europe, low RoB

# Rating 
predict(res_mod$model[c("rating")][[1]], newscale= c(1,4), digits = 3, transf = exp) # self-report, low RoB
predict(res_mod$model[c("rating")][[1]], newscale= c(0,4), digits = 3, transf = exp) # clinician, low RoB


# 7 permutation tests -----------------------------------------------------
load("results/res_mod.rda")

res_permutest <- map(res_mod$model, ~permutest(.)) 

permutest_instrument_red


# save
save(res_permutest, file= "results/res_permutest.rda")

# 8 save results --------------------------------------------------------

# location estimates + omnibus test
tibble(
  moderator = res_mod$model_name,
  N = map_dbl(res_mod$model, "k"),
  loc_coef = map_dbl(res_mod$model, "b")%>% round(3),
  loc_coef_CI_lower = map_dbl(res_mod$model, "ci.lb")%>% round(3),
  loc_coef_CI_upper = map_dbl(res_mod$model, "ci.ub")%>% round(3),
  pval_omnibus_test_scale_coefs = map_dbl(res_mod$model, "QSp")%>% round(3),
  pval_likelihoodRatio = unlist(map(res_mod$res_likelihoodRatio, ~round(.x$pval, digits = 3))),
  pvals_permutest = map_dbl(res_permutest, "QSp")%>% round(3)
)-> res_main_location

res_main_location


# scale estimates
tibble(
  moderator = rep(res_mod$model_name,
                  times= map(res_mod$model, "q") %>% unlist),
  coefficients = unlist(map(res_mod$model, ~dimnames(.x$Z)[[2]])),
  scale_coef = unlist(map(res_mod$model, "alpha"))%>% round(3),
  scale_coef_CI_lower= unlist(map(res_mod$model, "ci.lb.alpha"))%>% round(3),
  scale_coef_CI_lower_profileLikelihod = map(res_mod$profileLikelihoodCI,
                                             ~map_chr(.x, ~unlist(.x)["random2"]) %>% 
                                               .[-length(.)] %>% as.numeric()) %>% 
    unlist() %>% round(3),
  scale_coef_CI_upper= unlist(map(res_mod$model, "ci.ub.alpha"))%>% round(3),
  scale_coef_CI_upper_profileLikelihod = map(res_mod$profileLikelihoodCI,
                                             ~map_chr(.x, ~unlist(.x)["random3"]) %>% 
                                               .[-length(.)] %>% as.numeric()) %>% 
    unlist() %>% round(3),
  scale_coef_se = unlist(map(res_mod$model, "se.alpha"))%>% round(3),
  tau2 =  res_mod$pred_tau2 %>% unlist %>% .[grep("pred", names(.))] %>% as.numeric() %>% na.omit %>% round(3),
  tau2_CI_lower =  res_mod$pred_tau2 %>% unlist %>% .[grep("ci.lb", names(.))] %>% 
    as.numeric() %>% na.omit %>% round(3),
  tau2_CI_upper =  res_mod$pred_tau2 %>% unlist %>% .[grep("ci.ub", names(.))] %>% 
    as.numeric() %>% na.omit %>% round(3),
  F_t_z_val = unlist(map(res_mod$model, "zval.alpha"))%>% round(3),
  pval = unlist(map(res_mod$model, "pval.alpha")) %>% round(3),
  pvals_permutest = unlist(map(res_permutest, "pval.alpha")) %>% round(3)
) -> res_main_scale
res_main_scale

res_mod$pred_tau2$rob %>% length()


# save in excel sheet
write.xlsx(res_main_location, file = "results/results.xlsx",
           sheetName = "Main analysis (location +omnibus)", append = FALSE)

# save scale results in a new worksheet
write.xlsx(res_main_scale, file = "results/results.xlsx", 
           sheetName="Main analysis (scale)", append=TRUE)


# save as rda file
save(res_mod, file="results/res_mod.rda")


