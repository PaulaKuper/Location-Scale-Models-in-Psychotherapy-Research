# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#                  Location-Scale Models in Psychotherapy Research:         #
#                       Predictors of Heterogeneity                         #
#                2. Analysis- Multi-model selection (after revision)        #
#                                                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pacman::p_load(
  dplyr,
  glmulti,
  extrafont,
  metafor,
  xlsx,
  tidyr,
  purrr,
  viridis,
  corrplot,
  vcd
)

load("results/res_mod.rda")
load("data/data_agg.rda")


# load helper functions
source("utils/multimodel.selection.helpers.R")

# load modified coef function adapted to LSM models
source("utils/coef.glmulti_mod.R")

# 1 Main analysis ------------------------------------------------------------


## 1.1 prepare dataset -----------------------------------------------------

data <- data_agg %>% 
  mutate("yi"=.g, "sei"= .g_se,
         `Year (:100)` = year/100,
         `Baseline N`= totaln_bl/10,
         `Type of IG` = condition_arm1,
         `Type of CG` = condition_arm2,
        # "Instrument"= instrument,
         `N Sessions (IG)` = n_sessions_arm1,
         `IG format` = format,
         "Country"= country,
         "Rating"= rating,
         "RoB"= rob,
         `Target group`= target_group,
         `Age group`= age_group,
         `Mental comorbidity`= comorbid_mental,
         "Diagnosis"= diagnosis,
         `Percent women`= percent_women,
         "Recruitment"= recruitment) %>% 
  select(yi, sei, `Year (:100)`,  `Baseline N`,  `Type of IG`, `Type of CG`,
          `N Sessions (IG)`, `IG format`, "Country", "Rating", "RoB",
         `Target group`, `Age group`, `Mental comorbidity`, "Diagnosis", `Percent women`,
         "Recruitment", study) %>% 
  na.omit()

unique(data$study) %>% length()


## 1.2 correlation matrix (continuous vars) ------------------------------------

windows(pointsize=15, family="Times New Roman", width= 4000, height= 4000)

data %>% 
  select(c("Year (:100)", "Baseline N", "N Sessions (IG)", "RoB", "Percent women")) %>% 
  na.omit() %>% 
  cor %>% 
  corrplot(., method= "number", col = viridis(10),
           number.cex= 1, type= "lower", tl.col = "black") # other method: color


dev.off()


## 1.3 correlation matrix (categorical vars) -------------------------------


windows(pointsize=15, family="Times New Roman", width= 4000, height= 4000)

data %>% 
  select(c("Type of IG", "Type of CG","Instrument", "IG format", "Country", "Rating",
           "Target group", "Age group", "Mental comorbidity", "Diagnosis",  
           "Recruitment")) %>% 
  na.omit() %>% 
  compute_associations(.) %>% 
  corrplot(., method= "number", col = viridis(20),
           number.cex= 1, type= "lower", tl.col = "black",
           col.lim=c(0,1)) # other method: color

dev.off()

## 1.4 combine moderators --------------------------------------------------

moderators <- c("`Year (:100)`", "`Baseline N`","`Type of IG`", "`Type of CG`",
                "`N Sessions (IG)`", "`IG format`", "Country", "Rating", 
                "RoB", "`Target group`", "`Age group`", "`Mental comorbidity`", 
                "Diagnosis", "`Percent women`", "Recruitment"
               )


# model all possible combinations between moderators (without interaction) 
sizes <- seq(1, 4) # allow combinations from 1 - 4 moderators
combinations <- unlist(lapply(sizes, function(size) combn(moderators, size, 
                                                          simplify = FALSE)), 
                       recursive = FALSE)

# add intercept only model
combinations <- c("1", combinations) 
combinations

# save as formulas in a list
formulas <- lapply(combinations, 
                   function(comb) as.formula(paste("~", paste(comb, collapse = "+"))))

formulas



## 1.5 fit all possible models ---------------------------------------------
res_mod_sel <- lapply(formulas, function(form) fit_all_models(form, data)) 


# save
#save(res_mod_sel, file= "results/res_mod_sel_complete.rda")
load("results/res_mod_sel_complete.rda")

## 1.6  check for errors -------------------------------------------------


# apply the getfit function to each model in coffee and store the results
model_coefs <- lapply(res_mod_sel, getfit)  

# separate location and scale intercept
model_coefs<- lapply(model_coefs, function(x) {
  dimnames(x)[[1]][1] <- c("intercept_l")
  return(x)
})

# check for missing values and identify problematic models
models_fit_issue <- map_lgl(model_coefs, ~ any(check_missing(.))) %>% which()


### 1.6.1 inspect models ----------------------------------------------------

res_mod_sel[models_fit_issue]

# check for multicollinearity
map(res_mod_sel[models_fit_issue], ~vif(., table= TRUE))



## 1.7 rank models  ---------------------------------------------------------
# determine the 'best' model by ranking the models based on a model selection 
# criterion (i.e. Akaike’s Information Criterion) using the R package ‘glmulti’


res_glmulti <- glmulti(res_mod_sel, level=1,  crit="aicc", 
                       includeobjects = TRUE)


# correct model formulas (needed for multimodel inference)
formulas_to_AICc=list(
  "AICc"=unlist(lapply(res_mod_sel, function(x) summary(x)$fit.stats["AICc","ML"])),
  "formulas"= eval(formulas)
)

# Get the indices that sort the values in ascending order
sorted_indices <- order(formulas_to_AICc$AICc)

# Reorder the formulas based on the sorted indices
res_glmulti@formulas <- formulas_to_AICc$formulas[sorted_indices]

# inspect
print(res_glmulti)

# get top models
top <- weightable(res_glmulti)
top <- top[top$aicc <= min(top$aicc) + 2,]
top

# get top 10 models
top10 <- weightable(res_glmulti)[1:10,] 
top10$weights <- top10$weights %>% round(3)
top10  

# examine best model
rma(yi=.g, sei=.g_se, scale= ~`Baseline N` + Country + RoB + Recruitment, 
    data=data, method= "REML",  test="knha", control= list(optimizer="BFGS")) 


write.xlsx(top10, file = "results/results.xlsx", 
           sheetName="Multimodel Selection (Top 10), Complete data", append=TRUE)
save(res_glmulti, file= "results/res_glmulti_complete.rda")



## 1.8 plots  -----------------------------------------------
# Export the plot as JPG with a resolution of at least 300 dpi

# importance plot: importance = sum of the weights for the models in which the variable appears
# 0.8 = cutoff to differentiate between (un-)important variables
jpeg("results/plots/varImp_multimodelSelection_complete.jpg", width = 2000, height = 2000, 
     quality = 100, units = "px", res = 300)
par(mar = c(5, 5, 4, 2) + 0.6)
plot(res_glmulti, type= "s", family = "Times_New_Roman")
dev.off()

# IC profile
jpeg("results/plots/ICprofile_multimodelSelection_complete.jpg", width = 2000, height = 2000, 
     quality = 100, units = "px", res = 300)
plot(res_glmulti, type= "p", family = "Times_New_Roman")
dev.off()


## 1.9 model-averaged coefs ---------------------------------------------
# The relative variable importance is calculated by adding the weights of 
# each model in which the variable occurs. Variables that appear in a lot of
# models with large weights have thus higher importance values.

  
# model-averaging (done through: coef)
coef(res_glmulti, varweighting="Johnson", select= "all", 
     models=res_mod_sel, formulas = formulas)

mmi <- as.data.frame(coef(res_glmulti, varweighting="Johnson", select= "all", 
                          models=res_mod_sel, formulas = formulas))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]

# inspect model-averaged parameter estimates
mmi <- round(mmi, 3)
mmi

write.xlsx(mmi, file = "results/results.xlsx", 
           sheetName="Multimodel Selection  full data", append=TRUE)


# 4 sensitivity analysis: outliers ------------------------------------------

## 4.1 prepare dataset -----------------------------------------------------

data <- data_agg %>% 
  mutate("yi"=.g, "sei"= .g_se,
         `Year (:100)` = year/100,
         `Baseline N`= totaln_bl/10,
         `Type of IG` = condition_arm1,
         `Type of CG` = condition_arm2,
         # "Instrument"= instrument,
         `N Sessions (IG)` = n_sessions_arm1,
         `IG format` = format,
         "Country"= country,
         "Rating"= rating,
         "RoB"= rob,
         `Target group`= target_group,
         `Age group`= age_group,
         `Mental comorbidity`= comorbid_mental,
         "Diagnosis"= diagnosis,
         `Percent women`= percent_women,
         "Recruitment"= recruitment) %>% 
  select(yi, sei, `Year (:100)`,  `Baseline N`,  `Type of IG`, `Type of CG`,
         `N Sessions (IG)`, `IG format`, "Country", "Rating", "RoB",
         `Target group`, `Age group`, `Mental comorbidity`, "Diagnosis", `Percent women`,
         "Recruitment", study, comparison_id) %>% na.omit()


## 4.2 combine moderators --------------------------------------------------

moderators <- c("`Year (:100)`", "`Baseline N`","`Type of IG`", "`Type of CG`",
                "`N Sessions (IG)`", "`IG format`", "Country", "Rating", 
                "RoB", "`Target group`", "`Age group`", "`Mental comorbidity`", 
                "Diagnosis", "`Percent women`", "Recruitment"
)


# model all possible combinations between moderators (without interaction) 
sizes <- seq(1, 4) # allow combinations from 1 - 4 moderators
combinations <- unlist(lapply(sizes, function(size) combn(moderators, size, 
                                                          simplify = FALSE)), 
                       recursive = FALSE)

# add intercept only model
combinations <- c("1", combinations) 
combinations

# save as formulas in a list
formulas <- lapply(combinations, 
                   function(comb) as.formula(paste("~", paste(comb, collapse = "+"))))

formulas[[1]]; formulas[[1941]]

## 4.3 fit all possible models ---------------------------------------------

# excluding outliers
data_outliers <- data %>% 
  filter(!(yi > 3 | yi < -3))

unique(data_outliers$study) %>% length()

#fit 
res_mod_sel_outliers <- lapply(formulas, function(form) fit_all_models(form, data_outliers)) 
res_mod_sel_outliers


# save
#save(res_mod_sel_outliers, file= "results/res_mod_sel_outliers_fullData.rda")
load("results/res_mod_sel_outliers_fullData.rda")

## 4.3  check for errors -------------------------------------------------


# apply the getfit function to each model in coffee and store the results
model_coefs <- lapply(res_mod_sel_outliers, getfit)  

# separate location and scale intercept
model_coefs<- lapply(model_coefs, function(x) {
  dimnames(x)[[1]][1] <- c("intercept_l")
  return(x)
})

# check for missing values and identify problematic models
models_fit_issue <- map_lgl(model_coefs, ~ any(check_missing(.))) %>% which() # no issues


## 4.5 inspect models ----------------------------------------------------

res_mod_sel_outliers[models_fit_issue] 

# check for multicollinearity
map(res_mod_sel_outliers[models_fit_issue], ~vif(., table= TRUE))


# exclude models  
res_mod_sel_outlier_reduced <- res_mod_sel_outliers[-models_fit_issue]

## 4.6 rank models  ---------------------------------------------------------
# determine the 'best' model by ranking the models based on a model selection 
# criterion (i.e. Akaike’s Information Criterion) using the R package ‘glmulti’


res_glmulti_ouliers <- glmulti(res_mod_sel_outlier_reduced, level=1,  crit="aicc", 
                       includeobjects = TRUE)

summary(res_glmulti_ouliers@objects[[1]]); summary(res_glmulti_ouliers@objects[[1940]]) # sorted as initially given
res_glmulti_ouliers@formulas[[1]]; res_glmulti_ouliers@formulas[[100]]
res_glmulti_ouliers@crits[[1]]; res_glmulti_ouliers@crits[[1940]] # critical values sorted in ascending order 

# correct model formulas (needed for multimodel inference), before all formales queal ~1:
# store AICc values for each model 
# store formulas 
formulas_to_AICc=list(
  "AICc"=unlist(lapply(res_mod_sel_outlier_reduced, function(x) summary(x)$fit.stats["AICc","ML"])), # unordered AICc
  "formulas"= eval(formulas[-models_fit_issue]) # unordered formulas
)


# Get the indices that sort the values in ascending order
sorted_indices <- order(formulas_to_AICc$AICc)


# Reorder the formulas based on the sorted indices
# -> res_glmulti_ouliers formulas object is updated to store the model formulas sorted by AICc
res_glmulti_ouliers@formulas <- formulas_to_AICc$formulas[sorted_indices]
formulas_to_AICc$formulas[sorted_indices][1:10]


# inspect
print(res_glmulti_ouliers)

# get top models
top <- weightable(res_glmulti_ouliers)
top <- top[top$aicc <= min(top$aicc) + 2,]
top

# get top 10 models
top10 <- weightable(res_glmulti_ouliers)[1:10,] 
top10$weights <- top10$weights %>% round(3)
top10  

write.xlsx(top10, file = "results/results_outliers.xlsx", 
           sheetName="Multimodel Selection (full data)", append=TRUE)


## 4.7 plots  -----------------------------------------------
# Export the plot as JPG with a resolution of at least 300 dpi

# importance plot: importance = sum of the weights for the models in which the variable appears
# 0.8 = cutoff to differentiate between (un-)important variables
jpeg("results/plots/varImp_multimodelSelection_outliers_fullData.jpg", width = 2000, height = 2000, 
     quality = 100, units = "px", res = 300)
par(mar = c(5, 5, 4, 2) + 0.6)
plot(res_glmulti_ouliers, type= "s", family = "Times_New_Roman")
dev.off()

# IC profile
jpeg("results/plots/ICprofile_multimodelSelection_outliers_fullData.jpg", width = 2000, height = 2000, 
     quality = 100, units = "px", res = 300)
plot(res_glmulti_ouliers, type= "p", family = "Times_New_Roman")
dev.off()


## 4.8 model-averaged coefs ---------------------------------------------
# The relative variable importance is calculated by adding the weights of 
# each model in which the variable occurs. Variables that appear in a lot of
# models with large weights have thus higher importance values.


# model-averaging (done through: coef)
coef(res_glmulti_ouliers, varweighting="Johnson", select= "all", 
     models=res_mod_sel_outlier_reduced, formulas = formulas[-models_fit_issue])

mmi <- as.data.frame(coef(res_glmulti_ouliers, varweighting="Johnson", select= "all", 
                          models=res_mod_sel_outlier_reduced, 
                          formulas = formulas[-models_fit_issue]))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), 
                  Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]

# inspect model-averaged parameter estimates
mmi <- round(mmi, 3)
mmi

write.xlsx(mmi, file = "results/results_outliers.xlsx", 
           sheetName="Multimodel Selection full data", append=TRUE)


# 5 sensitivity analysis: influentials  ------------------------------------------

## 5.2 prepare dataset ---------------------------------------------------------

data <- data_agg %>% 
  mutate("yi"=.g, "sei"= .g_se,
         `Year (:100)` = year/100,
         `Baseline N`= totaln_bl/10,
         `Type of IG` = condition_arm1,
         `Type of CG` = condition_arm2,
         # "Instrument"= instrument,
         `N Sessions (IG)` = n_sessions_arm1,
         `IG format` = format,
         "Country"= country,
         "Rating"= rating,
         "RoB"= rob,
         `Target group`= target_group,
         `Age group`= age_group,
         `Mental comorbidity`= comorbid_mental,
         "Diagnosis"= diagnosis,
         `Percent women`= percent_women,
         "Recruitment"= recruitment,
         comparison_id) %>% 
  select(yi, sei, `Year (:100)`,  `Baseline N`,  `Type of IG`, `Type of CG`,
         `N Sessions (IG)`, `IG format`, "Country", "Rating", "RoB",
         `Target group`, `Age group`, `Mental comorbidity`, "Diagnosis", `Percent women`,
         "Recruitment", comparison_id, study) %>% na.omit()




## 5.2 combine moderators --------------------------------------------------

moderators <- c("`Year (:100)`", "`Baseline N`","`Type of IG`", "`Type of CG`",
                "`N Sessions (IG)`", "`IG format`", "Country", "Rating", 
                "RoB", "`Target group`", "`Age group`", "`Mental comorbidity`", 
                "Diagnosis", "`Percent women`", "Recruitment"
)


# model all possible combinations between moderators (without interaction) 
sizes <- seq(1, 4) # allow combinations from 1 - 4 moderators
combinations <- unlist(lapply(sizes, function(size) combn(moderators, size, 
                                                          simplify = FALSE)), 
                       recursive = FALSE)

# add intercept only model
combinations <- c("1", combinations) 
combinations

# save as formulas in a list
formulas <- lapply(combinations, 
                   function(comb) as.formula(paste("~", paste(comb, collapse = "+"))))

formulas

## 5.3 fit all possible models ---------------------------------------------

# excluding influentials
load("results/res_sensitivity_influentials_metafor.rda")

# inspect
data_agg %>% 
  filter(row_number() %in% which(res_influentials2$is.infl)) %>% 
  pull(comparison_id) -> exclude_influentials
exclude_influentials

# remove
data_influentials <- data %>% 
  filter(!comparison_id %in% exclude_influentials)

unique(data_influentials$study) %>% length()

# fit 
res_mod_sel_influentials <- lapply(formulas, function(form) fit_all_models(form, data_influentials)) 
res_mod_sel_influentials

# save
#save(res_mod_sel_influentials, file= "results/res_mod_sel_influentials_fullData2.rda")
load("results/res_mod_sel_influentials_fullData2.rda")

### 4.1.1  check for errors -------------------------------------------------


# apply the getfit function to each model in coffee and store the results
model_coefs <- lapply(res_mod_sel_influentials, getfit)  

# separate location and scale intercept
model_coefs<- lapply(model_coefs, function(x) {
  dimnames(x)[[1]][1] <- c("intercept_l")
  return(x)
})

# check for missing values and identify problematic models
models_fit_issue <- map_lgl(model_coefs, ~ any(check_missing(.))) %>% which() # no issues


### 4.1.2 inspect models ----------------------------------------------------

res_mod_sel_influentials[models_fit_issue] 

# check for multicollinearity
map(res_mod_sel_influentials[models_fit_issue], ~vif(., table= TRUE))


## 1.7 rank models  ---------------------------------------------------------
# determine the 'best' model by ranking the models based on a model selection 
# criterion (i.e. Akaike’s Information Criterion) using the R package ‘glmulti’


res_glmulti_influentials <- glmulti(res_mod_sel_influentials, level=1,  crit="aicc", 
                       includeobjects = TRUE)


# correct model formulas (needed for multimodel inference)
formulas_to_AICc=list(
  "AICc"=unlist(lapply(res_mod_sel_influentials, function(x) summary(x)$fit.stats["AICc","ML"])),
  "formulas"= eval(formulas)
)

# Get the indices that sort the values in ascending order
sorted_indices <- order(formulas_to_AICc$AICc)

# Reorder the formulas based on the sorted indices
res_glmulti_influentials@formulas <- formulas_to_AICc$formulas[sorted_indices]

# inspect
print(res_glmulti_influentials)

# get top models
top <- weightable(res_glmulti_influentials)
top <- top[top$aicc <= min(top$aicc) + 2,]
top

# get top 10 models
top10 <- weightable(res_glmulti_influentials)[1:10,] 
top10$weights <- top10$weights %>% round(3)
top10 



write.xlsx(top10, file = "results/results_influentials.xlsx", 
           sheetName="Multimodel Selection full data (2)", append=TRUE)
save(res_glmulti_influentials, file= "results/res_glmulti_influentials_fullData2.rda")



## 1.8 plots  -----------------------------------------------
# Export the plot as JPG with a resolution of at least 300 dpi

# importance plot: importance = sum of the weights for the models in which the variable appears
# 0.8 = cutoff to differentiate between (un-)important variables
jpeg("results/plots/varImp_multimodelSelection_influentials_fullData2.jpg", width = 2000, height = 2000, 
     quality = 100, units = "px", res = 300)
par(mar = c(5, 5, 4, 2) + 0.6)
plot(res_glmulti_influentials, type= "s", family = "Times_New_Roman")
dev.off()

# IC profile
jpeg("results/plots/ICprofile_multimodelSelection_influentials_fullData2.jpg", width = 2000, height = 2000, 
     quality = 100, units = "px", res = 300)
plot(res_glmulti_influentials, type= "p", family = "Times_New_Roman")
dev.off()


## 1.9 model-averaged coefs ---------------------------------------------
# The relative variable importance is calculated by adding the weights of 
# each model in which the variable occurs. Variables that appear in a lot of
# models with large weights have thus higher importance values.


# model-averaging (done through: coef)
coef(res_glmulti_influentials, varweighting="Johnson", select= "all", 
     models=res_mod_sel_influentials, formulas = formulas)

mmi <- as.data.frame(coef(res_glmulti_influentials, varweighting="Johnson", select= "all", 
                          models=res_mod_sel_influentials, formulas = formulas))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]

# inspect model-averaged parameter estimates
mmi <- round(mmi, 3)
mmi

write.xlsx(mmi, file = "results/results_influentials.xlsx", 
           sheetName="Model-averaged est. full data2", append=TRUE)
