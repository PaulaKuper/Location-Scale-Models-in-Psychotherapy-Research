# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#                  Location-Scale Models in Psychotherapy Research:         #
#                       Predictors of Heterogeneity                         #
#                   2. Analysis- Multi-model selection                      #
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

# 1 preparation ------------------------------------------------------------


## 1.1 prepare dataset -----------------------------------------------------

data <- data_agg %>% 
  mutate("sei"= .g_se,
         `Year (:100)` = year/100,
         `Baseline N`= totaln_bl/10,
         `Type of IG` = condition_arm1,
         `Type of CG` = condition_arm2,
         "Instrument"= instrument,
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
         "Recruitment"= recruitment)







## 1.2 correlation matrix (continuous vars) ------------------------------------

windows(pointsize=15, family="Times New Roman", width= 4000, height= 4000)

data %>% 
  select(c("year", "Baseline N", "N Sessions (IG)", "RoB", "Percent women")) %>% 
  mutate(
    "Year"=year)  %>% select(-year) %>% 
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

moderators <- c("`Year (:100)`", "`Baseline N`", "`Type of IG`", "`Type of CG`",
                "`N Sessions (IG)`", "`IG format`", "Country", "Rating", 
                "RoB", "`Target group`", "`Age group`", "`Mental comorbidity`", 
                "Diagnosis", "`Percent women`", "Recruitment")


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


# 2 multimodel selection --------------------------------------------------


## 2.1 fit all possible models ---------------------------------------------
res_mod_sel <- lapply(formulas, function(form) fit_all_models(form)) 

save(res_mod_sel, file= "results/res_mod_sel.rda")


### 2.1.1  check for errors -------------------------------------------------


# apply the getfit function to each model in coffee and store the results
model_coefs <- lapply(res_mod_sel, getfit)  

# separate location and scale intercept
model_coefs<- lapply(model_coefs, function(x) {
  dimnames(x)[[1]][1] <- c("intercept_l")
  return(x)
})

# check for missing values and identify problematic models
models_fit_issue <- map_lgl(model_coefs, ~ any(check_missing(.))) %>% which()


### 1.1.2 inspect models ----------------------------------------------------

res_mod_sel[models_fit_issue]

# check for multicollinearity
map(res_mod_sel[models_fit_issue], ~vif(., table= TRUE))


### 1.1.3 exclude models  ---------------------------------------------------

res_mod_sel_reduced <- res_mod_sel[-models_fit_issue]


## 2.2 rank models  ---------------------------------------------------------
# determine the 'best' model by ranking the models based on a model selection 
# criterion (i.e. Akaike’s Information Criterion) using the R package ‘glmulti’


res_glmulti <- glmulti(res_mod_sel_reduced, level=1,  crit="aicc", 
                       includeobjects = TRUE)


# correct model formulas (needed for multimodel inference)
formulas_to_AICc=list(
  "AICc"=unlist(lapply(res_mod_sel_reduced, function(x) summary(x)$fit.stats["AICc","ML"])),
  "formulas"= eval(formulas[-models_fit_issue])
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
rma(yi=.g, sei=.g_se, scale= ~`Baseline N` + Country + RoB + `Percent women`, 
    data=data, method= "REML",  test="knha", control= list(optimizer="BFGS")) 


write.xlsx(top10, file = "results/results.xlsx", 
           sheetName="Multimodel Selection (Top 10)", append=TRUE)
save(res_glmulti, file= "results/res_glmulti.rda")



## 2.3 plots  -----------------------------------------------
# Export the plot as JPG with a resolution of at least 300 dpi

# importance plot: importance = sum of the weights for the models in which the variable appears
# 0.8 = cutoff to differentiate between (un-)important variables
jpeg("results/plots/varImp_multimodelSelection.jpg", width = 2000, height = 2000, 
     quality = 100, units = "px", res = 300)
par(mar = c(5, 5, 4, 2) + 0.6)
plot(res_glmulti, type= "s", family = "Times_New_Roman")
dev.off()

# IC profile
jpeg("results/plots/ICprofile_multimodelSelection.jpg", width = 2000, height = 2000, 
     quality = 100, units = "px", res = 300)
plot(res_glmulti, type= "p", family = "Times_New_Roman")
dev.off()


# 3 model-averaged coefs ---------------------------------------------
# The relative variable importance is calculated by adding the weights of 
# each model in which the variable occurs. Variables that appear in a lot of
# models with large weights have thus higher importance values.

  
# model-averaging (done through: coef)
coef(res_glmulti, varweighting="Johnson", select= "all", 
     models=res_mod_sel_reduced, formulas = formulas[-models_fit_issue])

mmi <- as.data.frame(coef(res_glmulti, varweighting="Johnson", select= "all", 
                          models=res_mod_sel_reduced, formulas = formulas[-models_fit_issue]))
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
           sheetName="Multimodel Selection (Model-averaged Coefs)", append=TRUE)

