# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#                  Location-Scale Models in Psychotherapy Research:         #
#                       Predictors of Heterogeneity                         #
#                             CODE EXAMPLE                                  #
#                                                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# the code is orientated at the "metafor" package documentation by Wolfgang Viechtbauer:
# https://www.metafor-project.org/doku.php/tips:different_tau2_across_subgroups and
# https://osf.io/qkxcu

pacman::p_load(
  dplyr,
  purrr,
  metafor,  # to fit LS models
  glmulti # for multimodel selection and inference
)

# load example dataset that is integrated in the "metafor" package
# -> the dataset contains results from n= 48 studies that examine the
# effectiveness of writing-to-learn interventions on academic achievement.
dat <- dat.bangertdrowns2004


# 1. preparation -----------------------------------------------------------

# recode factor
within(dat, {
  grade = factor(grade)
}) -> dat

table(dat$grade)

#  1  2  3  4 
# 11  6 10 21 

unique(dat$id) %>% length()  # N = 48, no multi-arm trials/ multiple outcome measures


# 2. analysis -------------------------------------------------------------

# In the presence of dependent effect sizes (e.g. due to multiple trial arms or 
# outcome measures within one trial), an additional step would be required here.


## 2.1 location-scale model  -----------------------------------------------

# grade as scale moderator, tau^2 is log transformed
res <- rma(yi, vi, scale= ~ grade,  data=dat)
res

# Location-Scale Model (k = 48; tau^2 estimator: REML)
# 
# Test for Heterogeneity:
#   Q(df = 47) = 107.1061, p-val < .0001
# 
# Model Results (Location):
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# 0.2282  0.0456  5.0034  <.0001  0.1388  0.3176  *** 
#   
#   Test of Scale Coefficients (coefficients 2:4):
#   QS(df = 3) = 1.5353, p-val = 0.6741
# 
# Model Results (Scale):
#   
#           estimate      se     zval    pval    ci.lb    ci.ub      
# intrcpt   -3.3444  0.7451  -4.4882  <.0001  -4.8048  -1.8839  *** 
# grade2     1.2069  1.3890   0.8689  0.3849  -1.5154   3.9293      
# grade3     1.0028  1.0903   0.9197  0.3577  -1.1341   3.1397      
# grade4     0.0833  1.0674   0.0780  0.9378  -2.0087   2.1753      
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# -> no evidence that studies tended to yield more heterogeneous effects across 
# different school formats (elementary, middle, high-school, college)

predict(res, newscale= c(0,0,1), transf = exp, digits= 4) # predicted tau^2 for studies in college (grade4)

# pred   ci.lb  ci.ub 
# 0.0383 0.0082 0.1794    # pred = exp(-3.3444 + 0.0833)


 
# 2.2 likelihood ratio test -----------------------------------------------

# compare LS model to "unconditional" model
res_unconditional <- update(res, scale = ~ 1)
anova(res, res_unconditional)

#          df     AIC     BIC    AICc   logLik    LRT   pval       QE 
# Full     5 45.6271 54.8778 47.0905 -17.8135               107.1061 
# Reduced  2 40.9886 44.6889 41.2613 -18.4943 1.3615 0.7146 107.1061 


# 2.3 profile likelihood --------------------------------------------------

# Convergence will be examined for the scale coefficients in the model using 
# profile likelihood plots, where a particular parameter is fixed and then the 
# maximized (constrained) log-likelihood is calculated over the remaining 
# parameters of the scale model (α) (Viechtbauer, 2010). On this basis, a 
# profile of the (restricted) log-likelihood is constructed.

profile(res)

# 3. sensitivity analysis -------------------------------------------------
# to examine the robustness of the results regarding the test for differences
# in tau^2

## 3.1 permutation test ---------------------------------------------------
permutest(res, seed=1234)

# Running 1000 iterations for an approximate permutation test of the location model.
# |==================================================| 100% elapsed=48s  
# Running 1000 iterations for an approximate permutation test of the scale model.
# |==================================================| 100% elapsed=01m 33s
# 
# Model Results (Location):
#   
#   estimate      se    zval    pval¹   ci.lb   ci.ub      
# 0.2282  0.0456  5.0034  0.0010   0.1388  0.3176  *** 
#   
#   Test of Scale Coefficients (coefficients 2:4):¹
# QS(df = 3) = 1.5353, p-val = 0.5602
# 
# Model Results (Scale):
#   
#             estimate    se     zval    pval¹    ci.lb    ci.ub     
#   intrcpt   -3.3444  0.7451  -4.4882  0.0060   -4.8048  -1.8839  ** 
#   grade2     1.2069  1.3890   0.8689  0.2530   -1.5154   3.9293     
#   grade3     1.0028  1.0903   0.9197  0.2540   -1.1341   3.1397     
#   grade4     0.0833  1.0674   0.0780  0.9116   -2.0087   2.1753     
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 1) p-values based on permutation testing



# 4. multimodel selection and inference ------------------------------------
# perform multimodel inference to explore multivariable LS models by 
# modeling all possible combinations between scale moderators (without interactions)


## 4.1 model selection -----------------------------------------------------

# select relevant scale moderators
moderators <- c("year", "grade", "length")
dat$year <- dat$year/100 # transform to comparable scale

# generate all possible combinations and save as formulas in a list
sizes <- seq(1, length(moderators)) # allow combinations from 1 to N moderators
combinations <- unlist(lapply(sizes, function(size) combn(moderators, size, 
                                                          simplify = FALSE)), 
                       recursive = FALSE) %>% 
  append("1")# add intercept only model

formulas_vec <- lapply(combinations, function(comb) paste(comb, collapse = "+")) 
formulas <- lapply(formulas_vec, function(form) formula(paste("~", form))) # as formula
formulas  # 2^k = 2^3 = 8 possible models

# define function to fit all possible models
fit_all_models <- function(formula, data) {
  rma(yi=yi, vi=vi, mods=~1, scale=formula, method="ML", data=data)
}

# save all fitted models in a list to be passed to to the 'glmmulti' function
models <- lapply(formulas, function(form) fit_all_models(form, dat))


# glmulti
res <- glmulti(models, level=1,  crit="aicc", includeobjects = TRUE)

# correct model formulas (to print results correctly, needed for multimodel inference)
formulas_to_AICc <- list(
  "AICc" = unlist(lapply(models, function(x) summary(x)$fit.stats["AICc","ML"])),
  "formulas" = formulas
)

# get the indices that sort the values in ascending order
sorted_indices <- order(formulas_to_AICc$AICc) 

# reorder the formulas based on the sorted indices, save in res object
res@formulas <- formulas_to_AICc$formulas[sorted_indices] 


# inspect
print(res)

# glmulti.analysis
# Method: h / Fitting: glm / IC used: aicc
# Level: 1 / Marginality: FALSE
# From 8 models:
#   Best IC: 38.3117180679151
# Best model:
#   [1] "~length"
# Evidence weight: 0.493386028334476
# Worst IC: 48.5885508959637
# 2 models within 2 IC units.
# 4 models to reach 95% of evidence weight..

plot(res) # IC profile

# get top models
top <- weightable(res)
top <- top[top$aicc <= min(top$aicc) + 2,]
top

#            model     aicc   weights
# 1        ~length 38.31172 0.4933860
# 2 ~year + length 39.93725 0.2188803


# plot relative importance of the model terms
plot(res, type="s")

# equivalence check: calculate importance values manually by adding the weights 
# of each model in which the relevant variable occurs
weightable(res)[grepl("length", weightable(res)$model), "weights"] %>% sum() # length:  0.7799159
weightable(res)[grepl("year", weightable(res)$model), "weights"] %>% sum() # year:  0.3059529
weightable(res)[grepl("grade", weightable(res)$model), "weights"] %>% sum() # grade: 0.07847294


## 4.2 multimodel inference -----------------------------------------------

# load modified coef function (adapted to LS models)
source("utils/coef.glmulti_mod.R")

# model averaging (done through: coef, produces model-averaged estimates and
# unconditional CIs)
coef(res, varweighting="Johnson", select= "all", models=models, 
     formulas_as_vector = formula_vec) %>% 
  round(4)

# reorder and calculate standard errors, CIs and p-values
mmi <- as.data.frame(coef(res, varweighting="Johnson", models=models))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), ]

# inspect model-averaged parameter estimates, which are weighted averages of the
# model coefficients across the various models (with weights equal to the model 
# probabilities), conditional on n^k models that are fitted
round(mmi, 4)

#             Estimate Std. Error Importance z value Pr(>|z|)     ci.lb    ci.ub
# intercept_l   0.2012     0.0454     1.0000  4.4309   0.0000    0.1122   0.2903
# intrcpt     -26.7828    64.7613     1.0000 -0.4136   0.6792 -153.7127 100.1470
# length        0.0436     0.0713     0.7799  0.6113   0.5410   -0.0962   0.1834
# year          0.0117     0.0325     0.3060  0.3585   0.7200   -0.0521   0.0754
# grade2        0.1165     0.5921     0.0785  0.1967   0.8441   -1.0441   1.2770
# grade3        0.1155     0.5061     0.0785  0.2282   0.8195   -0.8764   1.1074
# grade4       -0.0071     0.3615     0.0785 -0.0195   0.9844   -0.7156   0.7015


# equivalence check of estimates
data.frame(
  model = weightable(res)[grepl("length", weightable(res)$model), "model"], 
  weight = weightable(res)[grepl("length", weightable(res)$model), "weights"],
  estimate= c(0.0530, 0.0546, 0.0804, 0.0817) # model coefficients for 'length'
) -> model_averaged_estimates 

model_averaged_estimates
#                    model     weight estimate
# 1                ~length 0.49338603   0.0530
# 2         ~year + length 0.21888032   0.0546
# 3        ~grade + length 0.04990259   0.0804
# 4 ~year + grade + length 0.01774696   0.0817


sum(model_averaged_estimates$weight * model_averaged_estimates$estimate) 
# length: 0.04356242

## 4.3 multimodel predictions ---------------------------------------------

# predicted value for studies with 

x <- c("length"=10, # treatment length = 10 weeks
       "year" = 1996/100, # publication year 1996
       "grade2"=0, "grade3"=0, "grade4"=1 # grade = 4 = college
       )

# loop through all 8 models, compute the predicted value based on each model
preds <- list()

for (j in 1:res@nbmods) {
  
  model <- res@objects[[j]]
  vars <- names(coef((model))$alpha)[-1]
  
  if (length(vars) == 0) {
    # intercept only model (location and scale)
    preds[[j]] <- list(-3.0567,0.4674, -3.9728,-2.1406) # models[[8]] -> scale estimate, SE, and CI 
    names(preds[[j]]) <- c("pred", "se", "ci.lb", "ci.ub")
  } else {
    preds[[j]] <- predict(model, newscale=x[vars])
  }
}

# compute a weighted average (using the model weights) of the predicted values 
# across all models
data.frame(
 yhat = sum(weightable(res)$weights* sapply(preds, function(x) x$pred)) # multimodel predicted value
 ) %>% 
 mutate(
   se = sqrt(sum(weightable(res)$weights * sapply(preds, function(x) x$se^2 + (x$pred - yhat)^2))),  # compute the corresponding (unconditional) standard error 
   CI_lb= yhat + c(-1)*qnorm(.975)*se, # 95% CI 
   CI_ub= yhat + c(1)*qnorm(.975)*se) %>% 
  round(4) 

#  yhat     se   CI_lb   CI_ub
# -3.0232 0.6802 -4.3563 -1.6901  # predicted effect 

# the multimodel predicted average tau^2 for studies with treatment length = 10 weeks, 
# conducted in college and published in 1996 is 0.05 [0.01, 0.18].






