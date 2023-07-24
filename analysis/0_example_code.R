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
  metafor  # to fit LS models
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

# grade as scale moderator, tau2 is log transformed
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
#   estimate      se     zval    pval    ci.lb    ci.ub      
# intrcpt   -3.3444  0.7451  -4.4882  <.0001  -4.8048  -1.8839  *** 
# grade2     1.2069  1.3890   0.8689  0.3849  -1.5154   3.9293      
# grade3     1.0028  1.0903   0.9197  0.3577  -1.1341   3.1397      
# grade4     0.0833  1.0674   0.0780  0.9378  -2.0087   2.1753      
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


 
# 2.2 likelihood ratio test -----------------------------------------------

# compare LS model to "unconditional" model
res_unconditional <- update(res, scale = ~ 1)
anova(res, res_unconditional)

# df     AIC     BIC    AICc   logLik    LRT   pval       QE 
# Full     5 45.6271 54.8778 47.0905 -17.8135               107.1061 
# Reduced  2 40.9886 44.6889 41.2613 -18.4943 1.3615 0.7146 107.1061 


# 2.3 profile likelihood --------------------------------------------------

# Convergence will be examined for the scale coefficients in the model using 
# profile likelihood plots, where a particular parameter is fixed and then the 
# maximized (constrained) log-likelihood is calculated over the remaining 
# parameters of the scale model (α) (Viechtbauer, 2010). On this basis, a profile of the 
# (restricted) log-likelihood is constructed.

profile(res)

# 3. sensitivity analysis -------------------------------------------------
# to examine the robustness of the results

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
#   estimate      se     zval    pval¹    ci.lb    ci.ub     
#   intrcpt   -3.3444  0.7451  -4.4882  0.0060   -4.8048  -1.8839  ** 
#   grade2     1.2069  1.3890   0.8689  0.2530   -1.5154   3.9293     
#   grade3     1.0028  1.0903   0.9197  0.2540   -1.1341   3.1397     
#   grade4     0.0833  1.0674   0.0780  0.9116   -2.0087   2.1753     
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 1) p-values based on permutation testing



# 4. multimodel inference -------------------------------------------------
# perform multimodel inference to explore multivariable LS models by 
# modeling all possible combinations between scale moderators
# for model fitting based on likelihood approaches


## 4.1 glmulti -------------------------------------------------------------
# list of possible location-scale models, fitted with a specified fitting 
# function for relevant moderators:

# select relevant scale moderators
moderators <- c("year", "grade", "length")

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
  rma(yi=yi, vi=vi, mods=~1, scale=formula, data=data)
}

# save all fitted models in a list to be passed to to the 'glmmulti' function
models <- lapply(formulas, function(form) fit_all_models(form, data=dat))


# glmulti
res <- glmulti(models, level=1,  crit="aicc", includeobjects = TRUE)


# correct model formulas (needed for multimodel inference)
formulas_to_AICc <- list(
  "AICc" = unlist(lapply(models, function(x) summary(x)$fit.stats["AICc","ML"])),
  "formulas" = formulas
)

# get the indices that sort the values in ascending order
sorted_indices <- order(formulas_to_AICc$AICc)

# reorder the formulas based on the sorted indices
res@formulas <- formulas_to_AICc$formulas[sorted_indices]


# inspect
print(res)

# glmulti.analysis
# Method: h / Fitting: glm / IC used: aicc
# Level: 1 / Marginality: FALSE
# From 8 models:
#   Best IC: 38.8542780233706
# Best model:
#   [1] "~length"
# Evidence weight: 0.496922554541912
# Worst IC: 49.1685594004916
# 2 models within 2 IC units.
# 4 models to reach 95% of evidence weight.

plot(res) # IC profile

# get top models
top <- weightable(res)
top <- top[top$aicc <= min(top$aicc) + 2,]
top

# model     aicc   weights
# 1        ~length 38.85428 0.4969226
# 2 ~year + length 40.54881 0.2129735


# plot relative importance of the model terms
# (= sum of the weights for the models in which the variable appears)
# 0.8 = cutoff to differentiate between (un-)important variables
plot(res, type="s")

# calculate importance values manually by adding the weights of each model 
# in which the relevant variable occurs
0.496922555 + 0.212973536 + 0.046702210 + 0.016088618  # length:  0.7726869
0.212973536 + 0.067221246 + 0.016088618 + 0.002861350  # year 0.2991447
0.046702210 + 0.016088618 + 0.008087513 + 0.002861350  # grade: 0.07373969




## 4.2 multimodel inference -----------------------------------------------

# load modified coef function (adapted to LS models)
source("utils/coef.glmulti_mod.R")


# model averaging (done through: coef, produces model-averaged estimates and
# unconditional CIs)
coef(res, varweighting="Johnson", select= "all", models=models, formulas_as_vector = formula_vec)
coef(res, model=models, formulas_as_vector = formula_vec) %>% round(4)

#               Estimate    Uncond. variance Nb models Importance +/- (alpha=0.05)
# grade2        0.1069           0.0617         4     0.0737           0.4868
# grade3        0.1047           0.0499         4     0.0737           0.4379
# grade4       -0.0040           0.0088         4     0.0737           0.1838
# year          0.0109           0.0006         4     0.2991           0.0473
# length        0.0412           0.0046         4     0.7727           0.1324
# intercept_l   0.2027           0.0021         8     1.0000           0.0905
# intrcpt     -25.1308        2309.8193         8     1.0000          94.1970


# reorder
mmi <- as.data.frame(coef(res, varweighting="Johnson", models=models))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]

# inspect model-averaged parameter estimates, which are weighted averages of the
# model coefficients across the various models (with weights equal to the model 
# probabilities), conditional on n^k models that we have fitted
round(mmi, 4)

#             Estimate Std. Error z value Pr(>|z|)     ci.lb   ci.ub Importance
# intercept_l   0.2027     0.0462  4.3895   0.0000    0.1122  0.2933     1.0000
# intrcpt     -25.1308    62.8479 -0.3999   0.6893 -148.3104 98.0487     1.0000
# length        0.0412     0.0690  0.5964   0.5509   -0.0941  0.1765     0.7727
# year          0.0109     0.0316  0.3445   0.7305   -0.0510  0.0727     0.2991
# grade2        0.1069     0.5604  0.1908   0.8487   -0.9915  1.2054     0.0737
# grade3        0.1047     0.4777  0.2191   0.8266   -0.8315  1.0409     0.0737
# grade4       -0.0040     0.3322 -0.0121   0.9904   -0.6552  0.6471     0.0737


# manual cross-check of estimates (estimates * model weights)
(0.0394*0.016088618 + 0.0364*0.212973536 + 0.0326*0.002861350 + 0.0357*0.067221246)    # year: 0.01087921
(0.0785*0.016088618 + 0.0768*0.046702210 + 0.0523*0.212973536 + 0.0507*0.496922555)    # length: 0.04118218
(1.4412*0.016088618 + 1.5146*0.046702210 + 1.1306*0.002861350 + 1.2069*0.008087513)    # grade2: 0.1069179
