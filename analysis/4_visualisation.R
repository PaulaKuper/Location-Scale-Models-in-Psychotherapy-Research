# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#                  Location-Scale Models in Psychotherapy Research:         #
#                       Predictors of Heterogeneity                         #
#                       4. Plots and visualisation                          #
#                                                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# load packages
pacman::p_load(
  metafor,
  dplyr,
  purrr,
  xlsx,
  robvis,
  extrafont,
  viridis,
  readxl,
  ggplot2,
  grid,
  gridExtra
)

load("results/res_mod.rda")
load("data/data_agg.rda")

data <- data_agg

source("utils/generate.plots.R")

# 1 RoB plots -------------------------------------------------------------

# rob = overall risk of bias score. ranging from 0 (high risk) to 5 (low risk)
# ac = allocation concealment (0= high risk; 1= low risk)
# ba = blinding of assessors (0= high risk; 1= low risk; sr= self-report)
# itt = intention-to-treat analyses (0= high risk; 1= low risk)
# sg = sequence generation (0= high risk; 1= low risk)

data$rob <- as.numeric(data$rob)
data %>% 
  dplyr::select(study, rob, ac, ba, itt, sg) %>% skim()

data %>% 
  dplyr::select(study, sg, ac, ba, itt, rob) %>% 
  dplyr::mutate(
    sg = recode(sg, "1" = "Low", "0" = "High"),
    ac= recode(ac, "1" = "Low", "0" = "High"),
    ba = recode(ba, "1" = "Low", "0" = "High", "sr"= "Low", "2"= "Low"),
    itt= recode(itt, "1" = "Low", "0" = "High"),
    rob = recode(rob, "4"= "Low", "0"= "High", "1"= "High", 
             "2"= "Unclear", "3"= "Unclear"),
    Weight= rep(1, nrow(.))
    ) %>% 
  rename(
    "Study" = study,
    "Random.sequence.generation."= sg,
    "Allocation.concealment."= ac,
    "Blinding.of.outcome.assessors." = ba,
    "Intention.to.treat.analyses." = itt,
    "Overall" = rob
  ) -> rob_data 

# summary plot
jpeg("results/plots/RoB_summary_26_04_2024.jpg", width = 2000, height = 2000, 
     quality = 100, units = "px", res = 300, family = "Times New Roman")
rob_summary("ROB1", overall=TRUE, data= rob_data)
dev.off()




# 2 Bubble plots -----------------------------------------------------------
# for continuous predictors

# preparation
data$year100 = data$year/100
data$percent_women <- data$percent_women * 100

dataSessions <- data %>% dplyr::select(.g, .g_se, yi, vi, n_sessions_arm1) %>% na.omit()
dataN <- data %>% dplyr::select(.g, .g_se, yi, vi, totaln_bl) %>% na.omit()
dataWomen <- data %>% dplyr::select(.g, .g_se, yi, vi, percent_women) %>% na.omit()
dataRoB <- data %>% dplyr::select(.g, .g_se, yi, vi, rob) %>% na.omit()

# run models as basis for plots
res_sessions <- rma(yi, vi, mods = ~ 1, scale = ~ n_sessions_arm1, data=dataSessions, test="knha")
res_n <- rma(yi, vi, mods = ~ 1, scale = ~ totaln_bl, data=dataN, test="knha")
res_perc_women <- rma(yi, vi, mods = ~ 1, scale = ~ percent_women, data=dataWomen, test="knha")
res_rob <- rma(yi, vi, mods = ~ 1, scale = ~ rob, data=dataRoB, test="knha")
res_year <- rma(yi, vi, mods = ~ 1, scale = ~ year, data=data, test="knha")

# plots
jpeg(filename = "results/plots/bubble_plots.jpg",  width = 3000, height = 2000, 
     res= 300, family = "Times New Roman")
layout(matrix(c(1:6), nrow = 2, ncol = 3, byrow = TRUE)) # Set up a custom layout

generate_bubble_plot(res_year, data, xlimMin = 1970, xlimMax = 2025, ylimMax=10, moderator = "year", 
              xTicks= seq(1970, 2025, by= 5), plotLabel = "Year and Scale", legendTextSize = 0.7,
              xlab = "year", legend.position= "topleft")


generate_bubble_plot(res_sessions, dataSessions, xlimMin =1 , xlimMax = 60, ylimMax = 10, moderator= "n_sessions_arm1", 
              xTicks = seq(0, 60, by=10), plotLabel = "N Sessions and Scale", legendTextSize = 0.7,
              xlab = "number of sessions")

generate_bubble_plot(res_n, dataN, xlimMin = 0, xlimMax = 1250, ylimMax=10, moderator= "totaln_bl", 
              xTicks =seq(0, 1250, by= 100), plotLabel= "Total Sample Size and Scale", legendTextSize = 0.7,
              xlab = "total sample size (baseline)")

generate_bubble_plot(res_perc_women, dataWomen, xlimMin = 0, xlimMax = 100, ylimMax=10, moderator= "percent_women", 
                     xTicks =seq(0, 100, by= 10), plotLabel= "Percent Women and Scale", legendTextSize = 0.7,
                     legend.position= "topleft",
                     xlab = "% women")

generate_bubble_plot(res_rob, dataRoB, xlimMin = 0, xlimMax = 4, ylimMax=10, moderator= "rob", 
                     xTicks =seq(0, 14, by= 1), plotLabel= "RoB and Scale", legendTextSize = 0.7,
                     xlab = "RoB")


dev.off()


# 3 plots for categorical predictors  ---------------------------------------


# rename variable levels for plot labels
data <- data %>%
  mutate(
    condition_arm1 = recode(condition_arm1, "cbt" = "CBT", "bat" = "BAT", "pst" = "PST", "ipt" = "IPT",
                            "dyn" = "Dyn", "lrt" = "LRT", "other psy" = "Other psy", "sup"= "Sup"),  
    condition_arm2 = recode(condition_arm2, "cau" = "Care-as-usual", "wl" = "Waitlist", "other ctr" = "Other control"),  
    format = recode(format, "ind" = "Individual", "grp" = "Group", "gsh" = "Guided self-help", 
                    "tel" = "Telephone", "oth" = "Other (mixed formats)"),
    instrument_red = recode(instrument_red, "bdi-1" = "BDI-I", "bdi-2"= "BDI-II", "ces-d"="CES-D", "epds"= "EPDS",
                            "hdrs"= "HDRS", "other"= "Other", "phq-9"= "PHQ-9") %>% as.factor,
    country = recode(country, "eu" = "Europe", "us" = "USA", "oth" = "Other", "can" = "Canada", "uk" = "UK", "au" = "Australia", 
                     "eas" = "East Asia"),
    diagnosis = recode(diagnosis, "mdd" = "Major depression", "mood" = "Mood disorder", "cut" = "Cut-off score",
                       "sub" = "Subclinical depression", "chr" = "Chronic depression"),
    comorbid_mental = recode(comorbid_mental, "y" = "Yes", "n" = "No"),   
    age_group = recode(age_group, "adul" = "Adults", "old" = "Older adults", "adol" = "Adolescents", "yadul" = "Young adults",
                       "child"= "Child"),
    target_group = recode(target_group, "adul" = "Adults", "old" = "Older adults", "stud" = "Student population", 
                          "ppd" = "Women with perinatal depression", "med" = "Comorbid medical disorder", "oth" = "Other",
                          "adol"= "Adolescents", "child"= "Child"),
    recruitment = recode(recruitment, "com" = "Community", "clin" = "Clinical", "oth" = "Other"),
    rating = recode(rating, "clinician"= "Clinician", "parent_report"= "Parent-reported", "self-report"= "Self-reported")
  )

# drop studies with parent-reported outcomes as n= 2
data_rating <- subset(data, data$rating !='Parent-reported')
data_rating$rating <- data_rating$rating %>% droplevels()

    
# plots 1/2
jpeg(filename = "results/plots/densityPlots1.jpg",  width = 3000, height = 2000, res= 300, family = "Times New Roman")
par(mfrow = c(2,3))

generate_density_plots(res_mod$model[c("condition_arm1")][[1]], data, "condition_arm1", "Type of IG and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 5, y_axis_length= 5)

generate_density_plots(res_mod$model[c("condition_arm2")][[1]], data, "condition_arm2", "Type of CG and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 2, y_axis_length= 4)

generate_density_plots(res_mod$model[c("format")][[1]], data, "format", "Format and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 2, y_axis_length= 5)

generate_density_plots(res_mod$model[c("country")][[1]], data, "country", "Country and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 2, y_axis_length= 5)

generate_density_plots(res_mod$model[c("instrument_red")][[1]], data, "instrument_red", "Instrument and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 2, y_axis_length= 5)

generate_density_plots(res_mod$model[c("rating")][[1]], data_rating, "rating", "Rating and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 2, y_axis_length= 4)

par(mfrow = c(1, 1))
dev.off()

# plots 2/2
jpeg(filename = "results/plots/densityPlots2.jpg",  width = 3000, height = 2000, res= 300, family = "Times New Roman")
par(mfrow = c(2,3))

generate_density_plots(res_mod$model[c("target_group")][[1]], data, "target_group", "Target group and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 2, y_axis_length= 5)

generate_density_plots(res_mod$model[c("age_group")][[1]], data, "age_group", "Age group and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 2, y_axis_length= 4)

generate_density_plots(res_mod$model[c("comorbid_mental")][[1]], data, "comorbid_mental", "Mental comorbidity and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 2, y_axis_length= 4)

generate_density_plots(res_mod$model[c("diagnosis")][[1]], data, "diagnosis", "Diagnosis and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 2, y_axis_length= 4)

generate_density_plots(res_mod$model[c("recruitment")][[1]], data, "recruitment", "Recruitment and Scale",  
                       legendTextSize = 0.7,  x_axis_length = 2, y_axis_length= 5)

par(mfrow = c(1, 1))
dev.off()





# 4 Plot for main manuscript --------------------------------------------
load("data/data_agg.rda")

## 4.1 Bubbles ---------------------------------------------------------------

# fit models to create plots
dataTotalN_bl<- data %>% dplyr::select(.g, .g_se, yi, vi, totaln_bl) %>% na.omit()

res_year <- rma(yi, vi, mods = ~ 1, scale = ~ year, data=data, test="knha")
res_totalnbl <- rma(yi, vi, mods = ~ 1, scale = ~ totaln_bl, data=dataTotalN_bl, test="knha")



# plots
windows(pointsize=24, family="Times New Roman", width= 6000, height= 3000)
par(mfrow = c(1,2))

generate_bubble_plot(res_totalnbl, dataTotalN_bl, xlimMin = 0, xlimMax = 1250, ylimMax=8, moderator= "totaln_bl", 
                     xTicks =seq(0, 1250, by= 50), plotLabel= "a) Sample Size and Scale", legendTextSize = 0.4,
                     xlab = "Sample size", legend.box.size = 0.6, color= "#287D8EFF")

generate_bubble_plot(res_year, data, xlimMin = 1970, xlimMax = 2025, ylimMax=8, moderator = "year", 
                     xTicks= seq(1970, 2025, by= 5), plotLabel = "b) Year and Scale", legendTextSize = 0.4,
                     xlab = "Year", legend.box.size = 0.6, legend.position = "topleft", color= "#39558CFF")


dev.off()


## 4.2 Raincloud plot ------------------------------------------------------
load("data/data_agg.rda")

# Load fonts from the extrafont database
loadfonts()

data <- within(data_agg, {
  country = fct_relevel(country, "eu", "au", "can", "eas", "uk", "us", "oth")
  rob= as.factor(rob)
})

# recode factor levels for visualisation
data <- data %>% 
  mutate(country = fct_recode(country,
                              "Europe" = "eu",
                              "Australia" = "au",
                              "Canada" = "can",
                              "East Asia" = "eas",
                              "UK" = "uk",
                              "US" = "us",
                              "Other" = "oth")) %>% 
  mutate(rob = case_when(
    rob == 4 ~ "Low",
    rob %in% c(2, 3) ~ "Moderate",
    rob %in% c(0, 1) ~ "High"
  ))

# Convert 'rob' back to factor with specified levels and order
data$rob <- as.factor(data$rob)
data$rob <- fct_relevel(data$rob, "High", "Moderate", "Low")
data$country <- fct_relevel(data$country,  "Other" ,"US" ,"UK"  , "East Asia" ,"Canada","Australia" ,  "Europe"  )

# plot data for RoB
data_rob<-data.frame(
  moderator_level= data$rob,
  tau2i= pmax(0, resid(res_mod$model[c("rob")][[1]])^2/1-
                hatvalues(res_mod$model[c("rob")][[1]])-data$vi) %>% round(2)
)

# plot data for region
data_country <-data.frame(
  moderator_level= data$country,
  tau2i= pmax(0, resid(res_mod$model[c("country")][[1]])^2/1-hatvalues(res_mod$model[c("country")][[1]])-data$vi)
)



# create the raincloud plot for region
p_country <- ggplot(data_country, aes(x = tau2i, y = moderator_level, fill = moderator_level, color = moderator_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.5), size = 3, alpha = 0.5, shape = 21) +
  stat_halfeye(adjust = 19, width = 0.4, justification = -0.15, .width = 0.5, point_color = NA, scale=0.7) +  # Adjust width and adjust parameters
  scale_fill_viridis_d(option = "D", begin = 0.1, end = 1, direction = 1) +  # Viridis color palette
  scale_color_viridis_d(option = "D", begin = 0.1, end = 1, direction = 1) +
  labs(
    y = "",
    x = "Estimate of \u03C4\u00B2"
  ) +
  theme_minimal()+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size= 24),
    axis.title.x = element_text(size = 20, family = "Times New Roman"),  # Adjust x-axis title size
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 20), 
    axis.text.x = element_text(size = 20, family = "Times New Roman"),  # Adjust x-axis label size
    axis.text.y = element_text(size = 20, family = "Times New Roman")  # Adjust y-axis label size# Adjust y-axis title margin
  ) +
  scale_x_continuous(limits = c(0, 10))    +
  ggtitle("c) Region and Scale")             


# create the raincloud plot for RoB
p_rob <- ggplot(data_rob, aes(x = tau2i, y = moderator_level, fill = moderator_level, color = moderator_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.02), 
             size = 3, alpha = 0.5, shape = 21) +
  stat_halfeye(adjust = 19, width = 0.4, justification = -0.15, .width = 0.5, point_color = NA,
               position= "dodge", scale=0.7) +
  scale_fill_viridis_d(option = "D", begin = 0, end = 1, direction = 1) +  # Viridis color palette
  scale_color_viridis_d(option = "D", begin = 0, end = 1, direction = 1) +
  labs(
    y = "",
    x = "Estimate of \u03C4\u00B2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, family = "Times New Roman", face = "bold", size = 24),  # Adjust plot title properties
    axis.title.x = element_text(size = 20, family = "Times New Roman"),  # Adjust x-axis title size
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 20),  # Adjust y-axis title margin and size
    axis.text.x = element_text(size = 20, family = "Times New Roman"),  # Adjust x-axis label size
    axis.text.y = element_text(size = 20, family = "Times New Roman")  # Adjust y-axis label size
  ) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5)))+
  ggtitle("d) RoB and Scale")



windows(family= "Times New Roman", width = 6000, height= 3000)
grid.arrange(p_country, p_rob, ncol= 2, widths= c(20,20))


