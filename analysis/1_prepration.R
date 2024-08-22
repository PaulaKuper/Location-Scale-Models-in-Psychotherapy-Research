# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#                  Location-Scale Models in Psychotherapy Research:         #
#                       Predictors of Heterogeneity                         #
#                          1. Preparation                                   #
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
  readxl,
  tidyr,
  stringr,
  purrr,
  quanteda,
  forcats,
  xlsx
)



# the big “comparisons maps
depression_database <- read_excel("data/The depression database 1 May 2023_Paula.xlsx")
depression_database_child <- read_excel("data/The depression database 1 May 2023_Paula.xlsx", 2)
The_depression_database_1_Jan_2024 <- read_excel("data/The depression database 1 Jan 2024.xlsx")
The_depression_database_1_Jan_2024_child <- read_excel("data/The depression database 1 Jan 2024.xlsx", 2)

# psy vs ctr (adults-outpatients) 
psy_adults <- read_excel("data/psy vs ctr_25_04_24.xlsx")


# psy vs ctr (children and asolescents)
psy_child_adol <- read_excel("data/ChildAdol_8_04_2024 full database.xlsx")




# 1 prepare databases -------------------------------------------------------

# check study flow from larger database to filtered adults/ child and adolescents
# datasets
psy_adults %>%  pull(study) %>% unique %>% length # N= 929 RCTs included in larger database, N= 489 included here
psy_child_adol %>%  pull(study) %>% unique %>% length # N= 84 RCTs included in larger database, N= 76 included here

# check reasons for exclusion (children), n= 8
The_depression_database_1_Jan_2024_child %>% 
  filter(!(Study %in% psy_child_adol$study)) -> excluded_child_studies 

# check reasons for exclusion (adults), n= 440
The_depression_database_1_Jan_2024 %>% 
  filter(!(Study %in% psy_adults$study)) %>% 
  select(Study, PSYvsPSY, COMvsPHA, Inpatients, Unguided, PSYvsPHA, FRMvsFRM, DISMANTLING, BUT, OTHER) %>% 
  mutate_at(vars(-Study), as.factor) %>% 
  skim()


## 1.1 adults --------------------------------------------------------------

# The dataset contains one effect size type per study (primary_calc=1, with priority for 
# 1) means and SDs  2) dichotomous outcomes 3) change, and 4) other statistics


# filter comparisons at post intervention; only one effect size per comparison
# and exclude unguided interventions
table(psy_adults$CHANGES)

psy_adults <- psy_adults %>% 
  filter(primary_calc=="1" & time == "post" & format != "ush") %>% 
  filter(CHANGES != "change to unguided"| is.na(CHANGES))


# select relevant columns
psy_adults %>% 
  select(c("study", ".g", ".g_se", "n_arm1", "n_arm2",
         "condition_arm1", "condition_arm2","format","n_sessions_arm1", 
         "multi_arm1","multi_arm2","instrument",
         
         "year","baseline_n_arm1", "baseline_n_arm2","target_group","recruitment",
         "rating","country","rob",  "ac", "sg", "itt", "ba",
         
         "age_group","percent_women",  "comorbid_mental", "diagnosis",
         "full_ref")) %>% 
  mutate(
    database = rep("adults_outpatients", nrow(.))) -> psy_adults



# compare with depression database file
exclude_studies <- depression_database %>%
  filter(str_detect(`Exclusion/change`, regex("exclude", ignore_case = TRUE)) |
           `Exclusion/change` %in% c("add in dataset  COMB vs PHA",
                                     "change to unguided", "changed to unguided (1/05/23)",
                                     "needs to be changed to unguided (01-09-23)",
                                     "update 1-05-23// NEEDS to be classified as dismantling instead of psy vs ctr")) %>% 
  filter(Study != "Rohde, 2014b")

psy_adults$study %in% exclude_studies$Study %>% table()

# exclude studies from psy_adults database
psy_adults <- psy_adults %>% 
  filter(!(study %in% exclude_studies$Study))


## 1.2  child-adol -----------------------------------------------------------

# filter relevant comparisons and select relevant columns
psy_child_adol <- psy_child_adol %>% 
  filter(time == "post" & comparison == "PSY vs CTR") %>% 
  select("study", ".g", ".g_se", "n_arm1", "n_arm2",
         "condition_arm1", "condition_arm2","format","n_sessions_ig", 
         "multi_arm1","multi_arm2","instrument",
         
         "year","baseline_n_arm1", "baseline_n_arm2","target_group","recruitment",
         "rating","country","rob",   "ac", "sg", "itt", "ba",
         
         "age_group","percent_women",  "comorbid_mental", "diagnosis", "full_ref")%>% 
  rename("n_sessions_arm1" = n_sessions_ig) %>% 
  mutate(
  database = rep("child_adol", nrow(.))
  )

# compare with depression database file
table(depression_database_child$Exclusion) # -> no further exclusion necessary



## 1.3 merge databases ------------------------------------------------------

# merge
data <- rbind(psy_adults, psy_child_adol)
skim(data)

# correct variable formats
within(data, {
  condition_arm1  = as.factor(condition_arm1)
  condition_arm2 = as.factor(condition_arm2)
  format = as.factor(format)
  n_sessions_arm1 = as.numeric(n_sessions_arm1)
  multi_arm1 = as.factor(multi_arm1)
  multi_arm2 = as.factor(multi_arm2)
  instrument = as.factor(instrument)
  recruitment = as.factor(recruitment)
  rating = as.factor(rating)
  country = as.factor(country)
  ac= as.factor(ac)
  sg = as.factor(sg)
  itt= as.factor(itt)
  ba= as.factor(ba)
  age_group = as.factor(age_group)
  target_group = as.factor(target_group)
  percent_women = as.numeric(percent_women)
  comorbid_mental = as.factor(comorbid_mental)
  diagnosis = as.factor(diagnosis)
  database = as.factor(database)
}) -> data


# harmonize variable levels
data <- data %>%
  mutate(condition_arm1 = recode(condition_arm1, "oth psy" = "other psy"),
         condition_arm2 = recode(condition_arm2, "oth" = "other ctr", "wlc" = "wl"),
         format = recode(format, "cpl" = "oth", "mixed" = "oth", tsf= "oth"),
         target_group = recode(target_group, "older adult" = "old", "4 & 5" = "old", "adol&yadul" = "oth", 
                               "young adults" = "adul", "stu" = "stud", "med+ppd" = "oth", 
                               other= "oth", "yadul"= "adul"),
         age_group = recode(age_group, "olderold" = "old", "adol&yadul" = "yadul", "adol&yadul" = "yadul", 
                            "young adults" = "yadul"),
         ba = recode(ba, "2" = "sr"),
         diagnosis = recode(diagnosis, "modd" = "mood", "MDD" = "mdd"),
         country = recode(country, "US" = "us", "other"= "oth"),
         percent_women = case_when(
           percent_women > 1 & percent_women <= 100 ~ percent_women / 100,
           percent_women > 100 ~ percent_women / 1000,
           TRUE ~ percent_women
         ))

# check for missing values in "age group" and set them to NA
data$age_group[data$age_group %in% c("Not specified", NA)] <- NA
data$age_group <- data$age_group %>% droplevels()

# create variable for total N at baseline
data$totaln_bl <- data$baseline_n_arm1 + data$baseline_n_arm2
skim(data$totaln_bl)

# check number of studies per category/ consistent category names
lapply(
  c("condition_arm1", "condition_arm2", "format", "target_group", "recruitment",
    "rating", "country", "age_group", "comorbid_mental", "diagnosis", "database",
    "ba"),
  function(var) table(data[[var]]))

# check histogram after rescaling
hist(data$percent_women)


# 2 create unique identifier ---------------------------------------------

# The number of rows of the dataset equals the number of comparisons.
# For multi-arm studies and studies using multiple measurement instruments, 
# there are multiple rows 

unique(data$study) %>% length() # N studies= 541

# create unique identifier for each study
data <- data %>% 
  unite("id", study, condition_arm1, condition_arm2, multi_arm1, multi_arm2,
        format, instrument, rating, remove= FALSE, sep="_") 
unique(data$id) %>% length() 



# 3 handle multiple outcome measures -----------------------------------------

# inspect number of instruments per study
data %>% 
  select(study,instrument) %>% 
  group_by(study) %>% 
  summarise(n_instr= n()) %>% 
  pull(n_instr) %>% 
  table()

#   1   2   3   4   6   8   9  12  16  20 
# 275 176  38  30  15   3   1   1   1   1  


# check instruments that are used as clinician rated and self-rated 
intersect(
  data %>% filter(rating== "self-report") %>% pull(instrument) %>% unique(.),    
  data %>% filter(rating== "clinician") %>% pull(instrument) %>% unique(.)       
)     
# -> "phq-9"  "qids"   "ids"  "cdrs-r"


## 3.1 harmonize instruments -----------------------------------------------
data <- data %>%
  mutate(
    instrument = case_when(
    str_detect(instrument, pattern=regex(".*sf.*12.*m.*", ignore_case = TRUE)) ~ "sf12mh",
    str_detect(instrument, pattern=regex("bdi-ii", ignore_case = TRUE)) ~ "bdi-2",
    str_detect(instrument, pattern=regex(".*bsi.*d.*", ignore_case = TRUE)) ~ "bsi-depression",
    str_detect(instrument, pattern=regex("cds", ignore_case = TRUE)) ~ "cds",
    str_detect(instrument, pattern=regex("hrsd", ignore_case = TRUE)) ~ "hdrs",
    str_detect(instrument, pattern=regex("hrds", ignore_case = TRUE)) ~ "hdrs",
    str_detect(instrument, pattern=regex("hrsd-17", ignore_case = TRUE)) ~ "hdrs-17",
    str_detect(instrument, pattern=regex("hrds-17", ignore_case = TRUE)) ~ "hdrs-17",
    str_detect(instrument, pattern=regex("hrsd-24", ignore_case = TRUE)) ~ "hdrs-24",
    str_detect(instrument, pattern=regex("k.*sads", ignore_case = TRUE)) ~ "k-sads",
    str_detect(instrument, pattern=regex("who.*qol-bref.*", ignore_case = TRUE)) ~ "who-qol-bref",
    TRUE ~ instrument
  )) %>% 
  mutate(
    instrument= recode(
      instrument,  bid= "bdi", "DASS-21-D"= "dass21-d", DSM= "dsm", 
      "hads-d-severity"= "hads-d", "qids-sr16"= "qids-sr-16", 
      "qids16sr"= "qids-sr-16", "qids-sr"= "qids-sr-16",  "dsm-iii"= "dsm-3", 
      "Acholi Psychosocial Assessment Instrument (APAI): a local instrument"= "APAI")
  )

table(data$instrument) %>% sort(., decreasing = TRUE) %>% names()


## 3.2 implement priority rule -------------------------------------------------

# implement priority rule for effect size dependency due to multiple outcome 
# measures (within studies)

# sort by frequency of use of the instrument
data <- data %>% 
filterPriorityRule(
  instrument = c(table(.$instrument) %>% sort(., decreasing = TRUE) %>%
                   names()) %>% as.factor()
  ) 

unique(data$study) %>% length() # N = 541/541 studies (no study omitted)
unique(data$id) %>% length()  # N=  674/1032 comparisons


# create new instrument variable with reduced number of levels
data <- data %>%
  group_by(instrument) %>%
  mutate(count = n()) %>%
  ungroup() %>% 
  mutate(
    instrument_red= ifelse(count>30, instrument, "other")
  ) %>% 
  select(-count)

table(data$instrument_red)


# 4 multi-arm studies ----------------------------------------------------

## 4.1 inspect -------------------------------------------------------------

# numer of multi-arm studies and arms per study
data %>% 
  select(study, multi_arm1, multi_arm2) %>% 
  na.omit() %>% group_by(study) %>% 
  summarise(n= n()+1) %>% 
  pull(n) %>% table() 
#  number of arms:               2  3  4  5  9 
#  number of multi-arm studies: 10 89 12  4  1 = 116



## 4.2 aggregation ---------------------------------------------------------
# Three-level models are not yet possible for location-scale models in ‘metafor’. 
# Therefore, effect sizes are aggregated on an arm-level.


# create ID for comparisons on an arm level for aggregation
data$comparison_id <- paste(data$study, data$condition_arm1, data$condition_arm2, sep= "")
unique(data$comparison_id) %>% length() # n= 609 unique comparisons (on study arm level)

# check if the predictor values are identical for different arms that will get aggregated
data[!is.na(data$multi_arm1),  
     c("condition_arm1", "condition_arm2",  "year", "target_group","rating",
       "country", "rob", "totaln_bl", "n_sessions_arm1", "format", "age_group", 
       "percent_women", "comorbid_mental", "diagnosis", "recruitment", "comparison_id")] %>% 
  group_by(comparison_id) %>% 
  # summarize the grouped data by applying a function to each specified column:
  # check if the length of unique values in each column is equal to 1
  # If true: means that all values in that column are the same for the grouped data
  summarize(across(all_of(
    c("year", "target_group","rating","country", "rob","totaln_bl", "n_sessions_arm1", "format",
      "age_group", "percent_women", "comorbid_mental", "diagnosis", "recruitment")
  ), ~length(unique(.x)) == 1)) -> aggregate_check

# check number of NAs per column
map(aggregate_check, ~sum(is.na(.x)))
map(aggregate_check, ~table(.x)) 
#-> totaln_bl and n_sessions_arm1 get averaged; format would get aggregated


# save studies that are affected by an aggregation of categorical values (format)
aggregate_check %>% 
  filter(format == FALSE) %>% 
  pull(comparison_id) -> comparisons_aggregate_format

save(comparisons_aggregate_format, file= "results/comparisons_aggregate_format.rda")


# convert to escalc object for aggregation
data <- escalc(measure="SMD", yi=.g, sei= .g_se, data= data)

# aggregate effect sizes 
data_agg <- aggregate(data, cluster = comparison_id, rho = 0.5)

# sanity check
data_agg %>% pull(study) %>% unique() %>% length # 538/541 studies (3 studies omitted)
data_agg %>% pull(comparison_id) %>% unique() %>% length # 609/679 comparisons

# check omitted studies
anti_join(data, data_agg, by= "comparison_id") # studies with missing effect size data

# check if predictors were aggregated
which(aggregate_check == FALSE, arr.ind = TRUE) # number participants, session and format

# save
save(data_agg, file= "data/data_agg.rda")


# save list with full references in excel sheet
write.xlsx(data_agg %>% select(study, full_ref) %>% unique(), 
           file = "results/included_studies.xlsx",
           sheetName = "References", append = FALSE)


# 5 study characteristics ------------------------------------------------

load("data/data_agg.rda")

## 5.1 N -------------------------------------------------------------------
unique(data_agg$study) %>% length()

# n participants across comparisons
sum(data_agg$totaln_bl, na.rm = T)

# studies and comparisons per database
data_agg %>% 
  filter(database=="adults_outpatients") %>% 
  pull(study) %>% 
  unique()

data_agg %>% 
  filter(database=="child_adol") %>% 
  pull(study) %>% 
  unique()



## 5.2 other characteristics -----------------------------------------------

# intervention and control
table(data_agg$condition_arm1)/nrow(data_agg)
table(data_agg$condition_arm2)/nrow(data_agg)
table(data_agg$format)/nrow(data_agg)
skim(data_agg$n_sessions_arm1)

# diagnosis
table(data_agg$diagnosis)/nrow(data_agg)
table(data_agg$rating)/nrow(data_agg)
table(data_agg$comorbid_mental)/nrow(data_agg)


# RoB
table(data_agg$rob)/nrow(data_agg)
table(data_agg$sg)/nrow(data_agg)
table(data_agg$ac)/nrow(data_agg)
table(data_agg$ba)/nrow(data_agg)
table(data_agg$itt)/nrow(data_agg)


## 5.3 create HTML table -------------------------------------------------
# prepare data structure
data_agg$instrument <- as.factor(data_agg$instrument)
data_agg$percent_women100 <- data_agg$percent_women*100
skim(data_agg)

createStudyTable(
  data_agg,
  
  ### 5.3.1 select columns --------------------------------------
  study, .g, .g_se, 
  condition_arm1= c("CBT"= "cbt", "BAT"= "bat", "PST"= "pst", "IPT"= "ipt",
                    "Dyn"= "dyn", "LRT"= "lrt", "Other psy"= "other psy"),	
  condition_arm2 = c("CAU"= "cau", "Wl"= "wl", "Other ctr"= "other ctr"),	
  baseline_n_arm1, baseline_n_arm2,
  format = c("individual"= "ind", "group"= "grp", "guided self-help"= "gsh", 
             "telephone"="tel", "other (mixed formats)"="oth"),
  n_sessions_arm1, 
  instrument = c("APAI"="Acholi Psychosocial Assessment Instrument (APAI): a local instrument"),
  rating,
  country= c("Europe"= "eu", "USA"="us", "Other"= "oth", "Canada"= "can", "UK"= "uk", "Australia"= "au", "East Asia"= "eas"),
  diagnosis= c("major depression"="mdd", "mood disorder"="mood", "cut-off score"="cut",
                "subclinical depression"="sub", "chronic depression"="chr"),
  comorbid_mental = c("yes"= "y", "no"= "n"),	 
  age_group= c("adults"="adul", "older adults"="old", "adolescents"="adol", "young adults"="yadul"),
  target_group = c("adults"="adul", "older adults"="old", "student population"="stud", 
                   "women with perinatal depression"="ppd", "comorbid medical disorder"="med","other"="oth"),
  percent_women100, 
  recruitment = c("community"="com", "clinical"="clin", "other"="oth"),

  ac = c("high risk"="0", "low risk"= "1"),
  sg = c("high risk"="0", "low risk"= "1"),
  itt = c("high risk"="0", "low risk"= "1"),
  ba = c("high risk"="0", "low risk"= "1", "self-report"= "sr"),
  rob,
  
  ### 5.3.2 specifications -------------------------------------
  # .round.by.digits controls the number of rounded digits 
  .round.by.digits = list(
    .g= 2,
    .g_se=2,
    n_sessions_arm1=0,
    baseline_n_arm1 = 0,
    baseline_n_arm2 = 0,
    percent_women100=0
    ),
  
  # .column.names allows to rename columns
  .column.names = list(
    study= "Study",
    .g = "Hedge's g",
    .g_se = "Hedge's g (se)",
    condition_arm1 = "Intervention",
    condition_arm2 = "Control",
    format = "IG format",
    n_sessions_arm1 = "Number of sessions (IG)",
    baseline_n_arm1 = "N (IG)", 
    baseline_n_arm2 = "N (CG)",
    instrument= "Instrument",
    rating = "Rating",
    country = "Country",
    diagnosis = "Diagnosis",
    comorbid_mental= "Mental comorbidity",
    age_group= "Age group",
    percent_women100 = "% Female",
    recruitment = "Recruitment",
    target_group= "Target group",
    ac = "Allocation concealment",
    sg = "Sequence \n generation",
    itt= "Intention-to-treat analyses",
    ba= "Blinding of assessors",
    rob= "Overall RoB (0 = high risk to 4= low risk)")
 )


# 6 flowchart -------------------------------------------------------------
load("data/data_agg.rda")

table(data_agg$database)

# N overall RCTs included
length(unique(data_agg$study))

# N studies in adults
data_agg %>% filter(database== "adults_outpatients") %>%  pull(study) %>% unique %>% length

# N studies in children/adolescents
data_agg %>% filter(database== "child_adol") %>%  pull(study) %>% unique %>% length


