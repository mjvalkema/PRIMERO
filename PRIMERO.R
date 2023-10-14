# Code written by M.J. Valkema, March 2023

################## Setup ################## 
rm(list=ls()) # clear global environment
setwd(dir="/Users/Maartje/repos/PRIMERO")
library(tableone) # baseline tables
library(dplyr) # select
library(PairedData) # paired data comparison
library(tidyr)
library(tibble)
library(tidyselect)
library(car) # recode variables
library(ggpubr)
library(ggplot2)
library(irr) # inter-observer agreement
library(caret)

################## Load data ################## 
# Export from Castor with names(labels) instead of numbers (values)
data <- read.csv(file="PRIMERO_study_csv_export_20230214195152/PRIMERO_study_export_20230214.csv",
                         header=TRUE, sep = ";", na.strings = "")

mPETMRI <- read.csv(file="PRIMERO_study_csv_export_20230202172105/PRIMERO_study_Quantitative_measurements_PET-MRI_export_20230202.csv",
                 header=TRUE, sep = ";", na.strings = "")
mPETCT <- read.csv(file="PRIMERO_study_csv_export_20230202172105/PRIMERO_study_Quantitative_measurements_PET-CT_export_20230202.csv",
                         header=TRUE, sep = ";", na.strings = "")

measurements <- merge(mPETCT, mPETMRI, by = "Participant.Id", all = TRUE)

################## Dataprep: rewrite columns containing "stressful ################## 
# Define outcomes primary tumor area
data$outcome_conclusion_tu <- coalesce(data$specify_other_no_surgery, data$specify_no_surgery, data$surgery_cre, data$biopsies_tumor_cre4, data$biopsies_tumor_cre3)

data$outcome_conclusion_tu_final <- car::recode(data$outcome_conclusion_tu, "c('distant metastases', 'CRE-1', 'CRE-2', 'CRE-3', 'CRE-4', 'refusal')='tumor'; 
                                       c('No', 'active surveillance')='cCR'") # data$outcome_conclusion_tu[19] has cCR

data$outcome_conclusion_tu_final[15] <- "cCR" # this patient had HGD in biopsies but TRG 1 in resection specimen, therefore consider as cCR
data$outcome_conclusion_tu_final <- factor(data$outcome_conclusion_tu_final, ordered = TRUE, levels = c("tumor", "cCR"))

# LN outcomes
data$outcome_conclusion_ln <- car::recode(data$ypn,
                                          "c('N2 Metastasis in 3-6 regional lymph nodes', 'N1 Metastasis in 1-2 regional lymph nodes')='ypN+'; 
                                          c('N0 No regional lymph node metastasis')='ypN0'")
data$outcome_conclusion_ln_final <- coalesce(data$outcome_conclusion_ln, data$outcome_conclusion_tu_final)
data$outcome_conclusion_ln_final <- car::recode(data$outcome_conclusion_ln_final,
                                          "c('refusal', 'tumor')='unknown';
                                          c('cCR')='ypN0'")
# the patient noted as "tumor" has M+ and is no longer in follow-up:  LN status is therefore unknown

data$outcome_conclusion_ln[7] <- "ypN+" #LN 1R induction therapy, later curative resection


data$outcome_conclusion_ln_final <- factor(data$outcome_conclusion_ln_final, ordered = TRUE, levels = c("ypN+", "ypN0", "unknown"))

# Combining conclusion T + N
data$outcome_conclusion_TN <- paste(data$outcome_conclusion_tu_final, data$outcome_conclusion_ln_final)
data$outcome_conclusion_TN_final <-car::recode(data$outcome_conclusion_TN,
                                          "c('cCR ypN0')='ycT0N0';
                                          c('tumor ypN+', 'tumor ypN0', 'tumor unknown', 'cCR ypN+')='ypycT+/N+'")

# Dates
data$interval_ncrt_pet <- data$interval_ncrt_pet/7 # interval between end of nCRT and date of scan in weeks
data$interval_ncrt_surgery <- data$interval_ncrt_surgery/7
data$date_diagnosis <- as.Date(data$date_diagnosis, format = '%d-%m-%Y')
data$end_ncrt <- as.Date(data$end_ncrt, format = '%d-%m-%Y')
data$date_scan <- as.Date(data$date_scan, format = '%d-%m-%Y')
data$date_tumor_CRE <- as.Date(data$date_tumor_CRE, format = '%d-%m-%Y')
data$date_tumor_CRE <- replace(data$date_tumor_CRE, which(data$date_tumor_CRE == "2996-01-01"), NA) # remove N/A values
data$surgery_date <- as.Date(data$surgery_date, format = '%d-%m-%Y')
data$interval_scan_surgery <- as.numeric(difftime(data$surgery_date, data$date_scan)/7)
data$interval_tumor_surgery <- as.numeric(difftime(data$surgery_date, data$date_tumor_CRE)/7)
data$interval_ncrt_tumor <- as.numeric(difftime(data$date_tumor_CRE, data$end_ncrt)/7)

# Rearrange columns stressful
petct_stress_columnnames <- names(data)[grep('petct_stressful\\.', names(data))]
petct_stress_wide <- data[,c("Participant.Id", petct_stress_columnnames)]
petct_stress_long <- petct_stress_wide %>% 
  pivot_longer(
    cols = starts_with('petct'), 
    names_to = "petct_stressful",
    values_to = "value"
  ) %>%
  filter(value == 1) %>%
  mutate(petct_stressful = gsub('petct_stressful\\.', '', petct_stressful)) %>%
  mutate(petct_stressful = gsub('_', ' ', petct_stressful)) %>%
  summarize(Participant.Id=Participant.Id, petct_stressful=petct_stressful)

petmri_stress_columnnames <- names(data)[grep('petmri_stressful\\.', names(data))]
petmri_stress_wide <- data[,c("Participant.Id", petmri_stress_columnnames)]
petmri_stress_long <- petmri_stress_wide %>% 
  pivot_longer(
    cols = starts_with('petmri'), 
    names_to = "petmri_stressful",
    values_to = "value"
  ) %>%
  filter(value == 1) %>%
  mutate(petmri_stressful = gsub('petmri_stressful\\.', '', petmri_stressful)) %>%
  mutate(petmri_stressful = gsub('_', ' ', petmri_stressful)) %>%
  summarize(Participant.Id=Participant.Id, petmri_stressful=petmri_stressful)

data <- left_join(left_join(data, petct_stress_long, by='Participant.Id'), petmri_stress_long, by='Participant.Id') %>% 
  dplyr::select(-c(petct_stress_columnnames, petmri_stress_columnnames))

data$petct_stressful <- tolower(data$petct_stressful)
data$petmri_stressful <- tolower(data$petmri_stressful)

# make LBM variable
data$LBM <- ifelse(data$sex == "male", (1.1 * data$weight - 128 * (data$weight/data$height)^2), (1.07 * data$weight - 148 * (data$weight/data$height)^2)) # first argument for men, second for women (James equation)

measurements <- merge(measurements, data[,c('Participant.Id', 'LBM', 'weight', 'outcome_conclusion_tu_final', 'trg')], by = 'Participant.Id', all.x = TRUE)

# Calculate SULs
measurements$tu_sulmax_petct <- (measurements$tu_suvmax_petct / measurements$weight)* measurements$LBM
measurements$oes_sulmax_petct <- (measurements$oes_suvmax_petct / measurements$weight)* measurements$LBM
measurements$liv_sulmean_petct <- (measurements$liv_suvmean_petct / measurements$weight)* measurements$LBM
measurements$bp_sulmean_petct <- (measurements$bp_suvmean_petct / measurements$weight)* measurements$LBM

measurements$tu_sulmax_petmri <- (measurements$tu_suvmax_petmri / measurements$weight)* measurements$LBM
measurements$oes_sulmax_petmri <- (measurements$oes_suvmax_petmri / measurements$weight)* measurements$LBM
measurements$liv_sulmean_petmri <- (measurements$liv_suvmean_petmri / measurements$weight)* measurements$LBM
measurements$bp_sulmean_petmri <- (measurements$bp_suvmean_petmri / measurements$weight)* measurements$LBM


# Missing ADC for patient 4 (since DWI of esophagus was not acquired)
measurements$tu_adc_group1[4] <- NA
measurements$tu_adc_group2[4] <- NA

################## Results - patients and scans ################## 
# Baseline table
baseVars <- c("age", "sex", names(data)[23:26], "diabetes", "trg", "radicality", "ypt", "ypn", "prept", "prepn", "interval_ncrt_tumor", "interval_tumor_surgery", "interval_scan_surgery", "interval_ncrt_surgery")
catVars <- c("sex", "diabetes", "histology", "tumor_differentiation", "cT_stage", "cN_stage", "trg", "radicality", "ypt", "ypn", "prept", "prepn")
tableBase <- CreateTableOne(vars = baseVars, factorVars = catVars, data = data)
tableBase <- print(tableBase, nonnormal = c("age", "height", "interval_ncrt_tumor", "interval_tumor_surgery", "interval_scan_surgery", "interval_ncrt_surgery"), quote = FALSE, noSpaces = TRUE, digits=NULL)
#write.csv(tableBase, "output/Table1.csv", row.names = TRUE, na = "")

# Describe range of diagnosis
range(data$date_diagnosis)
data$activity_kg <- data$dose / data$weight
summary(data$activity_kg)

# Count of CRE-1 and CRE-2 scans
table(data$cre)

# Describe protocol violations
table(data$petct_performed)
table(data$diagnosticct_performed)
table(data$petmri_performed)
table(data$specify_nopetct)
table(data$specify_nodct)
table(data$specify_nopetmri)
table(data$protocoldev_petct)
table(data$specify_protocoldev_petct)
table(data$protocoldev_petmri)
table(data$specify_protocoldev_petmri)

# Descriptives of scans
scanVars <- c("interval_ncrt_pet", "duration_petct", "interval_inj_petct", "fasting", "prehydration", "dose", "glucose", "duration_petmri", "interval_petmri")
catVars <- c("fasting", "prehydration")
tableScan <- CreateTableOne(vars = scanVars, factorVars = catVars, data = data)
tableScan <- print(tableScan, nonnormal = c("interval_ncrt_pet", "duration_petct", 'interval_inj_petct',
                                            "glucose", "duration_petmri", "interval_petmri"), quote = FALSE, noSpaces = TRUE, digits=NULL)
#write.csv(tableScan, "output/ScanParams.csv", row.names = TRUE, na = "")

##################  Questionnaires: experience with PET-CT vs PET-MRI ################## 
petctQuestionnaires <- c("Participant.Id", "petct_uncomfortable", "petct_anxiety", "petct_painful", "petct_embarrassing")
petmriQuestionnaires <- c("Participant.Id", "petmri_uncomfortable", "petmri_anxiety", "petmri_painful", "petmri_embarrassing")

dfQuestionnaires1 <- data %>% dplyr::select(all_of(petctQuestionnaires))
dfQuestionnaires1$scan <- "petct"
names(dfQuestionnaires1) <- sub("petct_", "", names(dfQuestionnaires1)) # make names of both dataframes comparable
dfQuestionnaires2 <- data %>% dplyr::select(all_of(petmriQuestionnaires))
dfQuestionnaires2$scan <- "petmri"
names(dfQuestionnaires2) <- sub("petmri_", "", names(dfQuestionnaires2)) # make names of both dataframes comparable
dfQuestionnaires <- merge(dfQuestionnaires1, dfQuestionnaires2, by = names(dfQuestionnaires1), all = TRUE) # merge

LikertLevels <- c("absolutely not", "yes, a little", "quite a bit", "a lot", "extremely")
dfQuestionnaires$uncomfortable <- factor(dfQuestionnaires$uncomfortable, ordered = TRUE, levels = LikertLevels)
dfQuestionnaires$anxiety <- factor(dfQuestionnaires$anxiety, ordered = TRUE, levels = LikertLevels)
dfQuestionnaires$painful <- factor(dfQuestionnaires$painful, ordered = TRUE, levels = LikertLevels)
dfQuestionnaires$embarrassing <- factor(dfQuestionnaires$embarrassing, ordered = TRUE, levels = LikertLevels)

table(dfQuestionnaires$painful, dfQuestionnaires$scan) # show answer options for this question

catVars <- names(dfQuestionnaires[2:5])
tableQu <- CreateTableOne(vars = catVars, factorVars = catVars, strata = "scan", data = dfQuestionnaires)
tableQu <- print(tableQu, nonnormal = c(""), quote = FALSE, noSpaces = TRUE, digits=NULL)
#write.csv(tableQu, "output/TableBurden.csv", row.names = TRUE, na = "")
# ignore P-values, these should be tested pair-wise

# Assumption 1: are two-samples paired: yes
# Assumption 2: large sample: no, <30
# Assumption 3: normality?

dfQuestionnaires$uncomfortable_num <- as.numeric(car::recode(dfQuestionnaires$uncomfortable, "c('absolutely not')= 1 ; c('yes, a little')= 2; c('quite a bit')= 3; c('a lot')= 4"))
dfQuestionnaires$anxiety_num <- as.numeric(car::recode(dfQuestionnaires$anxiety, "c('absolutely not')= 1 ; c('yes, a little')= 2; c('quite a bit')= 3; c('a lot')= 4"))
dfQuestionnaires$painful_num <- as.numeric(car::recode(dfQuestionnaires$painful, "c('absolutely not')= 1 ; c('yes, a little')= 2; c('quite a bit')= 3; c('a lot')= 4"))
dfQuestionnaires$embarrassing_num <- as.numeric(car::recode(dfQuestionnaires$embarrassing, "c('absolutely not')= 1 ; c('yes, a little')= 2; c('quite a bit')= 3; c('a lot')= 4"))

# Summary statistics
catVars <- c("uncomfortable_num", "anxiety_num", "painful_num", "embarrassing_num")
tableQu <- CreateTableOne(vars = catVars, strata = "scan", data = dfQuestionnaires)
tableQu <- print(tableQu, nonnormal = c(""), quote = FALSE, noSpaces = TRUE, digits=NULL)
#write.csv(tableQu, "output/TableBurden_num.csv", row.names = TRUE, na = "")


# Paired data visualization
before <- subset(dfQuestionnaires, scan == "petct", uncomfortable_num, drop = TRUE)
after <- subset(dfQuestionnaires, scan == "petmri", uncomfortable_num, drop = TRUE)
pd <- paired(before, after)
plot(pd, type = "profile") + theme_bw()

before <- subset(dfQuestionnaires, scan == "petct", anxiety_num, drop = TRUE)
after <- subset(dfQuestionnaires, scan == "petmri", anxiety_num, drop = TRUE)
pd <- paired(before, after)
plot(pd, type = "profile") + theme_bw()

before <- subset(dfQuestionnaires, scan == "petct", painful_num, drop = TRUE)
after <- subset(dfQuestionnaires, scan == "petmri", painful_num, drop = TRUE)
pd <- paired(before, after)
plot(pd, type = "profile") + theme_bw()
# all values stay the same, so significance testing is NA

before <- subset(dfQuestionnaires, scan == "petct", embarrassing_num, drop = TRUE)
after <- subset(dfQuestionnaires, scan == "petmri", embarrassing_num, drop = TRUE)
pd <- paired(before, after)
plot(pd, type = "profile") + theme_bw()

d <- with(dfQuestionnaires, uncomfortable_num[scan == "petct"] - uncomfortable_num[scan == "petmri"])
shapiro.test(d) # significant, so normality cannot be assumed
d <- with(dfQuestionnaires, anxiety_num[scan == "petct"] - anxiety_num[scan == "petmri"])
shapiro.test(d) # significant, so normality cannot be assumed
d <- with(dfQuestionnaires, embarrassing_num[scan == "petct"] - embarrassing_num[scan == "petmri"])
shapiro.test(d) # significant, so normality cannot be assumed

# Wilcoxon test has to be used in all cases
wilcox.test(uncomfortable_num ~ scan, data = dfQuestionnaires, paired = TRUE, exact = FALSE)
wilcox.test(anxiety_num ~ scan, data = dfQuestionnaires, paired = TRUE, exact = FALSE)
wilcox.test(painful_num ~ scan, data = dfQuestionnaires, paired = TRUE, exact = FALSE) # all values stay the same, so significance testing is NA
wilcox.test(embarrassing_num ~ scan, data = dfQuestionnaires, paired = TRUE, exact = FALSE)

# Other parameters describing burden
vars <- c("petct_stressful", "petmri_stressful", "most_burdensome", "petmri_unpleasant", "petmri_willing_future")

data$petct_stressful <- factor(data$petct_stressful, ordered = TRUE, levels = c("non specifically", "scan duration", "insertion of the intravenous line", "body position in the scanner"))
data$petmri_stressful <- factor(data$petmri_stressful, ordered = TRUE, levels = c("non specifically", "scan duration", "noise of the scanner", "being in a small room"))
data$most_burdensome <- factor(car::recode(data$most_burdensome, "c('##USER_MISSING_96##')='none'; c('##USER_MISSING_98##')='missing'"), ordered = TRUE, levels = c("PET/CT", "PET/MRI", "none", "missing"))
data$petmri_unpleasant <- factor(data$petmri_unpleasant, ordered = TRUE, levels = c("not at all", "not so unpleasant", "unpleasant", "very unpleasant"))
data$petmri_willing_future <- factor(data$petmri_willing_future, ordered = TRUE, levels = c("absolutely not", "probably not", "neutral", "probably yes", "absolutely yes"))

tableQu <- CreateTableOne(vars = vars, factorVars = vars, includeNA = TRUE, data = data)
tableQu <- print(tableQu, nonnormal = c(""), quote = FALSE, noSpaces = TRUE, digits=NULL)
#write.csv(tableQu, "output/TableBurden_other.csv", row.names = TRUE, na = "")


################## Qualitative observations between teams ################## 
namesconclusion_qualitative1 <- names(data)[grep('conclusion', names(data))]
drop <-  names(data)[grep('clinic', names(data))] # remove clinical report conclusions
namesconclusion_qualitative <- c(namesconclusion_qualitative1[!(namesconclusion_qualitative1 %in% drop)])

# object for summarizing all qualitative scores in relation to reference standard
qualconclusions <- data[,c("Participant.Id", "date_scan", "quality_scan_petct", "quality_scan_petmri", "quality_scan_petct_team2", "quality_scan_petmri_team2", namesconclusion_qualitative,
                           "resection", "interval_ncrt_surgery", "trg", "ypt", "ypn", "prept", "prepn", "radicality")]
write.csv(qualconclusions, "output/Qualitative_conclusion.csv", row.names = FALSE, na = "")

# Quality scans
table(qualconclusions$quality_scan_petct)
table(qualconclusions$quality_scan_petmri)

table(qualconclusions$quality_scan_petct_team2)
table(qualconclusions$quality_scan_petmri_team2)

# Qualitative assessment tumor
table(qualconclusions$outcome_conclusion_tu_final)
table(qualconclusions$outcome_conclusion_ln_final)

## Qualitative conclusions (without dichotomization)
# team 1
table(qualconclusions$tu_conclusion_petct, qualconclusions$outcome_conclusion_tu_final) #petct
table(qualconclusions$tu_conclusion_petmri, qualconclusions$outcome_conclusion_tu_final) #petmri

table(qualconclusions$ln_conclusion_petct, qualconclusions$outcome_conclusion_ln_final) #petct
table(qualconclusions$ln_conclusion_petmri, qualconclusions$outcome_conclusion_ln_final) #petmri

# team 2
table(qualconclusions$tu_conclusion_petct_team2, qualconclusions$outcome_conclusion_tu_final) #petct
table(qualconclusions$tu_conclusion_petmri_team2, qualconclusions$outcome_conclusion_tu_final) #petmri

table(qualconclusions$ln_conclusion_petct_team2, qualconclusions$outcome_conclusion_ln_final) #petct
table(qualconclusions$ln_conclusion_petmri_team2, qualconclusions$outcome_conclusion_ln_final) #petmri


# Recode conclusions qualitative tumor assessment
qualconclusions$tu_conclusion_petct_final_t1 <- factor(car::recode(data$tu_conclusion_petct,
                                          "c('benign / tumor absent')='cCR'; 
                                          c('malignant / tumor present', 'equivocal')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))
qualconclusions$tu_conclusion_petct_final_t2 <- factor(car::recode(data$tu_conclusion_petct_team2,
                                                            "c('benign / tumor absent')='cCR'; 
                                          c('malignant / tumor present', 'equivocal')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))

qualconclusions$tu_conclusion_petmri_final_t1 <- factor(car::recode(data$tu_conclusion_petmri,
                                                            "c('benign / tumor absent')='cCR'; 
                                          c('malignant / tumor present', 'equivocal')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))
qualconclusions$tu_conclusion_petmri_final_t2 <- factor(car::recode(data$tu_conclusion_petmri_team2,
                                                            "c('benign / tumor absent')='cCR'; 
                                          c('malignant / tumor present', 'equivocal')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))

# Recode conclusions qualitative LN assessment
qualconclusions$ln_conclusion_petct_final_t1 <- factor(car::recode(data$ln_conclusion_petct,
                                                            "c('benign / tumor absent')='cCR'; 
                                          c('malignant / tumor present', 'equivocal')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))
qualconclusions$ln_conclusion_petct_final_t2 <- factor(car::recode(data$ln_conclusion_petct_team2,
                                                            "c('benign / tumor absent')='cCR'; 
                                          c('malignant / tumor present', 'equivocal')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))

qualconclusions$ln_conclusion_petmri_final_t1 <- factor(car::recode(data$ln_conclusion_petmri,
                                                            "c('benign / tumor absent')='cCR'; 
                                          c('malignant / tumor present', 'equivocal')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))
qualconclusions$ln_conclusion_petmri_final_t2 <- factor(car::recode(data$ln_conclusion_petmri_team2,
                                                            "c('benign / tumor absent')='cCR'; 
                                          c('malignant / tumor present', 'equivocal')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))


### Variable for integrated conclusion any residual tumor (T or N)
# recode variables Team 1
qualconclusions$overall_conclusion_petct_t1 <- paste(qualconclusions$tu_conclusion_petct_final_t1, qualconclusions$ln_conclusion_petct_final_t1)
qualconclusions$overall_conclusion_petct_final_t1 <- factor(car::recode(qualconclusions$overall_conclusion_petct_t1,
                                               "c('cCR cCR')='cCR';
                                               c('tumor cCR', 'tumor tumor', 'cCR tumor')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))

qualconclusions$overall_conclusion_petmri_t1 <- paste(qualconclusions$tu_conclusion_petmri_final_t1, qualconclusions$ln_conclusion_petmri_final_t1)
qualconclusions$overall_conclusion_petmri_final_t1 <- factor(car::recode(qualconclusions$overall_conclusion_petmri_t1,
                                                                        "c('cCR cCR')='cCR';
                                               c('tumor cCR', 'tumor tumor', 'cCR tumor')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))
# recode variables Team 2
qualconclusions$overall_conclusion_petct_t2 <- paste(qualconclusions$tu_conclusion_petct_final_t2, qualconclusions$ln_conclusion_petct_final_t2)
qualconclusions$overall_conclusion_petct_final_t2 <- factor(car::recode(qualconclusions$overall_conclusion_petct_t2,
                                                                        "c('cCR cCR')='cCR';
                                               c('tumor cCR', 'tumor tumor', 'cCR tumor')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))

qualconclusions$overall_conclusion_petmri_t2 <- paste(qualconclusions$tu_conclusion_petmri_final_t2, qualconclusions$ln_conclusion_petmri_final_t2)
qualconclusions$overall_conclusion_petmri_final_t2 <- factor(car::recode(qualconclusions$overall_conclusion_petmri_t2,
                                                                         "c('cCR cCR')='cCR';
                                               c('tumor cCR', 'tumor tumor', 'cCR tumor')='tumor'"), ordered = TRUE, levels = c("tumor", "cCR"))


## Overview qualitative conclusions dichotomized
# team 1, primary tumor
table(qualconclusions$tu_conclusion_petct_final_t1, qualconclusions$outcome_conclusion_tu_final) #petct
table(qualconclusions$tu_conclusion_petmri_final_t1, qualconclusions$outcome_conclusion_tu_final) #petmri

# team 2, primary tumor
table(qualconclusions$tu_conclusion_petct_final_t2, qualconclusions$outcome_conclusion_tu_final) #petct
table(qualconclusions$tu_conclusion_petmri_final_t2, qualconclusions$outcome_conclusion_tu_final) #petmri

# team 1, LN
table(qualconclusions$ln_conclusion_petct_final_t1, qualconclusions$outcome_conclusion_ln_final) #petct
table(qualconclusions$ln_conclusion_petmri_final_t1, qualconclusions$outcome_conclusion_ln_final) #petmri

# team 2, LN
table(qualconclusions$ln_conclusion_petct_final_t2, qualconclusions$outcome_conclusion_ln_final) #petct
table(qualconclusions$ln_conclusion_petmri_final_t2, qualconclusions$outcome_conclusion_ln_final) #petmri


# Additionaly, team 1 any tumor (T + N)
table(qualconclusions$overall_conclusion_petct_final_t1, qualconclusions$outcome_conclusion_TN_final) # petct
table(qualconclusions$overall_conclusion_petmri_final_t1, qualconclusions$outcome_conclusion_TN_final) # petmri

# Additionaly, team 2 any tumor (T + N)
table(qualconclusions$overall_conclusion_petct_final_t2, qualconclusions$outcome_conclusion_TN_final) # petct
table(qualconclusions$overall_conclusion_petmri_final_t2, qualconclusions$outcome_conclusion_TN_final) # petmri

# sens, spec, ppv and npv with 95% CI calculated with https://www.medcalc.org/calc/diagnostic_test.php

# Cohen's kappa primary tumor
irr::agree(ratings=cbind(qualconclusions$tu_conclusion_petct_final_t1, qualconclusions$tu_conclusion_petct_final_t2), tolerance = 0)
irr::agree(ratings=cbind(qualconclusions$tu_conclusion_petmri_final_t1, qualconclusions$tu_conclusion_petmri_final_t2), tolerance = 0)

irr::kappa2(ratings=cbind(qualconclusions$tu_conclusion_petct_final_t1, qualconclusions$tu_conclusion_petct_final_t2))
irr::kappa2(ratings=cbind(qualconclusions$tu_conclusion_petmri_final_t1, qualconclusions$tu_conclusion_petmri_final_t2))

# Cohen's kappa LN
irr::agree(ratings=cbind(qualconclusions$ln_conclusion_petct_final_t1, qualconclusions$ln_conclusion_petct_final_t2), tolerance = 0)
irr::agree(ratings=cbind(qualconclusions$ln_conclusion_petmri_final_t1, qualconclusions$ln_conclusion_petmri_final_t2), tolerance = 0)

irr::kappa2(ratings=cbind(qualconclusions$ln_conclusion_petct_final_t1, qualconclusions$ln_conclusion_petct_final_t2))
irr::kappa2(ratings=cbind(qualconclusions$ln_conclusion_petmri_final_t1, qualconclusions$ln_conclusion_petmri_final_t2))


################## Quantitative SUV measurements ################## 
SUL_PETCT <- c("Participant.Id", "tu_sulmax_petct", "oes_sulmax_petct", "liv_sulmean_petct", "bp_sulmean_petct")
SUL_PETMRI <- c("Participant.Id", "tu_sulmax_petmri", "oes_sulmax_petmri", "liv_sulmean_petmri", "bp_sulmean_petmri")

dfSUL_PETCT <- measurements %>% dplyr::select(all_of(SUL_PETCT))
dfSUL_PETCT$scan <- "petct"
names(dfSUL_PETCT) <- sub("_petct", "", names(dfSUL_PETCT)) # make names of both dataframes comparable
dfSUL_PETMRI <- measurements %>% dplyr::select(all_of(SUL_PETMRI))
dfSUL_PETMRI$scan <- "petmri"
names(dfSUL_PETMRI) <- sub("_petmri", "", names(dfSUL_PETMRI)) # make names of both dataframes comparable
dfSUL <- merge(dfSUL_PETCT, dfSUL_PETMRI, by = names(dfSUL_PETCT), all = TRUE) # merge
dfSUL <- merge(dfSUL, measurements[,c("Participant.Id", "outcome_conclusion_tu_final")], by = "Participant.Id", all.x = TRUE) # merge


# Check normality without log transformation
qqPlot(measurements$tu_sulmax_petct)
shapiro.test(measurements$tu_sulmax_petct) # non-significant, so normality can be assumed
qqPlot(measurements$tu_sulmax_petmri)
shapiro.test(measurements$tu_sulmax_petmri) # significant, so normality cannot be assumed !!!

qqPlot(measurements$oes_sulmax_petct)
shapiro.test(measurements$oes_sulmax_petct) # non-significant, so normality can be assumed
qqPlot(measurements$oes_sulmax_petmri)
shapiro.test(measurements$oes_sulmax_petmri) # non-significant, so normality can be assumed

qqPlot(measurements$liv_sulmean_petct)
shapiro.test(measurements$liv_sulmean_petct) # significant, so normality cannot be assumed !!!
qqPlot(measurements$liv_sulmean_petmri)
shapiro.test(measurements$liv_sulmean_petmri) # non-significant, so normality can be assumed

qqPlot(measurements$bp_sulmean_petct)
shapiro.test(measurements$bp_sulmean_petct) # significant, so normality cannot be assumed  !!!
qqPlot(measurements$bp_sulmean_petmri)
shapiro.test(measurements$bp_sulmean_petmri) # non-significant, so normality can be assumed


plot(measurements$tu_sulmax_petct, measurements$tu_sulmax_petmri)

ggscatter(measurements, x = "tu_sulmax_petct", y = "tu_sulmax_petmri", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "SULmax FDG-PET/CT", ylab = "SULmax FDG-PET/MRI")
cor.test(measurements$tu_sulmax_petct, measurements$tu_sulmax_petmri, method = "pearson") # provides CI for the correlation coefficient


# Density plots
ggplot(data=measurements, aes(x=tu_sulmax_petct, group=outcome_conclusion_tu_final, fill=outcome_conclusion_tu_final)) +
  geom_density(adjust=1.5, alpha=.4)

ggplot(data=measurements, aes(x=tu_sulmax_petmri, group=outcome_conclusion_tu_final, fill=outcome_conclusion_tu_final)) +
  geom_density(adjust=1.5, alpha=.4)


# Bland Altman plot SULmax
measurements$avg <- rowMeans(measurements[,c("tu_sulmax_petct", "tu_sulmax_petmri")])
measurements$diff <- measurements$tu_sulmax_petct - measurements$tu_sulmax_petmri

mean_diff <- mean(measurements$diff)
lower <- mean_diff - 1.96*sd(measurements$diff)
upper <- mean_diff + 1.96*sd(measurements$diff)

png(filename = "/output/figures/BlandAltmanSULmax.png", units = "cm", width=25, height=15, res=1200) # keep these sizes
fit <- (ggplot(measurements, aes(x = avg, y = diff)) +
  geom_point(size=2) +
  geom_hline(yintercept = mean_diff) +
  geom_hline(yintercept = lower, color = "red", linetype="dashed") +
  geom_hline(yintercept = upper, color = "red", linetype="dashed") +
  ggtitle("") +
  ylab("difference in SULmax measurements") +
  xlab("average of SULmax measurements")+theme_bw())
dev.off()

# Wilcoxon test for SULmax
summary(measurements$tu_sulmax_petct)
summary(measurements$tu_sulmax_petmri)
wilcox.test(tu_sulmax ~ scan, data = dfSUL, paired = TRUE, exact = FALSE)

################## Quantitative measurements overview in Table ################## 
hist(measurements$tu_adc_group1)
hist(measurements$tu_adc_group2)
hist(measurements$tu_adc_sd_group1)
hist(measurements$tu_adc_sd_group2)

qqPlot(measurements$tu_adc_group1)
shapiro.test(measurements$tu_adc_group1) # non-significant, so normality can be assumed
qqPlot(measurements$tu_adc_group2)
shapiro.test(measurements$tu_adc_group2) # non-significant, so normality can be assumed

# Quantitative measurements table
Quantvars <- c("tu_sulmax_petct", "tu_sulmax_petmri", "tu_adc_group1", "tu_adc_group2")
catVars <- c("outcome_conclusion_tu_final")
tableQuant <- CreateTableOne(vars = Quantvars, factorVars = catVars, strata = "outcome_conclusion_tu_final", data = measurements)
tableQuant <- print(tableQuant, nonnormal = c("tu_sulmax_petmri"), quote = FALSE, noSpaces = TRUE, digits=NULL)
#write.csv(tableQuant, "output/tableQuant.csv", row.names = TRUE, na = "")

# Intraclass correlation coefficient ADC
icc(measurements[, c("tu_adc_group1", "tu_adc_group2")], model = "twoway", type = "agreement", unit = "single")
