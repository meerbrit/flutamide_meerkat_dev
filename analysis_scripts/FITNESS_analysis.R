#### Analysis script for Flutamide study: Weight condition & Survival ##########
############ also includes group size descriptives #############################
############## BWalkenhorst 2023 ###############################################

#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(tidyverse)
library(brms)
library(lubridate)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(bayestestR)
library(ggdist)
library(ggeffects)
library(ggokabeito)   # Neat accessible color palette
library(tidybayes)
library(marginaleffects)
library(ggnewscale) # change colour within graph)
library(emmeans) # emtrends
library(scales) # survival percentage
library(psych) # descriptives for GS

# set working directory
setwd("G:/My Drive/Uni/UZH/Projects/KMP_PUPS/Flutamide/sheets/")

# add needed function (from external file)
source("../scripts/HighstatLib.R")

# load data
survival_data <- read_excel(
  "Flutamide_Offspring_Data_BW.xlsx", 
)

weight_data <- read_excel(
  "Flutamide_Offspring_Metadata_FULL_BW.xlsx",
)

weight_data$ID <- as.factor(weight_data$ID)
weight_data$TREATMENT <- as.factor(weight_data$TREATMENT)
weight_data$SEX <- as.factor(weight_data$SEX)
weight_data$LITTER_CODE <- as.factor(weight_data$LITTER_CODE)
weight_data$Helper_Pup_REC <- as.numeric(weight_data$Helper_Pup_REC)
weight_data$WEIGHT_DIFF_PER <- as.numeric(weight_data$WEIGHT_DIFF_PER)

survival_data$ID <- as.factor(survival_data$ID)
survival_data$LITTER_CODE <- as.factor(survival_data$LITTER_CODE)
survival_data$GROUP <- as.factor(survival_data$GROUP)

## add info about survival age
survived <- function(x) {
  return(interval(ymd(x['DOB']), ymd(x['DATE_LAST_SEEN'])) %/% days(1))
}

SURVIVED_D <- apply(survival_data, 1, survived)
SURVIVED_M <- SURVIVED_D/30
SURVIVED_01 <- SURVIVED_M >=1
SURVIVED_02 <- SURVIVED_M >=2
SURVIVED_03 <- SURVIVED_M >=3
SURVIVED_04 <- SURVIVED_M >=4
SURVIVED_06 <- SURVIVED_M >=6
SURVIVED_08 <- SURVIVED_M >=8
SURVIVED_10 <- SURVIVED_M >=10
SURVIVED_12 <- SURVIVED_M >=12
SURVIVED_14 <- SURVIVED_M >=14
SURVIVED_18 <- SURVIVED_M >=18
SURVIVED_24 <- SURVIVED_M >=24
SURVIVED_30 <- SURVIVED_M >=30
SURVIVED_36 <- SURVIVED_M >=36
survival_data <- cbind(survival_data, SURVIVED_D, SURVIVED_M, SURVIVED_01,
                       SURVIVED_02, SURVIVED_03, SURVIVED_04, SURVIVED_06,
                       SURVIVED_08, SURVIVED_10, SURVIVED_12, SURVIVED_14, 
                       SURVIVED_18, SURVIVED_24, SURVIVED_30, SURVIVED_36)
rm(SURVIVED_D, SURVIVED_M, SURVIVED_01,
   SURVIVED_02, SURVIVED_03, SURVIVED_04, SURVIVED_06,
   SURVIVED_08, SURVIVED_10, SURVIVED_12, SURVIVED_14, 
   SURVIVED_18, SURVIVED_24, SURVIVED_30, SURVIVED_36)

## remove 2nd Gen data
survival_data <- subset(survival_data, (TREATMENT == 'CTRL' | TREATMENT == 'FLUT'))
# remove P as sex
survival_data <- subset(survival_data, (SEX == 'M' | SEX == 'F'))
survival_data$TREATMENT <- as.factor(survival_data$TREATMENT)
survival_data$SEX <- as.factor(survival_data$SEX)
# recode for cens --> only needed for full survival models (Weibull) 
survival_data <- survival_data %>% 
  mutate(CURRENT_STATUS = recode(CURRENT_STATUS, 
                                 "ALIVE" = 1, 
                                 "DEAD/LAST SEEN" = 0))
survival_data$DATE_LAST_SEEN <- as.Date(survival_data$DATE_LAST_SEEN)

### check if calls have an impact on weight 
MIN_BOUTS = 3 

# select only needed data
BEG_data <- weight_data[, c(1:18, 79, 105:113, 118, 120, 122,
                            125:126, 128, 147, 162, 166:169)]
DIG_data <- weight_data[, c(19:33, 79, 105:113, 118, 120,  122, 
                          125:126, 128, 147, 162, 166:169)]
CC_data <- weight_data[, c(34:48, 79, 105:113, 118, 120,  122, 
                         125:126, 128, 147, 162, 166:169)]

BEG_data <- na.omit(BEG_data) 
DIG_data <- na.omit(DIG_data) 
CC_data <- na.omit(CC_data) 

#create the different datasets for each response variable
BEG_len_data <- subset(BEG_data, BEG_avg_Len >0 & BEG_num >= MIN_BOUTS) # call length
BEG_len_data$Avg_BEG_mil <- (BEG_len_data$BEG_avg_Len)*1000 #convert avg_BEG_LEN to milliseconds
BEG_int_data <- subset(BEG_data, BEG_avg_Int >0 & BEG_num >= MIN_BOUTS) # call interval length
BEG_int_data$Avg_BEG_mil <- (BEG_int_data$BEG_avg_Int)*1000 #convert avg_BEG_LEN to milliseconds
BEG_rate_data <- subset(BEG_data, BEG_rate >0 & BEG_num >= MIN_BOUTS) # call rate
DIG_len_data <- subset(DIG_data, DIG_avg_Len >0 & DIG_num >= MIN_BOUTS) # call length
DIG_len_data$Avg_DIG_mil <- (DIG_len_data$DIG_avg_Len)*1000 #convert avg_DIG_LEN to milliseconds
DIG_int_data <- subset(DIG_data, DIG_avg_Int >0 & DIG_num >= MIN_BOUTS) # call interval length
DIG_int_data$Avg_DIG_mil <- (DIG_int_data$DIG_avg_Int)*1000 #convert avg_DIG_LEN to milliseconds
DIG_rate_data <- subset(DIG_data, DIG_rate >0 & DIG_num >= MIN_BOUTS) # call rate
CC_len_data <- subset(CC_data, CC_avg_Len >0 & CC_num >= MIN_BOUTS) # call length
CC_len_data$Avg_CC_mil <- (CC_len_data$CC_avg_Len)*1000 #convert avg_CC_LEN to milliseconds
CC_int_data <- subset(CC_data, CC_avg_Int >0 & CC_num >= MIN_BOUTS) # call interval length
CC_int_data$Avg_CC_mil <- (CC_int_data$CC_avg_Int)*1000 #convert avg_CC_INT to milliseconds
CC_rate_data <- subset(CC_data, CC_rate >0 & CC_num >= MIN_BOUTS) # call rate

# switch working directory to models folder
setwd("G:/My Drive/Uni/UZH/Projects/KMP_PUPS/Flutamide/scripts/FITNESS_models/weight/")

# Custom ggplot theme to make pretty plots
theme_clean <- function() {
  theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_text(face = "bold"))
}

# check number of observations
table(weight_data$TREATMENT, weight_data$SEX) 

### GROUP SIZE #####
############################## ALL
# GS DOB
describeBy(weight_data$GS_all_DOB, weight_data$TREATMENT)


# GS REC
describeBy(weight_data$GS_all_REC, weight_data$TREATMENT)

# GS Emergence
describeBy(weight_data$GS_all_EMG, weight_data$TREATMENT)


#################### ADULTS
# GS DOB
describeBy(weight_data$GS_adults_DOB, weight_data$TREATMENT)

# GS REC
describeBy(weight_data$GS_adults_REC, weight_data$TREATMENT)

# GS Emergence
describeBy(weight_data$GS_adults_EMG, weight_data$TREATMENT)

############### LITTERMATES
describeBy(weight_data$GS_pups_REC_litter, weight_data$TREATMENT)

# could GS be included in models with Helpers?
corvif(sapply(weight_data[, c("REC_AGE_D", "SEX", "TREATMENT",
                                "Helper_Pup_REC", "WEIGHT_DIFF_PER",
                                "GS_all_REC", "GS_all_DOB", "GS_adults_REC", 
                                "GS_adults_DOB", "GS_pups_REC_litter"
                                )], as.numeric))

corvif(sapply(weight_data[, c("REC_AGE_D", "SEX", "TREATMENT",
                              "Helper_Pup_REC", "WEIGHT_DIFF_PER",
                              "GS_adults_DOB","GS_pups_REC_litter")], as.numeric))

corvif(sapply(weight_data[, c("REC_AGE_D", "SEX", "TREATMENT",
                              "Helper_Pup_REC", "WEIGHT_DIFF_PER",
                             "GS_pups_REC_litter")], as.numeric))


#### WEIGHT CONDITION #######
corvif(sapply(weight_data[, 
                       c("REC_AGE_D", "SEX", "TREATMENT",
                         "Helper_Pup_REC", "WEIGHT_DIFF_PER")], as.numeric))

ggqqplot(weight_data$WEIGHT_DIFF_PER)

ggplot(weight_data, aes(WEIGHT_DIFF_PER)) +
  geom_histogram(fill= "white", col= "black", binwidth = 10) +
  facet_wrap(~TREATMENT) +
  theme_bw()

ggplot(weight_data, aes(REC_AGE_D, WEIGHT_DIFF_PER)) +
  geom_point() +
  facet_grid(TREATMENT~SEX) +
  theme_bw()

#check if there is an obvious treatment effect
fit <- lm(WEIGHT_DIFF_PER~TREATMENT * REC_AGE_D * SEX, weight_data) 
summary(fit)

plot_weight <- ggplot(weight_data, aes( REC_AGE_D, WEIGHT_DIFF_PER)) +
  geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
  geom_smooth(method = "lm", col= 1) +
  labs(y= "Weight offset (%)") +
  facet_wrap(~ TREATMENT:SEX) +
  theme_bw()
annotate_figure(plot_weight, 
                top = text_grob("Weight offset (%)", 
                                color = "red", face = "bold", size = 14))
rm(plot_weight)

#### WEIGHT MODELS ####               
#### 1) Test random effects and nested level using intercept only models ####
B_WEIGHT_ran_00 <- brms::brm(formula =WEIGHT_DIFF_PER ~ 1 + (1|ID), data = weight_data, family = student(link='identity'),
                              chains = 3, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = 'B_WEIGHT_ran_00')
B_WEIGHT_ran_01 <- brms::brm(formula =WEIGHT_DIFF_PER ~ 1 + (1|LITTER_CODE), data = weight_data, family = student(link='identity'),
                              chains = 3, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = 'B_WEIGHT_ran_01')
B_WEIGHT_ran_02 <- brms::brm(formula =WEIGHT_DIFF_PER ~ 1 + (1|REC_GROUP_REF), data = weight_data, family = student(link='identity'),
                              chains = 3, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = 'B_WEIGHT_ran_02')
B_WEIGHT_ran_03 <- brms::brm(formula =WEIGHT_DIFF_PER ~ 1 + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                              chains = 3, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = 'B_WEIGHT_ran_03')
B_WEIGHT_ran_04 <- brms::brm(formula =WEIGHT_DIFF_PER ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                              chains = 3, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99), 
                              save_pars = save_pars(all=T), cores=3, file = 'B_WEIGHT_ran_04')

loo(B_WEIGHT_ran_00, B_WEIGHT_ran_01, B_WEIGHT_ran_02, B_WEIGHT_ran_03, B_WEIGHT_ran_04, moment_match = F)

performance::variance_decomposition(B_WEIGHT_ran_00)
performance::r2_bayes(B_WEIGHT_ran_00)

performance::variance_decomposition(B_WEIGHT_ran_01)
performance::r2_bayes(B_WEIGHT_ran_01)

performance::variance_decomposition(B_WEIGHT_ran_02)
performance::r2_bayes(B_WEIGHT_ran_02)

performance::variance_decomposition(B_WEIGHT_ran_03)
performance::r2_bayes(B_WEIGHT_ran_03)

performance::variance_decomposition(B_WEIGHT_ran_04)
performance::r2_bayes(B_WEIGHT_ran_04)

rm(B_WEIGHT_ran_00, B_WEIGHT_ran_01, B_WEIGHT_ran_02, B_WEIGHT_ran_03, B_WEIGHT_ran_04)

#### 2) Define all (biologically plausible) models ####
# intercept only
B_WEIGHT_00 <- brms::brm(formula =WEIGHT_DIFF_PER ~ +1 + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                          chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_00')
# one predictor
B_WEIGHT_01 <- brms::brm(formula =WEIGHT_DIFF_PER ~ TREATMENT + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                         chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                         save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01')

B_WEIGHT_01_a <- brms::brm(formula =WEIGHT_DIFF_PER ~ REC_AGE_D + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                         chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                         save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_a')

B_WEIGHT_01_b <- brms::brm(formula = WEIGHT_DIFF_PER~ SEX + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                          chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_b')

B_WEIGHT_01_c <- brms::brm(formula = WEIGHT_DIFF_PER~ Helper_Pup_REC + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                           chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_c')

# not included just checkup for impact
# B_WEIGHT_01_d <- brms::brm(formula = WEIGHT_DIFF_PER~ GS_pups_REC_litter + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
#                            chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
#                            save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_d')
# 
# B_WEIGHT_01_e <- brms::brm(formula = WEIGHT_DIFF_PER~ GS_adults_DOB + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
#                            chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
#                            save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_e')

# BEG impact ################################
# BEG rate
B_WEIGHT_01_f <- brms::brm(formula = WEIGHT_DIFF_PER~ BEG_rate + (1|LITTER_CODE/ID), data = BEG_rate_data, family = student(link='identity'),
                        chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                        save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_f')
# BEG length
B_WEIGHT_01_g <- brms::brm(formula = WEIGHT_DIFF_PER~ Avg_BEG_mil + (1|LITTER_CODE/ID), data = BEG_len_data, family = student(link='identity'),
                           chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_g')
# BEG int
B_WEIGHT_01_h <- brms::brm(formula = WEIGHT_DIFF_PER~ Avg_BEG_mil + (1|LITTER_CODE/ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_h')

# DIG impact ################################
# DIG rate
B_WEIGHT_01_i <- brms::brm(formula = WEIGHT_DIFF_PER~ DIG_rate + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_i')
# DIG length
B_WEIGHT_01_j <- brms::brm(formula = WEIGHT_DIFF_PER~ Avg_DIG_mil + (1|LITTER_CODE/ID), data = DIG_len_data, family = student(link='identity'),
                           chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_j')
# DIG int
B_WEIGHT_01_k <- brms::brm(formula = WEIGHT_DIFF_PER~ Avg_DIG_mil + (1|LITTER_CODE/ID), data = DIG_int_data, family = student(link='identity'),
                           chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_k')
# CC impact ################################
# CC rate
B_WEIGHT_01_l <- brms::brm(formula = WEIGHT_DIFF_PER~ CC_rate + (1|LITTER_CODE/ID), data = CC_rate_data, family = student(link='identity'),
                           chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_l')
# CC length
B_WEIGHT_01_m <- brms::brm(formula = WEIGHT_DIFF_PER~ Avg_CC_mil + (1|LITTER_CODE/ID), data = CC_len_data, family = student(link='identity'),
                           chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_m')
# CC int
B_WEIGHT_01_n <- brms::brm(formula = WEIGHT_DIFF_PER~ Avg_CC_mil + (1|LITTER_CODE/ID), data = CC_int_data, family = student(link='identity'),
                           chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_01_n')

#additive predictors
B_WEIGHT_02 <- brms::brm(formula =WEIGHT_DIFF_PER ~ TREATMENT + REC_AGE_D + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                         chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                         save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_02')

B_WEIGHT_03 <- brms::brm(formula =WEIGHT_DIFF_PER ~ TREATMENT + REC_AGE_D +  SEX + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                         chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                         save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_03')

B_WEIGHT_04 <- brms::brm(formula =WEIGHT_DIFF_PER ~ TREATMENT + REC_AGE_D + SEX + Helper_Pup_REC + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                         chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                         save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_04')

# interactions
B_WEIGHT_05 <- brms::brm(formula =WEIGHT_DIFF_PER ~ TREATMENT * REC_AGE_D + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                         chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                         save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_05')

B_WEIGHT_06 <- brms::brm(formula =WEIGHT_DIFF_PER ~ TREATMENT * REC_AGE_D * SEX + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                         chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                         save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_06')

B_WEIGHT_07 <- brms::brm(formula =WEIGHT_DIFF_PER ~ TREATMENT * REC_AGE_D * SEX * Helper_Pup_REC + (1|LITTER_CODE/ID), data = weight_data, family = student(link='identity'),
                         chains = 4, iter = 10000, warmup = 3000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                         save_pars = save_pars(all=T), cores=4, file = 'B_WEIGHT_07')

#### CHECK MODELS ####
summary(B_WEIGHT_00) 
plot(B_WEIGHT_00) 
pp_check(B_WEIGHT_00, ndraws = 100) 
# Pareto-k diagnostic can be useful to identify problematic point(s)
model_loo <- loo(B_WEIGHT_00, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_WEIGHT_00, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_00, type = 'acf') # autocorrelation

summary(B_WEIGHT_01) 
plot(B_WEIGHT_01) 
pp_check(B_WEIGHT_01, ndraws = 100) 
mcmc_plot(B_WEIGHT_01, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01, type = 'acf') 

summary(B_WEIGHT_01_a) 
plot(B_WEIGHT_01_a) 
pp_check(B_WEIGHT_01_a, ndraws = 100) 
mcmc_plot(B_WEIGHT_01_a, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_a, type = 'acf') 

summary(B_WEIGHT_01_b) 
plot(B_WEIGHT_01_b) #ok
pp_check(B_WEIGHT_01_b, ndraws = 100) 
mcmc_plot(B_WEIGHT_01_b, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_b, type = 'acf')

summary(B_WEIGHT_01_c) 
plot(B_WEIGHT_01_c) #ok
pp_check(B_WEIGHT_01_c, ndraws = 100) 
mcmc_plot(B_WEIGHT_01_c, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_c, type = 'acf')

summary(B_WEIGHT_01_d)
plot(B_WEIGHT_01_d)
pp_check(B_WEIGHT_01_d, ndraws = 100)
mcmc_plot(B_WEIGHT_01_d, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_d, type = 'acf')

summary(B_WEIGHT_01_e)
plot(B_WEIGHT_01_e)
pp_check(B_WEIGHT_01_e, ndraws=100)
mcmc_plot(B_WEIGHT_01_e, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_e, type = 'acf')

summary(B_WEIGHT_01_f) # BEG rate
plot(B_WEIGHT_01_f)
pp_check(B_WEIGHT_01_f, ndraws=100)
mcmc_plot(B_WEIGHT_01_f, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_f, type = 'acf')

summary(B_WEIGHT_01_g) #BEG length
plot(B_WEIGHT_01_g)
pp_check(B_WEIGHT_01_g, ndraws=100)
mcmc_plot(B_WEIGHT_01_g, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_g, type = 'acf')

summary(B_WEIGHT_01_h) # BEG int
plot(B_WEIGHT_01_h)
pp_check(B_WEIGHT_01_h, ndraws=100)
mcmc_plot(B_WEIGHT_01_h, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_h, type = 'acf')

summary(B_WEIGHT_01_i) # DIG rate
plot(B_WEIGHT_01_i)
pp_check(B_WEIGHT_01_i, ndraws=100)
mcmc_plot(B_WEIGHT_01_i, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_i, type = 'acf')

summary(B_WEIGHT_01_j) # DIG length
plot(B_WEIGHT_01_j)
pp_check(B_WEIGHT_01_j, ndraws=100)
mcmc_plot(B_WEIGHT_01_j, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_j, type = 'acf')

summary(B_WEIGHT_01_k)# DIG interval
plot(B_WEIGHT_01_k)
pp_check(B_WEIGHT_01_k, ndraws=100)
mcmc_plot(B_WEIGHT_01_k, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_k, type = 'acf')

summary(B_WEIGHT_01_l)# CC rate 
plot(B_WEIGHT_01_l)
pp_check(B_WEIGHT_01_l, ndraws=100)#super tight
mcmc_plot(B_WEIGHT_01_l, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_l, type = 'acf')

summary(B_WEIGHT_01_m) #CC length
plot(B_WEIGHT_01_m)
pp_check(B_WEIGHT_01_m, ndraws=100)
mcmc_plot(B_WEIGHT_01_m, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_m, type = 'acf')

summary(B_WEIGHT_01_n) # CC int
plot(B_WEIGHT_01_n)
pp_check(B_WEIGHT_01_n, ndraws=100)
mcmc_plot(B_WEIGHT_01_n, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_01_n, type = 'acf')

summary(B_WEIGHT_02) 
plot(B_WEIGHT_02) #ok
pp_check(B_WEIGHT_02, ndraws = 100) 
mcmc_plot(B_WEIGHT_02, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_02, type = 'acf')

summary(B_WEIGHT_03) 
plot(B_WEIGHT_03) #ok
pp_check(B_WEIGHT_03, ndraws = 100) 
mcmc_plot(B_WEIGHT_03, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_03, type = 'acf')

summary(B_WEIGHT_04) 
plot(B_WEIGHT_04) #ok
pp_check(B_WEIGHT_04, ndraws = 100)
mcmc_plot(B_WEIGHT_04, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_04, type = 'acf')

summary(B_WEIGHT_05) 
plot(B_WEIGHT_05) #ok
pp_check(B_WEIGHT_05, ndraws = 100) 
mcmc_plot(B_WEIGHT_05, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_05, type = 'acf')

summary(B_WEIGHT_06)
plot(B_WEIGHT_06) #ok
pp_check(B_WEIGHT_06, ndraws = 100) 
mcmc_plot(B_WEIGHT_06, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_06, type = 'acf')

summary(B_WEIGHT_07) 
plot(B_WEIGHT_07) #ok
pp_check(B_WEIGHT_07, ndraws = 100) 
mcmc_plot(B_WEIGHT_07, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_WEIGHT_07, type = 'acf')

rm(model_loo, df, k_rintercept)

#### MODEL COMPARISON ####
# leave-one-out-cross-validation (Vehtrai et al, 2017)
# without call models:
loo(B_WEIGHT_00, B_WEIGHT_01, B_WEIGHT_01_a, B_WEIGHT_01_b, B_WEIGHT_01_c,
    B_WEIGHT_01_d, B_WEIGHT_01_e, B_WEIGHT_02, B_WEIGHT_03, B_WEIGHT_04, 
    B_WEIGHT_05, B_WEIGHT_06, B_WEIGHT_07, moment_match=T)

#without GS info and Call models:
loo(B_WEIGHT_00, B_WEIGHT_01, B_WEIGHT_01_a, B_WEIGHT_01_b, B_WEIGHT_01_c, 
    B_WEIGHT_02, B_WEIGHT_03, B_WEIGHT_04, 
    B_WEIGHT_05, B_WEIGHT_06, B_WEIGHT_07, moment_match=T)


#### BEST (complex) WEIGHT MODEL ####
summary(B_WEIGHT_07)

diagnostic_posterior(B_WEIGHT_07)
(cred_int <- ci(B_WEIGHT_07, ci=c(.89, .95), effects=c('fixed')))

describe_posterior(
  B_WEIGHT_07,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_WEIGHT_07),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

equivalence_test(B_WEIGHT_07)

#hdi(B_WEIGHT_08, ci = 0.95)

pd <- p_direction(B_WEIGHT_07)
plot(pd)

loo_R2(B_WEIGHT_07)
loo_R2(B_WEIGHT_01_d) 
performance::variance_decomposition(B_WEIGHT_07)

##### CALCULATE ESTIMATES OF INTERACTIONS ####
# TREATMENT*AGE*HELPER 
helper_vars <- c(mean(weight_data$Helper_Pup_REC)-sd(weight_data$Helper_Pup_REC), 
                 mean(weight_data$Helper_Pup_REC), 
                 mean(weight_data$Helper_Pup_REC) + sd(weight_data$Helper_Pup_REC))


(treat_age_helper <- emtrends(B_WEIGHT_07, specs = c('TREATMENT', 'Helper_Pup_REC'), var = "REC_AGE_D",
                                at = list(Helper_Pup_REC = helper_vars)))

p_direction(treat_age_helper)

emmip(B_WEIGHT_07, TREATMENT ~ REC_AGE_D*Helper_Pup_REC, cov.reduce = range, type='response')

# plot all predictors
cond_eff <- conditional_effects(B_WEIGHT_07, re_formula = NULL, robust=TRUE)
plot(cond_eff)

### PLOTS: BEST MODEL ####
#TREATMENT*AGE*HELPER
min(weight_data$REC_AGE_D)#31
max(weight_data$REC_AGE_D)#130

helper_vars <- c(mean(weight_data$Helper_Pup_REC)-sd(weight_data$Helper_Pup_REC), 
                 mean(weight_data$Helper_Pup_REC), 
                 mean(weight_data$Helper_Pup_REC) + sd(weight_data$Helper_Pup_REC))
treatment <- c(
  "CTRL" = "Control",
  "FLUT" = "Flutamide"
)

age_dist <- B_WEIGHT_07 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(weight_data$TREATMENT),
                                    SEX = levels(weight_data$SEX),
                                    Helper_Pup_REC = helper_vars,
                                    REC_AGE_D = seq(30, 130, by=5)), 
              re_formula = NULL, allow_new_levels=T)

age_dist$WEIGHT_DIFF_PER <- age_dist$.epred
age_dist$Helper_Pup_REC <- as.factor(age_dist$Helper_Pup_REC)
age_dist$TREATMENT <- as.factor(age_dist$TREATMENT)

ggplot(age_dist, aes(x = REC_AGE_D, y = WEIGHT_DIFF_PER, color= Helper_Pup_REC, fill=Helper_Pup_REC)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = weight_data, size = 2, color='black',fill='black', alpha=0.7, shape='+') +# raw data
  scale_color_okabe_ito(name = "Helper/Pup ratio", labels=c('low', 'medium', 'high'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Helper/Pup ratio", labels=c('low', 'medium', 'high'))+
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=1)+
  geom_hline(yintercept=20, linetype="dashed", color = "green", linewidth=1)+
  geom_hline(yintercept=-20, linetype="dashed", color = "red", linewidth=1)+
  labs(x = "Age (days)", y = "Weight offset (%)") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  facet_wrap(~TREATMENT, labeller = as_labeller(treatment))+
  theme_clean()

### BEST (BASIC) MODEL ####
summary(B_WEIGHT_01_d)

describe_posterior(
  B_WEIGHT_01_d,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_WEIGHT_01_d),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

equivalence_test(B_WEIGHT_01_d)

pd <- p_direction(B_WEIGHT_01_d)
plot(pd)

loo_R2(B_WEIGHT_07)

loo_R2(B_WEIGHT_01_d) 

performance::variance_decomposition(B_WEIGHT_01_d)

cond_eff <-conditional_effects(B_WEIGHT_01_d)
plot(cond_eff)

### best basic model plots ####
min(weight_data$GS_pups_REC_litter)#1
max(weight_data$GS_pups_REC_litter)#8

gs_dist <- B_WEIGHT_01_d %>% 
  epred_draws(newdata = expand_grid(GS_pups_REC_litter = seq(1,8, by=1)), 
              re_formula = NULL, allow_new_levels=T)

gs_dist$WEIGHT_DIFF_PER <- gs_dist$.epred

ggplot(gs_dist, aes(x = GS_pups_REC_litter, y = WEIGHT_DIFF_PER)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = weight_data, size = 2, color='black',fill='black', alpha=0.7, shape='+') +# raw data
  scale_color_okabe_ito()+
  scale_fill_okabe_ito(alpha=0.2)+
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=1)+
  geom_hline(yintercept=20, linetype="dashed", color = "green", linewidth=1)+
  geom_hline(yintercept=-20, linetype="dashed", color = "red", linewidth=1)+
  labs(x = "Littersize at time of REC", y = "Weight offset (%)") +
  guides(fill='none')+
  theme_clean()

################################################################################
################################################################################
### WEIGHT - CALL models ####
# BEG calls ####
loo(B_WEIGHT_01_f, B_WEIGHT_01_g, B_WEIGHT_01_h)

########### call rate #
describe_posterior(
  B_WEIGHT_01_f,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_WEIGHT_01_f),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

equivalence_test(B_WEIGHT_01_f)

pd <- p_direction(B_WEIGHT_01_f)
plot(pd)

loo_R2(B_WEIGHT_01_f)
performance::variance_decomposition(B_WEIGHT_01_f)

# PLOTS
cond_eff <-conditional_effects(B_WEIGHT_01_f)
eff_plot <- plot(cond_eff)

min(BEG_rate_data$BEG_rate)
max(BEG_rate_data$BEG_rate)

BEG_rate_dist <- B_WEIGHT_01_f %>% 
  epred_draws(newdata = expand_grid(BEG_rate = seq(0.5, 2.5, by=0.1)), 
              re_formula = NULL, allow_new_levels=T)

BEG_rate_dist$WEIGHT_DIFF_PER <- BEG_rate_dist$.epred

ggplot(BEG_rate_dist, aes(x = BEG_rate, y = WEIGHT_DIFF_PER)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = BEG_rate_data, size = 2, color='black',fill='black', alpha=0.7, shape='+') +# raw data
  scale_color_okabe_ito()+
  scale_fill_okabe_ito(alpha=0.2)+
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=1)+
  geom_hline(yintercept=20, linetype="dashed", color = "green", linewidth=1)+
  geom_hline(yintercept=-20, linetype="dashed", color = "red", linewidth=1)+
  labs(x = "Begging call rate", y = "Weight offset (%)") +
  guides(fill='none')+
  theme_clean()

################ call length #
describe_posterior(
  B_WEIGHT_01_g,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_WEIGHT_01_g),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

equivalence_test(B_WEIGHT_01_g)

pd <- p_direction(B_WEIGHT_01_g)
plot(pd)

loo_R2(B_WEIGHT_01_g)

performance::variance_decomposition(B_WEIGHT_01_g)

### PLOTS
cond_eff <-conditional_effects(B_WEIGHT_01_g)
plot(cond_eff)

min(BEG_len_data$Avg_BEG_mil)
max(BEG_len_data$Avg_BEG_mil)

BEG_len_dist <- B_WEIGHT_01_g %>% 
  epred_draws(newdata = expand_grid(Avg_BEG_mil = seq(120, 520, by=20)), 
              re_formula = NULL, allow_new_levels=T)

BEG_len_dist$WEIGHT_DIFF_PER <- BEG_len_dist$.epred

ggplot(BEG_len_dist, aes(x = Avg_BEG_mil, y = WEIGHT_DIFF_PER)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = BEG_len_data, size = 2, color='black',fill='black', alpha=0.7, shape='+') +# raw data
  scale_color_okabe_ito()+
  scale_fill_okabe_ito(alpha=0.2)+
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=1)+
  geom_hline(yintercept=20, linetype="dashed", color = "green", linewidth=1)+
  geom_hline(yintercept=-20, linetype="dashed", color = "red", linewidth=1)+
  labs(x = "Begging call length (ms)", y = "Weight offset (%)") +
  guides(fill='none')+
  theme_clean()

################ call interval #
describe_posterior(
  B_WEIGHT_01_h,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_WEIGHT_01_h),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

equivalence_test(B_WEIGHT_01_h)

pd <- p_direction(B_WEIGHT_01_h)
plot(pd)

loo_R2(B_WEIGHT_01_h)

performance::variance_decomposition(B_WEIGHT_01_h)

#PLOTS
cond_eff <-conditional_effects(B_WEIGHT_01_h)
plot(cond_eff)

min(BEG_int_data$Avg_BEG_mil)
max(BEG_int_data$Avg_BEG_mil)

BEG_int_dist <- B_WEIGHT_01_h %>% 
  epred_draws(newdata = expand_grid(Avg_BEG_mil = seq(280, 980, by=20)), 
              re_formula = NULL, allow_new_levels=T)

BEG_int_dist$WEIGHT_DIFF_PER <- BEG_int_dist$.epred

ggplot(BEG_int_dist, aes(x = Avg_BEG_mil, y = WEIGHT_DIFF_PER)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = BEG_int_data, size = 2, color='black',fill='black', alpha=0.7, shape='+') +# raw data
  scale_color_okabe_ito()+
  scale_fill_okabe_ito(alpha=0.2)+
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=1)+
  geom_hline(yintercept=20, linetype="dashed", color = "green", linewidth=1)+
  geom_hline(yintercept=-20, linetype="dashed", color = "red", linewidth=1)+
  labs(x = "Begging call interval length (ms)", y = "Weight offset (%)") +
  guides(fill='none')+
  theme_clean()

#### DIG calls ####
########### call rate #
describe_posterior(
  B_WEIGHT_01_i,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_WEIGHT_01_i),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

equivalence_test(B_WEIGHT_01_i)

pd <- p_direction(B_WEIGHT_01_i)
plot(pd)

loo_R2(B_WEIGHT_01_i)

performance::variance_decomposition(B_WEIGHT_01_i)

cond_eff <-conditional_effects(B_WEIGHT_01_i)
eff_plot <- plot(cond_eff)

min(DIG_rate_data$DIG_rate)
max(DIG_rate_data$DIG_rate)

DIG_rate_dist <- B_WEIGHT_01_i %>% 
  epred_draws(newdata = expand_grid(DIG_rate = seq(0.5, 4.5, by=0.125)), 
              re_formula = NULL, allow_new_levels=T)

DIG_rate_dist$WEIGHT_DIFF_PER <- DIG_rate_dist$.epred

ggplot(DIG_rate_dist, aes(x = DIG_rate, y = WEIGHT_DIFF_PER)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = DIG_rate_data, size = 2, color='black',fill='black', alpha=0.7, shape='+') +# raw data
  scale_color_okabe_ito()+
  scale_fill_okabe_ito(alpha=0.2)+
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=1)+
  geom_hline(yintercept=20, linetype="dashed", color = "green", linewidth=1)+
  geom_hline(yintercept=-20, linetype="dashed", color = "red", linewidth=1)+
  labs(x = "Digging call rate", y = "Weight offset (%)") +
  guides(fill='none')+
  theme_clean()

##clean-up
rm(B_WEIGHT_00, B_WEIGHT_01, B_WEIGHT_01_a, B_WEIGHT_01_b, B_WEIGHT_01_c, B_WEIGHT_01_d,
   B_WEIGHT_01_e,B_WEIGHT_01_f, B_WEIGHT_01_g, B_WEIGHT_01_h, B_WEIGHT_01_i, B_WEIGHT_01_j, 
   B_WEIGHT_01_k, B_WEIGHT_01_l, B_WEIGHT_01_m, B_WEIGHT_01_n, B_WEIGHT_02, B_WEIGHT_03, 
   B_WEIGHT_04, B_WEIGHT_05, B_WEIGHT_06, weight_data, BEG_data, BEG_int_data, BEG_len_data,
   BEG_rate_data, DIG_data, DIG_int_data, DIG_len_data, DIG_rate_data, CC_data, 
   CC_int_data, CC_len_data, CC_rate_data, BEG_int_dist, BEG_len_dist, BEG_rate_dist,
   DIG_rate_dist, gs_dist, MIN_BOUTS)

################################################################################
################################################################################
################################################################################
################################################################################

#### SURVIVAL #### 
#################################################################
# Descriptives - overview #### 
# check number of observations
table(survival_data$TREATMENT, survival_data$SEX) # 

# Plot survival data
simple_surv_df <- survival_data[,c(7,58:64)]
# reshape the data into long format
simple_surv_df <- gather(simple_surv_df, key = "survival_time", value = "survived", 2:8)

plot_surv_data <- simple_surv_df %>%                                   
  group_by(TREATMENT, survival_time) %>%                                 
  arrange(TREATMENT) %>%                                        
  count(survived, TREATMENT) %>%                                     
  mutate(percent = n/sum(n))  

#kick out false
plot_surv_data <- subset(plot_surv_data,survived == TRUE)

plot1_survival <- 
  plot_surv_data %>%
  ggplot(data = ., mapping = aes(x = survival_time, y = percent, fill = TREATMENT)) +
  geom_col(position = 'dodge') +
  geom_text(mapping = aes(label = percent(percent)),              
            size = 3,                                             
            position = position_dodge2(width=0.5)) +             
  labs(x = "Survival stage",                                         
       y = "Percentage",                                       
       title = "Percentage of survival by treatment",        
       fill = "Treatment") +                               
  scale_y_continuous(labels = scales::percent_format()) +       
  scale_fill_okabe_ito()+
  theme_clean()

plot1_survival

# quick t-test
t.test(SURVIVED_D~TREATMENT, data=survival_data)

model_lm <-lm(SURVIVED_D~TREATMENT + SEX , data=survival_data)
summary(model_lm)
rm(model_lm)

# model for quick check
quick_survival <- brms::brm(formula = SURVIVED_D ~ TREATMENT * SEX , data=survival_data, family = gaussian(link="identity"),
          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99), 
          save_pars = save_pars(all=T), cores=4, file = 'B_Survival_AGE_Gauss')
summary(quick_survival)

pp_check(quick_survival, ndraws = 100) #ok
plot(p_direction(quick_survival))
conditional_effects(quick_survival)

quick_survival_gamma <- brms::brm(formula = SURVIVED_D ~ TREATMENT * SEX , data=survival_data, family = Gamma(link = "inverse"),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_AGE_GAMMA')
summary(quick_survival_gamma)

pp_check(quick_survival, ndraws=100)

#################################################################
#### SURVIVAL BY MONTH ####

# switch working directory to models folder
setwd("G:/My Drive/Uni/UZH/Projects/KMP_PUPS/Flutamide/scripts/FITNESS_models/survival/")

# Test group level effects
B_Survival_ran_00 <- brms::brm(formula = as.numeric(SURVIVED_04) ~ TREATMENT * SEX  + (1|LITTER_CODE) , data=survival_data, family = bernoulli(link='logit'),
                               chains = 3, iter = 6000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                               save_pars = save_pars(all=T), cores=3, file = 'B_Survival_ran_00')
B_Survival_ran_01 <- brms::brm(formula = as.numeric(SURVIVED_04) ~ TREATMENT * SEX  + (1|GROUP) , data=survival_data, family = bernoulli(link='logit'),
                               chains = 3, iter = 6000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                               save_pars = save_pars(all=T), cores=3, file = 'B_Survival_ran_01')
B_Survival_ran_02 <- brms::brm(formula = as.numeric(SURVIVED_04) ~ TREATMENT * SEX  + (1|GROUP/LITTER_CODE) , data=survival_data, family = bernoulli(link='logit'),
                               chains = 3, iter = 6000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                               save_pars = save_pars(all=T), cores=3, file = 'B_Survival_ran_02')

loo(B_Survival_ran_00, B_Survival_ran_01, B_Survival_ran_02, moment_match=F)


# Survival models by month ####
# all survived months 01, otherwise they would not be in study!
# B_Survival_M01 <- brms::brm(formula = as.numeric(SURVIVED_01) ~ TREATMENT * SEX  + (1|GROUP) , data=survival_data, family = bernoulli(link='logit'),
#                             chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
#                             save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M01')

# switch working directory to models folder
setwd("G:/My Drive/Uni/UZH/Projects/KMP_PUPS/Flutamide/scripts/FITNESS_models/survival/")

B_Survival_M02 <- brms::brm(formula = as.numeric(SURVIVED_02) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M02')
B_Survival_M03 <- brms::brm(formula = as.numeric(SURVIVED_03) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M03')
B_Survival_M04 <- brms::brm(formula = as.numeric(SURVIVED_04) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M04')
B_Survival_M06 <- brms::brm(formula = as.numeric(SURVIVED_06) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M06')
B_Survival_M08 <- brms::brm(formula = as.numeric(SURVIVED_08) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M08')
B_Survival_M10 <- brms::brm(formula = as.numeric(SURVIVED_10) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M10')
B_Survival_M12 <- brms::brm(formula = as.numeric(SURVIVED_12) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M12')
B_Survival_M14 <- brms::brm(formula = as.numeric(SURVIVED_14) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M14')
B_Survival_M18 <- brms::brm(formula = as.numeric(SURVIVED_18) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M18')
B_Survival_M24 <- brms::brm(formula = as.numeric(SURVIVED_24) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M24')
B_Survival_M30 <- brms::brm(formula = as.numeric(SURVIVED_30) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M30')
B_Survival_M36 <- brms::brm(formula = as.numeric(SURVIVED_36) ~ TREATMENT * SEX  + (1|GROUP), data=survival_data, family = bernoulli(link='logit'),
                            chains = 4, iter = 15000, warmup = 3000, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99), 
                            save_pars = save_pars(all=T), cores=4, file = 'B_Survival_M36')

# CHECK SURVIVAL (MONTHS) MODELS ####
summary(B_Survival_M02)
plot(B_Survival_M02)
pp_check(B_Survival_M02, ndraws = 100) 

summary(B_Survival_M03)
plot(B_Survival_M03)
pp_check(B_Survival_M03, ndraws = 100)

summary(B_Survival_M04) 
plot(B_Survival_M04)
pp_check(B_Survival_M04, ndraws = 100) 

summary(B_Survival_M06) ###
plot(B_Survival_M06)
pp_check(B_Survival_M06, ndraws = 100) 

summary(B_Survival_M08) # 
plot(B_Survival_M08)
pp_check(B_Survival_M08, ndraws = 100) 

summary(B_Survival_M10)
plot(B_Survival_M10)
pp_check(B_Survival_M10, ndraws = 100) 

summary(B_Survival_M12)  
plot(B_Survival_M12)
pp_check(B_Survival_M12, ndraws = 100)

summary(B_Survival_M14)
plot(B_Survival_M14)
pp_check(B_Survival_M14, ndraws = 100)

summary(B_Survival_M18)
plot(B_Survival_M18)
pp_check(B_Survival_M18, ndraws = 100) 

summary(B_Survival_M24)
plot(B_Survival_M24)
pp_check(B_Survival_M24, ndraws = 100) 

summary(B_Survival_M30)
plot(B_Survival_M30)
pp_check(B_Survival_M30, ndraws = 100) 

summary(B_Survival_M36)
plot(B_Survival_M36)
pp_check(B_Survival_M36, ndraws = 100) 

################################################################################
#### Check valid survival models ####
################################################################################

#### MONTHS: 4 #### TREATMENT
mcmc_plot(B_Survival_M04, type = 'intervals',  prob = 0.89, prob_outer=0.95)
mcmc_plot(B_Survival_M04, type='acf')
mcmc_plot(B_Survival_M04, 
         type = "areas",
         prob = 0.95)

diagnostic_posterior(B_Survival_M04)
(cred_int <- ci(B_Survival_M04, ci=c(.89, .95), effects=c('fixed')))

exp(fixef(B_Survival_M04)[,-2])

#probability
inv_logit_scaled(fixef(B_Survival_M04)[,-2])

describe_posterior(
  B_Survival_M04,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_Survival_M04),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

inv_logit_scaled(17.33)
inv_logit_scaled(0.18)
inv_logit_scaled(183.53)


equivalence_test(B_Survival_M04)

pd <- p_direction(B_Survival_M04)
plot(pd) 

bayes_R2(B_Survival_M04) 
loo_R2(B_Survival_M04)

performance::r2_loo(B_Survival_M04) 

# actual survival probabilites
(B_Survival_M04.emm <- emmeans(B_Survival_M04, pairwise~ TREATMENT))

#treatment effect = 
(treat_eff = 1- 0.9137258)
#0.0862742

p_direction(B_Survival_M04.emm)

### PLOTS ###
cond_eff <- conditional_effects(B_Survival_M04, re_formula = NULL, robust=TRUE) # include all random effects (vs NA), robust uses median
plot(cond_eff)

surv_04_dist <- B_Survival_M04 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = c('CTRL', 'FLUT'),
                                    SEX = c('F', 'M')), 
              re_formula = NULL, allow_new_levels=T)

surv_04_dist$EST <- surv_04_dist$.epred

ggplot(surv_04_dist, aes(x = TREATMENT, y = EST, color=TREATMENT)) +
  stat_boxplot(geom='errorbar', width=.8)+
  geom_boxplot(outlier.shape = NA)+
  #geom_boxplot(outlier.alpha = 0.1, outlier.colour = 'red')+
  scale_colour_okabe_ito(name="Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Treatment", y = "Probability of surviving 4 months") +
  scale_x_discrete(labels=c('Control', 'Flutamide'))+
  ylim(0.0, 1.0)+
  theme_clean()

#### MONTHS: 6 ####
mcmc_plot(B_Survival_M06, type = 'intervals',  prob = 0.89, prob_outer=0.95)
mcmc_plot(B_Survival_M06, type='acf')

diagnostic_posterior(B_Survival_M06)
(cred_int <- ci(B_Survival_M06, ci=c(.89, .95), effects=c('fixed')))

# odds ratio
exp(fixef(B_Survival_M06)[,-2])

#probability
inv_logit_scaled(fixef(B_Survival_M06)[,-2])

describe_posterior(
  B_Survival_M06,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_Survival_M06),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

equivalence_test(B_Survival_M06)

pd <- p_direction(B_Survival_M06)
plot(pd)

loo_R2(B_Survival_M06)

(B_Survival_M06.emm <- emmeans(B_Survival_M06, pairwise~ TREATMENT))

p_direction(B_Survival_M06.emm)

### PLOTS ##
cond_eff <- conditional_effects(B_Survival_M06, re_formula = NULL, robust=TRUE) # include all random effects (vs NA), robust uses median
plot(cond_eff)

surv_06_dist <- B_Survival_M06 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = c('CTRL', 'FLUT'),
                                    SEX = c('F', 'M')), 
              re_formula = NULL, allow_new_levels=T)

surv_06_dist$EST <- surv_06_dist$.epred

ggplot(surv_06_dist, aes(x = TREATMENT, y = EST, color=TREATMENT)) +
  stat_boxplot(geom='errorbar', width=.8)+
  geom_boxplot(outlier.shape = NA)+
  #geom_boxplot(outlier.alpha = 0.1, outlier.colour = 'red')+
  scale_colour_okabe_ito(name="Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Treatment", y = "Probability of surviving 6 months") +
  scale_x_discrete(labels=c('Control', 'Flutamide'))+
  ylim(0.0, 1.0)+
  theme_clean()

#### MONTHS: 8 ####
mcmc_plot(B_Survival_M08, type = 'intervals',  prob = 0.89, prob_outer=0.95)
mcmc_plot(B_Survival_M08, type='acf')

diagnostic_posterior(B_Survival_M08)
(cred_int <- ci(B_Survival_M08, ci=c(.89, .95), effects=c('fixed')))

describe_posterior(
  B_Survival_M08,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_Survival_M08),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

inv_logit_scaled(15.71)
inv_logit_scaled(0.02)
inv_logit_scaled(100.14)

equivalence_test(B_Survival_M08)

pd <- p_direction(B_Survival_M08)
plot(pd)

loo_R2(B_Survival_M08, moment_match=T)

# calculate interaction terms ####
emmip(B_Survival_M08, SEX~TREATMENT, type='response')+
  scale_colour_okabe_ito(name="Sex", labels=c('Female', 'Male'))+
  labs(x = "Treatment", y = "Probability of surviving 8 months") +
  scale_x_discrete(labels=c('Control', 'Flutamide'))+
  ylim(0.7, 1.0)+
  theme_clean()

(B_Survival_M08.emm <- emmeans(B_Survival_M08, pairwise~ TREATMENT*SEX))

p_direction(B_Survival_M08.emm)

### PLOTS ##
cond_eff <- conditional_effects(B_Survival_M08, re_formula = NULL, robust=TRUE) # include all random effects (vs NA), robust uses median
plot(cond_eff)

surv_08_dist <- B_Survival_M08 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = c('CTRL', 'FLUT'),
                                    SEX = c('F', 'M')), 
              re_formula = NULL, allow_new_levels=T)

surv_08_dist$EST <- surv_08_dist$.epred

ggplot(surv_08_dist, aes(x = TREATMENT, y = EST, color=SEX)) +
  stat_boxplot(geom='errorbar', width=.8)+
  geom_boxplot(outlier.shape = NA)+
  #geom_boxplot(outlier.alpha = 0.1, outlier.colour = 'red')+
  scale_colour_okabe_ito(name="Sex", labels=c('Female', 'Male'))+
  labs(x = "Treatment", y = "Probability of surviving 8 months") +
  scale_x_discrete(labels=c('Control', 'Flutamide'))+
  ylim(0.0, 1.0)+
  theme_clean()

# #### MONTHS: 12 ####
mcmc_plot(B_Survival_M12, type = 'intervals',  prob = 0.89, prob_outer=0.95)
mcmc_plot(B_Survival_M12, type='acf')

diagnostic_posterior(B_Survival_M12)
(cred_int <- ci(B_Survival_M12, ci=c(.89, .95), effects=c('fixed')))

 
describe_posterior(
  B_Survival_M12,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_Survival_M12),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

equivalence_test(B_Survival_M12)

pd <- p_direction(B_Survival_M12)
plot(pd)

loo_R2(B_Survival_M12, moment_match=T) 
performance::variance_decomposition(B_Survival_M12)

# calculate interaction terms ####
emmip(B_Survival_M12, SEX~TREATMENT, type='response')+
  scale_colour_okabe_ito(name="Sex", labels=c('Female', 'Male'))+
  labs(x = "Treatment", y = "Probability of surviving 12 months") +
  scale_x_discrete(labels=c('Control', 'Flutamide'))+
  ylim(0.0, 1.0)+
  theme_clean()

(B_Survival_M12.emm <- emmeans(B_Survival_M12, pairwise~ TREATMENT*SEX))

p_direction(B_Survival_M12.emm)

### PLOTS ##
cond_eff <- conditional_effects(B_Survival_M12, re_formula = NULL, robust=TRUE) # include all random effects (vs NA), robust uses median
plot(cond_eff)

surv_12_dist <- B_Survival_M12 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = c('CTRL', 'FLUT'),
                                    SEX = c('F', 'M')), 
              re_formula = NULL, allow_new_levels=T)

surv_12_dist$EST <- surv_12_dist$.epred

ggplot(surv_12_dist, aes(x = TREATMENT, y = EST, color=SEX)) +
  stat_boxplot(geom='errorbar', width=.8)+
  geom_boxplot(outlier.shape = NA)+
  #geom_boxplot(outlier.alpha = 0.1, outlier.colour = 'red')+
  scale_colour_okabe_ito(name="Sex", labels=c('Female', 'Male'))+
  labs(x = "Treatment", y = "Probability of surviving 12 months") +
  scale_x_discrete(labels=c('Control', 'Flutamide'))+
  ylim(0.0, 1.0)+
  theme_clean()

### cleanup ####
rm(B_Survival_ran_00, B_Survival_ran_01, B_Survival_ran_02, B_Survival_ran_03, 
   B_Survival_ran_04, B_Survival_M02, B_Survival_M03, B_Survival_M04, 
   B_Survival_M04_reloo, B_Survival_M06, B_Survival_M08, B_Survival_M10, 
   B_Survival_M12, B_Survival_M14, B_Survival_M18, B_Survival_M24, B_Survival_M30, 
   B_Survival_M36, surv_04_dist, surv_06_dist, surv_08_dist, surv_12_dist,
   survival_data)
