#### Analysis script for Flutamide study: DIG calls ###################
###### Bayesian Multilevel models, best model & Bayesian Model Averaging (BMA)##
############## BWalkenhorst 2022 ###############################################

#### SETUP ####

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(EnvStats)#rosnerTest for outliers
library(ggplot2) 
library(ggpubr) 
library(car)
library(tidyverse)
library(tidybayes) 
library(brms)
library(bayestestR) #e.g. diagnostic_posterior
library(sjPlot) #e.g. plot_model
library(ggokabeito)   #color palette
library(emmeans) # emtrends
library(ggnewscale) # change colour within graph
library(marginaleffects) #predictions

# set working directory
setwd("G:/My Drive/Uni/UZH/Projects/KMP_PUPS/Flutamide/")

# add needed function (from external file)
source("scripts/HighstatLib.R")

# load data
full_data <- read_excel(
  "Flutamide_Offspring_Metadata_FULL_BW.xlsx", 
)

# ensure correct format, create factors for all predictors
full_data$Code <- as.factor(full_data$Code)
full_data$REC_GROUP_REF <- as.factor(full_data$REC_GROUP_REF)
full_data$NATAL_GROUP_REF <- as.factor(full_data$NATAL_GROUP_REF)
full_data$ID <- as.factor(full_data$ID)
full_data$TREATMENT <- as.factor(full_data$TREATMENT)
full_data$SEX <- as.factor(full_data$SEX)
full_data$LITTER_CODE <- as.factor(full_data$LITTER_CODE)
full_data$Helper_Pup_REC <- as.numeric(full_data$Helper_Pup_REC)
full_data$CONDITION <- as.factor(full_data$CONDITION)
full_data$WEIGHT_DIFF <- as.numeric(full_data$WEIGHT_DIFF)
full_data$WEIGHT_DIFF_PER <- as.numeric(full_data$WEIGHT_DIFF_PER)
full_data$CONDITION_PER <- as.factor(full_data$CONDITION_PER)

# how many bouts are needed for data calculation/ included in analysis?
# eg. min 3/5 (5=maxBouts in data)
MIN_BOUTS = 3

# select only needed data
DIG_data <- full_data[, c(19:33, 79, 105:113, 118, 120,  122, 125:126, 128, 147, 162, 166:169)]
# remove NAs
DIG_data <- na.omit(DIG_data) 
# how many observations of each sex in FLUT & CTRL?
table(DIG_data$TREATMENT, DIG_data$SEX) 

#create the different datasets for each response variable
DIG_len_data <- subset(DIG_data, DIG_avg_Len >0 & DIG_num >= MIN_BOUTS) # call length
DIG_len_data$Avg_DIG_mil <- (DIG_len_data$DIG_avg_Len)*1000 #convert avg_DIG_LEN to milliseconds
DIG_int_data <- subset(DIG_data, DIG_avg_Int >0 & DIG_num >= MIN_BOUTS) # call interval length
DIG_int_data$Avg_DIG_mil <- (DIG_int_data$DIG_avg_Int)*1000 #convert avg_DIG_LEN to milliseconds
DIG_rate_data <- subset(DIG_data, DIG_rate >0 & DIG_num >= MIN_BOUTS) # call rate
DIG_prop_data <- subset(DIG_data, nBouts >= MIN_BOUTS) #call proportion
DIG_prop_data$Total_calls <- DIG_prop_data$Sum_BEG+
  DIG_prop_data$Sum_DIG+DIG_prop_data$Sum_CC+DIG_prop_data$Sum_OTHER #save the number of total calls for calc later

# switch working directory to models folder
setwd("G:/My Drive/Uni/UZH/Projects/KMP_PUPS/Flutamide/scripts/DIG_models/")

# Custom ggplot theme 
theme_clean <- function() {
  theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_text(face = "bold"))
}

#### PRECHECK ####
corvif(sapply(DIG_data[, 
                       c("REC_AGE_D", "SEX", "TREATMENT",
                         "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                         "CONDITION", "CONDITION_PER")], as.numeric))

corvif(sapply(DIG_data[, 
                       c("REC_AGE_D", "SEX", "TREATMENT",
                         "Helper_Pup_REC", "WEIGHT_DIFF_PER")], as.numeric))

################################################################################
######################## Call length ###########################################
################################################################################

#### EXPLORE DATA ####
# how many observations of each sex in FLUT & CTRL?
table(DIG_len_data$TREATMENT, DIG_len_data$SEX) #1:2 in F for sex, C balanced

# unique individuals
sample_df <- DIG_len_data[!duplicated(DIG_len_data$ID),] 
table(sample_df$SEX)
table(sample_df$TREATMENT)

# test for normality of data
ggqqplot(DIG_len_data$Avg_DIG_mil)
shapiro.test(DIG_len_data$Avg_DIG_mil) 

# use a histogram to check data distribution
ggplot(DIG_len_data, aes(Avg_DIG_mil)) +
  geom_histogram(fill= "white", col= "black", binwidth = 10) +
  facet_wrap(~TREATMENT) +
  theme_bw()
# by sex and treatment
ggplot(DIG_len_data, aes(Avg_DIG_mil)) +
  geom_histogram(fill= "white", col= "black", binwidth = 10) +
  facet_wrap(TREATMENT~SEX) +
  theme_bw()

# check for correlation of vars
cor(sapply(DIG_len_data[, c("Avg_DIG_mil", "REC_AGE_D", "SEX", "TREATMENT",
                            "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                            "CONDITION", "CONDITION_PER")], as.numeric))


# plot a simple overview for treatment/age/sex
ggplot(DIG_len_data, aes(REC_AGE_D, Avg_DIG_mil)) +
  geom_point() +
  facet_grid(TREATMENT~SEX) +
  theme_bw()

#check if there is an obvious treatment effect
plot(DIG_len_data$Avg_DIG_mil~DIG_len_data$TREATMENT)
fit <- lm(Avg_DIG_mil~TREATMENT * REC_AGE_D * SEX, DIG_len_data) 
coef(fit)
S(fit)
rm(fit)

# overview plot 
plot_DIG_len <- ggplot(DIG_len_data, aes( REC_AGE_D, Avg_DIG_mil)) +
  geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
  geom_smooth(method = "lm", col= 1) +
  labs(y= "DIG call length") +
  facet_wrap(~ TREATMENT:SEX) +
  theme_bw()
annotate_figure(plot_DIG_len, 
                top = text_grob("Avg_DIG_mil", 
                                color = "red", face = "bold", size = 14))

rm(plot_DIG_len)

####OUTLIER CHECK ####
ggplot(DIG_len_data) +
  aes(x = "", y = Avg_DIG_mil) +
  geom_boxplot(fill = "#0c4c8a") +
  scale_color_okabe_ito()+
  theme_clean()+
  facet_grid(~SEX:TREATMENT)
# get the outliers (including indices in df)
out <- boxplot.stats(DIG_len_data$Avg_DIG_mil)$out
(out_ind <- which(DIG_len_data$Avg_DIG_mil %in% c(out)))

# which datapoints seem to be outliers?
DIG_len_data[out_ind, ]
DIG_len_data[out_ind, ]$Avg_DIG_mil

# Rosner Test for outliers to double check
# k= number of suspected outliers
test <- rosnerTest(DIG_len_data$Avg_DIG_mil,
                   k = length(out))
test 

#### BAYES MODELS ##############################################################

#### 1) Test random effects and nested level
B_DIG_len_ran_00 <- brms::brm(formula =Avg_DIG_mil ~ 1 + (1|ID), data = DIG_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file="B_DIG_len_ran_00")
B_DIG_len_ran_01 <- brms::brm(formula =Avg_DIG_mil ~ 1 + (1|LITTER_CODE), data = DIG_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file="B_DIG_len_ran_01")
B_DIG_len_ran_02 <- brms::brm(formula =Avg_DIG_mil ~ 1 + (1|REC_GROUP_REF), data = DIG_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file="B_DIG_len_ran_02")
B_DIG_len_ran_03 <- brms::brm(formula =Avg_DIG_mil ~ 1 + (1|LITTER_CODE/ID), data = DIG_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file="B_DIG_len_ran_03")
B_DIG_len_ran_04 <- brms::brm(formula =Avg_DIG_mil ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = DIG_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file="B_DIG_len_ran_04")

loo(B_DIG_len_ran_00, B_DIG_len_ran_01, B_DIG_len_ran_02, B_DIG_len_ran_03, B_DIG_len_ran_04, moment_match = T)

###### ICC ####
performance::variance_decomposition(B_DIG_len_ran_00)
performance::r2_bayes(B_DIG_len_ran_00) 

performance::variance_decomposition(B_DIG_len_ran_01)
performance::r2_bayes(B_DIG_len_ran_01)

performance::variance_decomposition(B_DIG_len_ran_02)
performance::r2_bayes(B_DIG_len_ran_02) 

performance::variance_decomposition(B_DIG_len_ran_03)
performance::r2_bayes(B_DIG_len_ran_03)

performance::variance_decomposition(B_DIG_len_ran_04)
performance::r2_bayes(B_DIG_len_ran_04)

rm(B_DIG_len_ran_00,B_DIG_len_ran_01, B_DIG_len_ran_02, B_DIG_len_ran_03, B_DIG_len_ran_04)

#### 2) Define all (biologically plausible) models ####
B_DIG_len_00 <- brms::brm(formula =Avg_DIG_mil ~ +1 + (1|ID), data = DIG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_DIG_len_00")
# only one predictor
B_DIG_len_01 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT + (1|ID), data = DIG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_DIG_len_01")
# additive predictors - no interactions
B_DIG_len_02 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT +  REC_AGE_D + (1|ID), data = DIG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_DIG_len_02")
B_DIG_len_03 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT +  REC_AGE_D + SEX + (1|ID), data = DIG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_DIG_len_03")
B_DIG_len_04 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|ID), data = DIG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_DIG_len_04")
B_DIG_len_04_a <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT + REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|ID), data = DIG_len_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                            save_pars = save_pars(all=T), cores=4, file="B_DIG_len_04_a")
B_DIG_len_05 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER+  (1|ID), data = DIG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_DIG_len_05")
# interactions
B_DIG_len_06 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT *  REC_AGE_D + (1|ID), data = DIG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_DIG_len_06")
B_DIG_len_07 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT *  REC_AGE_D * SEX + (1|ID), data = DIG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_DIG_len_07")
B_DIG_len_08 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|ID), data = DIG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_DIG_len_08")
B_DIG_len_08_a <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT * REC_AGE_D * SEX  * WEIGHT_DIFF_PER + (1|ID), data = DIG_len_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                            save_pars = save_pars(all=T), cores=4, file="B_DIG_len_08_a")
B_DIG_len_09 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|ID), data = DIG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_DIG_len_09")

#### CHECK MODELS ####
summary(B_DIG_len_00) 
#posterior_summary(B_DIG_len_00)
plot(B_DIG_len_00) #ok
pp_check(B_DIG_len_00, ndraws = 100) #
# Pareto-k diagnostic can be useful to identify problematic point(s)
model_loo <- loo(B_DIG_len_00, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_01) 
#posterior_summary(B_DIG_len_01)
plot(B_DIG_len_01) #ok
pp_check(B_DIG_len_01, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_01, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_02) 
#posterior_summary(B_DIG_len_02)
plot(B_DIG_len_02) #ok
pp_check(B_DIG_len_02, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_02, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_03) 
#posterior_summary(B_DIG_len_03)
plot(B_DIG_len_03) #ok
pp_check(B_DIG_len_03, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_03, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_04) 
#posterior_summary(B_DIG_len_04)
plot(B_DIG_len_04) #ok
pp_check(B_DIG_len_04, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_04, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_04_a) 
#posterior_summary(B_DIG_len_04_a)
plot(B_DIG_len_04_a) #ok
pp_check(B_DIG_len_04_a, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_04_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_05)
#posterior_summary(B_DIG_len_05)
plot(B_DIG_len_05) #ok
pp_check(B_DIG_len_05, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_05, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_06) 
#posterior_summary(B_DIG_len_06)
plot(B_DIG_len_06) #ok
pp_check(B_DIG_len_06, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_06, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_07) 
#posterior_summary(B_DIG_len_07)
plot(B_DIG_len_07) #ok
pp_check(B_DIG_len_07, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_07, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_08) 
#posterior_summary(B_DIG_len_08)
plot(B_DIG_len_08) #ok
pp_check(B_DIG_len_08, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_08, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_08_a) 
#posterior_summary(B_DIG_len_08_a)
plot(B_DIG_len_08_a) #ok
pp_check(B_DIG_len_08_a, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_08_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_len_09) 
#posterior_summary(B_DIG_len_09)
plot(B_DIG_len_09) #ok
pp_check(B_DIG_len_09, ndraws = 100) #ok
model_loo <- loo(B_DIG_len_09, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

rm(model_loo, k_intercept, df)

#### MODEL COMPARISON ####
# leave-one-out-cross-validation (Vehtrai et al, 2017)
loo(B_DIG_len_00, B_DIG_len_01, B_DIG_len_02, B_DIG_len_03, B_DIG_len_04, 
    B_DIG_len_04_a, B_DIG_len_05, B_DIG_len_06, B_DIG_len_07, B_DIG_len_08, 
    B_DIG_len_08_a, B_DIG_len_09)

#### MOMENT MATCHING ####
loo(B_DIG_len_00, B_DIG_len_01, B_DIG_len_02, B_DIG_len_03, B_DIG_len_04, 
    B_DIG_len_04_a, B_DIG_len_05, B_DIG_len_06, B_DIG_len_07, B_DIG_len_08, 
    B_DIG_len_08_a, B_DIG_len_09, moment_match = T)

#### BEST MODEL ####
summary(B_DIG_len_04_a)

mcmc_plot(B_DIG_len_04_a, type = 'intervals')
mcmc_plot(B_DIG_len_04_a, type ='acf') #autocorrelation

pd <- p_direction(B_DIG_len_04_a)
plot(pd)

diagnostic_posterior(B_DIG_len_04_a)
(cred_int <- ci(B_DIG_len_04_a, ci=c(.89, .95), effects=c('fixed')))

describe_posterior(
  B_DIG_len_04_a,
  effects = "fixed", #'fixed vs all'
  component = "all",
  rope_range = rope_range(B_DIG_len_04_a),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

equivalence_test(B_DIG_len_04_a)

loo_R2(B_DIG_len_04_a)# adjusted R2
performance::variance_decomposition(B_DIG_len_04_a)

#### PLOTS BEST MODEL ####
# plot all predictors
cond_eff <- conditional_effects(B_DIG_len_04_a, re_formula = NULL, robust=TRUE)
plot(cond_eff)

#  AGE:
#min(DIG_len_data$REC_AGE_D)
#max(DIG_len_data$REC_AGE_D)

age_dist <- B_DIG_len_04_a %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_len_data$TREATMENT),
                                    SEX = levels(DIG_len_data$SEX),
                                    Helper_Pup_REC = mean(DIG_len_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_len_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5)), 
              re_formula = NULL, allow_new_levels=T)

age_dist$Avg_DIG_mil <- age_dist$.epred
#age only
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_DIG_mil)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = DIG_len_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Digging call length (ms)") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  guides(fill='none')+
  theme_clean()
# age & treatment
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_DIG_mil, color=TREATMENT, fill=TREATMENT)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = DIG_len_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Digging call length (ms)") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  theme_clean()

# CI crossing 0:
# treatment:
treatment_dist <- B_DIG_len_04_a %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_len_data$TREATMENT),
                                    SEX = levels(DIG_len_data$SEX),
                                    Helper_Pup_REC = mean(DIG_len_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_len_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = mean(DIG_len_data$REC_AGE_D)), 
              re_formula = NULL, allow_new_levels=T)

treatment_dist$Avg_DIG_mil <- treatment_dist$.epred

ggplot(treatment_dist, aes(x = TREATMENT, y = Avg_DIG_mil, color= TREATMENT), fill=TREATMENT) +
  stat_boxplot(geom='errorbar', width=.8)+
  geom_boxplot(outlier.shape = NA)+
  scale_color_okabe_ito(name = "Treatment", labels=c("Control", "Flutamide"))+
  scale_fill_okabe_ito(alpha=0.3)+
  geom_point(data = DIG_len_data, size = 2, aes(color=TREATMENT), alpha=0.5) +   # raw data
  labs(x = "Treatment", y = "Digging call length (ms)") +
  scale_x_discrete(labels=c('Control', 'Flutamide'))+
  theme_clean()+
  guides(fill='none')

# weight:
min(DIG_len_data$WEIGHT_DIFF_PER)#-35.90035
max(DIG_len_data$WEIGHT_DIFF_PER)#48.2253

weight_dist <- B_DIG_len_04_a %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_len_data$TREATMENT),
                                    SEX = levels(DIG_len_data$SEX),
                                    Helper_Pup_REC = mean(DIG_len_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = seq(floor(min(DIG_len_data$WEIGHT_DIFF_PER)),
                                                          ceiling(max(DIG_len_data$WEIGHT_DIFF_PER)), by=2),
                                    REC_AGE_D = mean(DIG_len_data$REC_AGE_D)), 
              re_formula = NULL, allow_new_levels=T)

weight_dist$Avg_DIG_mil <- weight_dist$.epred
#weight only
ggplot(weight_dist, aes(x = WEIGHT_DIFF_PER, y = Avg_DIG_mil)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = DIG_len_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Weight offset (%)", y = "Digging call length (ms)") +
  guides(fill='none')+
  theme_clean()
# weight & treatment
ggplot(weight_dist, aes(x = WEIGHT_DIFF_PER, y = Avg_DIG_mil, color=TREATMENT, fill=TREATMENT)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = DIG_len_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Weight offset (%)", y = "Digging call length (ms)") +
  theme_clean()

# sex:
sex_dist <- B_DIG_len_04_a %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_len_data$TREATMENT),
                                    SEX = levels(DIG_len_data$SEX),
                                    Helper_Pup_REC = mean(DIG_len_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_len_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = mean(DIG_len_data$REC_AGE_D)), 
              re_formula = NULL, allow_new_levels=T)

sex_dist$Avg_DIG_mil <- sex_dist$.epred

ggplot(sex_dist, aes(x = SEX, y = Avg_DIG_mil, color= SEX), fill=SEX) +
  stat_boxplot(geom='errorbar', width=.8)+
  geom_boxplot(outlier.shape = NA)+
  scale_color_okabe_ito(name = "Sex", labels=c("Female", "Male"))+
  scale_fill_okabe_ito(alpha=0.3)+
  geom_point(data = DIG_len_data, size = 2, aes(color=SEX), alpha=0.5) +   # raw data
  labs(x = "Sex", y = "Digging call length (ms)") +
  scale_x_discrete(labels=c('Female', 'Male'))+
  theme_clean()+
  guides(fill='none')

# AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA) ####
# no new data, just get the information of fitted values
post_avg <- posterior_average(B_DIG_len_04_a, B_DIG_len_05, B_DIG_len_02, 
                              B_DIG_len_03, B_DIG_len_04, B_DIG_len_06, 
                              B_DIG_len_07, weights = "stacking", 
                              missing = 0)

posterior_summary(post_avg)

# ci(post_avg, ci=c(.89, .95)) 

describe_posterior(
  post_avg,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_len_05),
  test = c("p_direction", "p_significance", "rope"),
 # centrality ='all',
  dispersion = TRUE
)

rm(post_avg)

## cleanup ####
rm(B_DIG_len_00, B_DIG_len_01, B_DIG_len_02, B_DIG_len_03, B_DIG_len_04, 
   B_DIG_len_04_a, B_DIG_len_05, B_DIG_len_06, B_DIG_len_07, B_DIG_len_08, 
   B_DIG_len_08_a, B_DIG_len_09, DIG_len_data)

################################################################################
######################## Call interval #########################################
################################################################################
# 1) explore data
# how many observations of each sex in FLUT & CTRL?
table(DIG_int_data$TREATMENT, DIG_int_data$SEX) 

# test for normality of data
ggqqplot(DIG_int_data$Avg_DIG_mil)
# or using Shapiro-Wilk normality test (p > 0.05 = normal distribution)
shapiro.test(DIG_int_data$Avg_DIG_mil) 

# use a histogram to check data distribution
ggplot(DIG_int_data, aes(Avg_DIG_mil)) +
  geom_histogram(fill= "white", col= "black", binwidth = 500) +
  facet_wrap(~TREATMENT) +
  theme_bw()
# by sex and treatment
ggplot(DIG_int_data, aes(Avg_DIG_mil)) +
  geom_histogram(fill= "white", col= "black", binwidth = 500) +
  facet_wrap(TREATMENT~SEX) +
  theme_bw()

# check for correlation of vars
cor(sapply(DIG_int_data[, c("Avg_DIG_mil", "REC_AGE_D", "SEX", "TREATMENT",
                            "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                            "CONDITION", "CONDITION_PER")], as.numeric))

# plot a simple overview for treatment/age/sex
ggplot(DIG_int_data, aes(REC_AGE_D, Avg_DIG_mil)) +
  geom_point() +
  facet_grid(TREATMENT~SEX) +
  theme_bw()

#check if there is an obvious treatment effect
plot(DIG_int_data$Avg_DIG_mil~DIG_int_data$TREATMENT)
fit <- lm(Avg_DIG_mil~TREATMENT * REC_AGE_D * SEX, DIG_int_data) 
coef(fit)
S(fit)
rm(fit)

# overview plot ########
plot_DIG_int <- ggplot(DIG_int_data, aes( REC_AGE_D, Avg_DIG_mil)) +
  geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
  geom_smooth(method = "lm", col= 1) +
  labs(y= "DIG call interval length") +
  facet_wrap(~ TREATMENT:SEX) +
  theme_bw()
annotate_figure(plot_DIG_int, 
                top = text_grob("DIG_avg_Int", 
                                color = "red", face = "bold", size = 14))

rm(plot_DIG_int)

####OUTLIER CHECK ####
ggplot(DIG_int_data) +
  aes(x = "", y = Avg_DIG_mil) +
  geom_boxplot(fill = "#0c4c8a") +
  scale_color_okabe_ito()+
  theme_clean()+
  facet_grid(~SEX:TREATMENT)
# get the outliers (including indices in df)
out <- boxplot.stats(DIG_int_data$Avg_DIG_mil)$out
(out_ind <- which(DIG_int_data$Avg_DIG_mil %in% c(out)))

# which datapoints seem to be outliers?
DIG_int_data[out_ind, ] 
DIG_int_data[out_ind, ]$Avg_DIG_mil

# Rosner Test for outliers to double check
# k= number of suspected outliers
test <- rosnerTest(DIG_int_data$Avg_DIG_mil,
                   k = length(out))
test 
#### BAYES MODELS ##############################################################
#### 1) Test random effects and nested level
B_DIG_int_ran_00 <- brms::brm(formula =Avg_DIG_mil ~ 1 + (1|ID), data = DIG_int_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_DIG_int_ran_00")
B_DIG_int_ran_01 <- brms::brm(formula =Avg_DIG_mil ~ 1 + (1|LITTER_CODE), data = DIG_int_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_DIG_int_ran_01")
B_DIG_int_ran_02 <- brms::brm(formula =Avg_DIG_mil ~ 1 + (1|REC_GROUP_REF), data = DIG_int_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_DIG_int_ran_02")
B_DIG_int_ran_03 <- brms::brm(formula =Avg_DIG_mil ~ 1 + (1|LITTER_CODE/ID), data = DIG_int_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_DIG_int_ran_03")
B_DIG_int_ran_04 <- brms::brm(formula =Avg_DIG_mil ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = DIG_int_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_DIG_int_ran_04")

loo(B_DIG_int_ran_00, B_DIG_int_ran_01, B_DIG_int_ran_02, B_DIG_int_ran_03, 
    B_DIG_int_ran_04, moment_match = F)

###### ICC ####
performance::variance_decomposition(B_DIG_int_ran_00)
performance::r2_bayes(B_DIG_int_ran_00) # marginal = fixed effects only


performance::variance_decomposition(B_DIG_int_ran_01)
performance::r2_bayes(B_DIG_int_ran_01)

performance::variance_decomposition(B_DIG_int_ran_02)
performance::r2_bayes(B_DIG_int_ran_02) 

performance::variance_decomposition(B_DIG_int_ran_03)
performance::r2_bayes(B_DIG_int_ran_03)

performance::variance_decomposition(B_DIG_int_ran_04)
performance::r2_bayes(B_DIG_int_ran_04)

rm(B_DIG_int_ran_00, B_DIG_int_ran_01, B_DIG_int_ran_02, B_DIG_int_ran_03, 
   B_DIG_int_ran_04)

#### 2) Define all (biologically plausible) models ####
B_DIG_int_00 <- brms::brm(formula =Avg_DIG_mil ~ +1 + (1|ID), data = DIG_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_00")
# only one predictor
B_DIG_int_01 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT + (1|ID), data = DIG_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_01")
# additive predictors - no interactions
B_DIG_int_02 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT +  REC_AGE_D + (1|ID), data = DIG_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_02")
B_DIG_int_03 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT +  REC_AGE_D + SEX + (1|ID), data = DIG_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_03")
B_DIG_int_04 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|ID), data = DIG_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_04")
B_DIG_int_04_a <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT + REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|ID), data = DIG_int_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                            save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_04_a")
B_DIG_int_05 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER+  (1|ID), data = DIG_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_05")
# interactions
B_DIG_int_06 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT *  REC_AGE_D + (1|ID), data = DIG_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_06")
B_DIG_int_07 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT *  REC_AGE_D * SEX + (1|ID), data = DIG_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_07")
B_DIG_int_08 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|ID), data = DIG_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_08")
B_DIG_int_08_a <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT * REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|ID), data = DIG_int_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                            save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_08_a")
B_DIG_int_09 <- brms::brm(formula =Avg_DIG_mil ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|ID), data = DIG_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_int_09")

#### CHECK MODELS ####
summary(B_DIG_int_00) 
#posterior_summary(B_DIG_int_00)
plot(B_DIG_int_00) #ok
pp_check(B_DIG_int_00, ndraws = 100) #ok
model_loo <- loo(B_DIG_int_00, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_01) 
#posterior_summary(B_DIG_int_01)
plot(B_DIG_int_01) #ok
pp_check(B_DIG_int_01, ndraws = 100) #
model_loo <- loo(B_DIG_int_01, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_02) 
#posterior_summary(B_DIG_int_02)
plot(B_DIG_int_02) #ok
pp_check(B_DIG_int_02, ndraws = 100) #
model_loo <- loo(B_DIG_int_02, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_03) 
#posterior_summary(B_DIG_int_03)
plot(B_DIG_int_03) #ok
pp_check(B_DIG_int_03, ndraws = 100) #
model_loo <- loo(B_DIG_int_03, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_04) 
#posterior_summary(B_DIG_int_04)
plot(B_DIG_int_04) #ok
pp_check(B_DIG_int_04, ndraws = 100) #
model_loo <- loo(B_DIG_int_04, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_04_a) 
#posterior_summary(B_DIG_int_04_a)
plot(B_DIG_int_04_a) #ok
pp_check(B_DIG_int_04_a, ndraws = 100) #
model_loo <- loo(B_DIG_int_04_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_05) 
#posterior_summary(B_DIG_int_05)
plot(B_DIG_int_05) #ok
pp_check(B_DIG_int_05, ndraws = 100) #
model_loo <- loo(B_DIG_int_05, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_06) 
#posterior_summary(B_DIG_int_06)
plot(B_DIG_int_06) #ok
pp_check(B_DIG_int_06, ndraws = 100) #
model_loo <- loo(B_DIG_int_06, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_07) 
#posterior_summary(B_DIG_int_07)
plot(B_DIG_int_07) #ok
pp_check(B_DIG_int_07, ndraws = 100) #
model_loo <- loo(B_DIG_int_07, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_08) 
#posterior_summary(B_DIG_int_08)
plot(B_DIG_int_08) #ok
pp_check(B_DIG_int_08, ndraws = 100) #
model_loo <- loo(B_DIG_int_08, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_08_a) 
#posterior_summary(B_DIG_int_08_a)
plot(B_DIG_int_08_a) #ok
pp_check(B_DIG_int_08_a, ndraws = 100) #
model_loo <- loo(B_DIG_int_08_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_int_09) 
#posterior_summary(B_DIG_int_09)
plot(B_DIG_int_09) #ok
pp_check(B_DIG_int_09, ndraws = 100) #
model_loo <- loo(B_DIG_int_09, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

rm(k_rintercept, df, model_loo)

#### MODEL COMPARISON ####
# leave-one-out-cross-validation (Vehtrai et al, 2017)
loo(B_DIG_int_00, B_DIG_int_01, B_DIG_int_02, B_DIG_int_03, B_DIG_int_04, 
    B_DIG_int_04_a, B_DIG_int_05, B_DIG_int_06, B_DIG_int_07, B_DIG_int_08, 
    B_DIG_int_08_a, B_DIG_int_09)

#### MOMENT MATCHING ####
loo(B_DIG_int_00, B_DIG_int_01, B_DIG_int_02, B_DIG_int_03, B_DIG_int_04,
    B_DIG_int_04_a, B_DIG_int_05, B_DIG_int_06, B_DIG_int_07, B_DIG_int_08,
    B_DIG_int_08_a, B_DIG_int_09, moment_match=T)

#### BEST MODEL ####
summary(B_DIG_int_06)

mcmc_plot(B_DIG_int_06, type = 'intervals')
mcmc_plot(B_DIG_int_06, type = 'acf')

pd <- p_direction(B_DIG_int_06)
plot(pd)

diagnostic_posterior(B_DIG_int_06)
(cred_int <- ci(B_DIG_int_06, ci=c(.89, .95), effects=c('fixed')))

describe_posterior(
  B_DIG_int_06,
  effects = "fixed", #'all' for random effects!
  component = "all",
  rope_range = rope_range(B_DIG_int_06),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

equivalence_test(B_DIG_int_06)

loo_R2(B_DIG_int_06) # loo adjusted R2

performance::variance_decomposition(B_DIG_int_06)

####PLOTS: BEST MODEL ####
# plot all predictors
cond_eff <- conditional_effects(B_DIG_int_06, re_formula = NULL,  robust=TRUE)
plot(cond_eff)

#AGE:
min(DIG_int_data$REC_AGE_D)
max(DIG_int_data$REC_AGE_D)

age_dist <- B_DIG_int_06 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_int_data$TREATMENT),
                                    SEX = levels(DIG_int_data$SEX),
                                    Helper_Pup_REC = mean(DIG_int_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_int_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5)), 
              re_formula = NULL, allow_new_levels=T)

age_dist$Avg_DIG_mil <- age_dist$.epred
#age only
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_DIG_mil)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = DIG_int_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Digging call interval length (ms)") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  guides(fill='none')+
  theme_clean()
# age & treatment
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_DIG_mil, color=TREATMENT, fill=TREATMENT)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = DIG_int_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Digging call interval length (ms)") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  theme_clean()

# CI crossing 0:
# TREATMENT:
treatment_dist <- B_DIG_int_06 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_int_data$TREATMENT),
                                    SEX = levels(DIG_int_data$SEX),
                                    Helper_Pup_REC = mean(DIG_int_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_int_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = mean(DIG_int_data$REC_AGE_D)), 
              re_formula = NULL, allow_new_levels=T)

treatment_dist$Avg_DIG_mil <- treatment_dist$.epred

ggplot(treatment_dist, aes(x = TREATMENT, y = Avg_DIG_mil, color= TREATMENT), fill=TREATMENT) +
  stat_boxplot(geom='errorbar', width=.8)+
  geom_boxplot(outlier.shape = NA)+
  scale_color_okabe_ito(name = "Treatment", labels=c("Control", "Flutamide"))+
  scale_fill_okabe_ito(alpha=0.3)+
  geom_point(data = DIG_int_data, size = 2, aes(color=TREATMENT), alpha=0.5) +   # raw data
  labs(x = "Treatment", y = "Digging call interval length (ms)") +
  scale_x_discrete(labels=c('Control', 'Flutamide'))+
  theme_clean()+
  guides(color='none', fill='none')
  
################################################################################
# AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA) ####
# no new data, just get the information of fitted values
post_avg <- posterior_average(B_DIG_int_06, B_DIG_int_01, B_DIG_int_02, B_DIG_int_03, 
                              B_DIG_int_04, B_DIG_int_04_a, B_DIG_int_05, B_DIG_int_07, 
                              weights = "stacking", 
                              missing = 0)

posterior_summary(post_avg)

# ci(post_avg, ci=c(.89, .95)) 

describe_posterior(
  post_avg,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_int_06),
  test = c("p_direction", "p_significance", "rope"),
#  centrality ='all',
  dispersion = TRUE
)

rm(post_avg)

## cleanup ####
rm(B_DIG_int_00, B_DIG_int_01, B_DIG_int_02, B_DIG_int_03, B_DIG_int_04, 
   B_DIG_int_04_a, B_DIG_int_05, B_DIG_int_06, B_DIG_int_07, B_DIG_int_08, 
   B_DIG_int_08_a, B_DIG_int_09, DIG_int_data)

################################################################################
######################## Call rate #############################################
################################################################################
# 1) explore data
# how many observations of each sex in FLUT & CTRL?
table(DIG_rate_data$TREATMENT, DIG_rate_data$SEX) #

# test for normality of data
ggqqplot(DIG_rate_data$DIG_rate)
shapiro.test(DIG_rate_data$DIG_rate) # 0.805

# use a histogram to check data distribution
ggplot(DIG_rate_data, aes(DIG_rate)) +
  geom_histogram(fill= "white", col= "black", binwidth = 0.5) +
  facet_wrap(~TREATMENT) +
  theme_bw()
# by sex and treatment
ggplot(DIG_rate_data, aes(DIG_rate)) +
  geom_histogram(fill= "white", col= "black", binwidth = 0.5) +
  facet_wrap(TREATMENT~SEX) +
  theme_bw()

# check for correlation of vars
cor(sapply(DIG_rate_data[, c("DIG_rate", "REC_AGE_D", "SEX", "TREATMENT",
                             "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                             "CONDITION", "CONDITION_PER")], as.numeric))

# plot a simple overview for treatment/age/sex
ggplot(DIG_rate_data, aes(REC_AGE_D, DIG_rate)) +
  geom_point() +
  facet_grid(TREATMENT~SEX) +
  theme_bw()

#check if there is an obvious treatment effect
plot(DIG_rate_data$DIG_rate~DIG_rate_data$TREATMENT)
fit <- lm(DIG_rate~TREATMENT * REC_AGE_D * SEX, DIG_rate_data) 
coef(fit)
S(fit)
rm(fit)

# overview plot ########
plot_DIG_rate <- ggplot(DIG_rate_data, aes( REC_AGE_D, DIG_rate)) +
  geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
  geom_smooth(method = "lm", col= 1) +
  labs(y= "DIG call rate") +
  facet_wrap(~ TREATMENT:SEX) +
  theme_bw()
annotate_figure(plot_DIG_rate, 
                top = text_grob("DIG_rate", 
                                color = "red", face = "bold", size = 14))

rm(plot_DIG_rate)

####OUTLIER CHECK ####
ggplot(DIG_rate_data) +
  aes(x = "", y = DIG_rate) +
  geom_boxplot(fill = "#0c4c8a") +
  scale_color_okabe_ito()+
  theme_clean()+
  facet_grid(~SEX:TREATMENT)
# get the outliers (including indices in df)
out <- boxplot.stats(DIG_rate_data$DIG_rate)$out
(out_ind <- which(DIG_rate_data$DIG_rate %in% c(out)))

# which datapoints seem to be outliers?
DIG_rate_data[out_ind, ]
DIG_rate_data[out_ind, ]$DIG_rate

# Rosner Test for outliers to double check
# k= number of suspected outliers
test <- rosnerTest(DIG_rate_data$DIG_rate,
                   k = length(out))
test

#### BAYES MODELS ##############################################################
#### 1) Test random effects and nested level
B_DIG_rat_ran_00 <- brms::brm(formula =DIG_rate ~ 1 + (1|ID), data = DIG_rate_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_DIG_rat_ran_00")
B_DIG_rat_ran_01 <- brms::brm(formula =DIG_rate ~ 1 + (1|LITTER_CODE), data = DIG_rate_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_DIG_rat_ran_01")
B_DIG_rat_ran_02 <- brms::brm(formula =DIG_rate ~ 1 + (1|REC_GROUP_REF), data = DIG_rate_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_DIG_rat_ran_02")
B_DIG_rat_ran_03 <- brms::brm(formula =DIG_rate ~ 1 + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_DIG_rat_ran_03")
B_DIG_rat_ran_04 <- brms::brm(formula =DIG_rate ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_DIG_rat_ran_04")

loo(B_DIG_rat_ran_00, B_DIG_rat_ran_01, B_DIG_rat_ran_02, B_DIG_rat_ran_03, B_DIG_rat_ran_04, moment_match = F)

###### ICC ####
performance::variance_decomposition(B_DIG_rat_ran_00)
performance::r2_bayes(B_DIG_rat_ran_00) # marginal = fixed effects only

performance::variance_decomposition(B_DIG_rat_ran_01)
performance::r2_bayes(B_DIG_rat_ran_01)

performance::variance_decomposition(B_DIG_rat_ran_02)
performance::r2_bayes(B_DIG_rat_ran_02) 

performance::variance_decomposition(B_DIG_rat_ran_03)
performance::r2_bayes(B_DIG_rat_ran_03)

performance::variance_decomposition(B_DIG_rat_ran_04)
performance::r2_bayes(B_DIG_rat_ran_04)

rm(B_DIG_rat_ran_00, B_DIG_rat_ran_01, B_DIG_rat_ran_02, B_DIG_rat_ran_03, 
   B_DIG_rat_ran_04)

#### 2) Define all (biologically plausible) models ####
B_DIG_rat_00 <- brms::brm(formula = DIG_rate~ +1 + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_00")
# only one predictor
B_DIG_rat_01 <- brms::brm(formula =DIG_rate ~ TREATMENT + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_01")
# additive predictors - no interactions
B_DIG_rat_02 <- brms::brm(formula =DIG_rate ~ TREATMENT +  REC_AGE_D + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_02")
B_DIG_rat_03 <- brms::brm(formula =DIG_rate ~ TREATMENT +  REC_AGE_D + SEX + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                          chains = 3, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=3, file = "B_DIG_rat_03")
B_DIG_rat_04 <- brms::brm(formula =DIG_rate ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_04")
B_DIG_rat_04_a <- brms::brm(formula =DIG_rate ~ TREATMENT + REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_04_a")
B_DIG_rat_05 <- brms::brm(formula =DIG_rate ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER+  (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_05")
# interactions
B_DIG_rat_06 <- brms::brm(formula =DIG_rate ~ TREATMENT *  REC_AGE_D + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_06")
B_DIG_rat_07 <- brms::brm(formula =DIG_rate ~ TREATMENT *  REC_AGE_D * SEX + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_07")
B_DIG_rat_08 <- brms::brm(formula =DIG_rate ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_08")
B_DIG_rat_08_a <- brms::brm(formula =DIG_rate ~ TREATMENT * REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_08_a")
B_DIG_rat_09 <- brms::brm(formula =DIG_rate ~ TREATMENT * REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|LITTER_CODE/ID), data = DIG_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_DIG_rat_09")

#### CHECK MODELS ####
summary(B_DIG_rat_00) 
#posterior_summary(B_DIG_rat_00)
plot(B_DIG_rat_00) #ok
pp_check(B_DIG_rat_00, ndraws = 100) #ok
model_loo <- loo(B_DIG_rat_00, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_01) 
#posterior_summary(B_DIG_rat_01)
plot(B_DIG_rat_01) #ok
pp_check(B_DIG_rat_01, ndraws = 100) #
model_loo <- loo(B_DIG_rat_01, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_02) 
#posterior_summary(B_DIG_rat_02)
plot(B_DIG_rat_02) #ok
pp_check(B_DIG_rat_02, ndraws = 100) #
model_loo <- loo(B_DIG_rat_02, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_03) 
#posterior_summary(B_DIG_rat_03)
plot(B_DIG_rat_03) #ok
pp_check(B_DIG_rat_03, ndraws = 100) #
model_loo <- loo(B_DIG_rat_03, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_04) 
#posterior_summary(B_DIG_rat_04)
plot(B_DIG_rat_04) #ok
pp_check(B_DIG_rat_04, ndraws = 100) #
model_loo <- loo(B_DIG_rat_04, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_04_a) 
#posterior_summary(B_DIG_rat_04_a)
plot(B_DIG_rat_04_a) #ok
pp_check(B_DIG_rat_04_a, ndraws = 100) #
model_loo <- loo(B_DIG_rat_04_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_05) 
#posterior_summary(B_DIG_rat_05)
plot(B_DIG_rat_05) #ok
pp_check(B_DIG_rat_05, ndraws = 100) #
model_loo <- loo(B_DIG_rat_05, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_06) 
#posterior_summary(B_DIG_rat_06)
plot(B_DIG_rat_06) #ok
pp_check(B_DIG_rat_06, ndraws = 100) #
model_loo <- loo(B_DIG_rat_06, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_07) 
#posterior_summary(B_DIG_rat_07)
plot(B_DIG_rat_07) #ok
pp_check(B_DIG_rat_07, ndraws = 100) #
model_loo <- loo(B_DIG_rat_07, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_08) 
#posterior_summary(B_DIG_rat_08)
plot(B_DIG_rat_08) #ok
pp_check(B_DIG_rat_08, ndraws = 100) #
model_loo <- loo(B_DIG_rat_08, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_08_a) 
#posterior_summary(B_DIG_rat_08_a)
plot(B_DIG_rat_08_a) #ok
pp_check(B_DIG_rat_08_a, ndraws = 100) #
model_loo <- loo(B_DIG_rat_08_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_rat_09) 
#posterior_summary(B_DIG_rat_09)
plot(B_DIG_rat_09) #ok
pp_check(B_DIG_rat_09, ndraws = 100) #
model_loo <- loo(B_DIG_rat_09, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

rm(k_rintercept, df, model_loo)

#### MODEL COMPARISON ####
# leave-one-out-cross-validation (Vehtrai et al, 2017)
loo(B_DIG_rat_00, B_DIG_rat_01, B_DIG_rat_02, B_DIG_rat_03, B_DIG_rat_04, 
    B_DIG_rat_04_a, B_DIG_rat_05, B_DIG_rat_06, B_DIG_rat_07, B_DIG_rat_08, 
    B_DIG_rat_08_a, B_DIG_rat_09)

#### MOMENT MATCHING ####
loo(B_DIG_rat_00, B_DIG_rat_01, B_DIG_rat_02, B_DIG_rat_03, B_DIG_rat_04, 
    B_DIG_rat_04_a, B_DIG_rat_05, B_DIG_rat_06, B_DIG_rat_07, B_DIG_rat_08, 
    B_DIG_rat_08_a, B_DIG_rat_09, moment_match = T)

#### BEST MODEL ####
summary(B_DIG_rat_06)

mcmc_plot(B_DIG_rat_06, type = 'intervals')
mcmc_plot(B_DIG_rat_06, type = 'acf')

pd <- p_direction(B_DIG_rat_06)
plot(pd)

diagnostic_posterior(B_DIG_rat_06)
(cred_int <- ci(B_DIG_rat_06, ci=c(.89, .95), effects=c('fixed')))

describe_posterior(
  B_DIG_rat_06,
  effects = "fixed", #all
  component = "all",
  rope_range= rope_range(B_DIG_rat_06),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
  )

equivalence_test(B_DIG_rat_06)

loo_R2(B_DIG_rat_06)# loo adjusted
performance::variance_decomposition(B_DIG_rat_06)

## calculate estimates for INTERACTION ####
(treat_age <- emtrends(B_DIG_rat_06, pairwise ~ TREATMENT, var="REC_AGE_D"))

emmip(B_DIG_rat_06, TREATMENT ~ REC_AGE_D, cov.reduce = range)
p_direction(treat_age)

#### PLOTS_ BEST MODEL ####
# plot all predictors
cond_eff <- conditional_effects(B_DIG_rat_06, re_formula = NULL, robust=TRUE) 
plot(cond_eff)

#TREATMENT*AGE: CI at 0!!
min(DIG_rate_data$REC_AGE_D)#34
max(DIG_rate_data$REC_AGE_D)#130

age_dist <- B_DIG_rat_06 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_rate_data$TREATMENT),
                                    SEX = levels(DIG_rate_data$SEX),
                                    Helper_Pup_REC = mean(DIG_rate_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_rate_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5)), 
              re_formula = NULL, allow_new_levels=T)

age_dist$DIG_rate <- age_dist$.epred
# age & treatment
ggplot(age_dist, aes(x = REC_AGE_D, y = DIG_rate, color=TREATMENT, fill=TREATMENT)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = DIG_rate_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Digging call rate") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  theme_clean()

# TREATMENT:
treatment_dist <- B_DIG_rat_06 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_rate_data$TREATMENT),
                                    SEX = levels(DIG_rate_data$SEX),
                                    Helper_Pup_REC = mean(DIG_rate_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_rate_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = mean(DIG_rate_data$REC_AGE_D)), 
              re_formula = NULL, allow_new_levels=T)

treatment_dist$DIG_rate <- treatment_dist$.epred

ggplot(treatment_dist, aes(x = TREATMENT, y = DIG_rate, color= TREATMENT), fill=TREATMENT) +
  stat_boxplot(geom='errorbar', width=.8)+
  geom_boxplot(outlier.shape = NA)+
  scale_color_okabe_ito(name = "Treatment", labels=c("Control", "Flutamide"))+
  scale_fill_okabe_ito(alpha=0.3)+
  geom_point(data = DIG_rate_data, size = 2, aes(color=TREATMENT), alpha=0.5) +   # raw data
  labs(x = "Treatment", y = "Digging call rate") +
  scale_x_discrete(labels=c('Control', 'Flutamide'))+
  theme_clean()+
  guides(color='none', fill='none')

#AGE:
min(DIG_rate_data$REC_AGE_D)#34
max(DIG_rate_data$REC_AGE_D)#130

age_dist <- B_DIG_rat_06 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_rate_data$TREATMENT),
                                    SEX = levels(DIG_rate_data$SEX),
                                    Helper_Pup_REC = mean(DIG_rate_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_rate_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5)), 
              re_formula = NULL, allow_new_levels=T)

age_dist$DIG_rate <- age_dist$.epred
# age
ggplot(age_dist, aes(x = REC_AGE_D, y = DIG_rate)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = DIG_rate_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Digging call rate") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  guides(fill='none')+
  theme_clean()

# AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA) ####
# no new data, just get the information of fitted values
post_avg <- posterior_average(B_DIG_rat_06, B_DIG_rat_07, B_DIG_rat_08_a,
                              B_DIG_rat_08, weights = "stacking", 
                              missing = 0)

posterior_summary(post_avg)

# ci(post_avg, ci=c(.89, .95)) 

describe_posterior(
  post_avg,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_rat_06),
  test = c("p_direction", "p_significance", "rope"),
  #centrality = 'all',
  dispersion = TRUE
)

rm(post_avg)

## cleanup ####
rm(B_DIG_rat_00, B_DIG_rat_01, B_DIG_rat_02, B_DIG_rat_03, B_DIG_rat_04, 
   B_DIG_rat_04_a, B_DIG_rat_05, B_DIG_rat_06, B_DIG_rat_07, B_DIG_rat_08, 
   B_DIG_rat_08_a, B_DIG_rat_09, DIG_rate_data)

################################################################################
######################## Call proportion #######################################
################################################################################

#### EXPLORE DATA ####
# how many observations of each sex in FLUT & CTRL?
table(DIG_prop_data$TREATMENT, DIG_prop_data$SEX) #

# unique individuals
sample_df <- DIG_prop_data[!duplicated(DIG_prop_data$ID),]
table(sample_df$SEX)
table(sample_df$TREATMENT)

# test for normality of data
# visually using qq-plots 
ggqqplot(DIG_prop_data$Avg_DIG_prop)
# or using Shapiro-Wilk normality test (p > 0.05 = normal distribution)
shapiro.test(DIG_prop_data$Avg_DIG_prop)

# use a histogram to check data distribution
ggplot(DIG_prop_data, aes(Avg_DIG_prop)) +
  geom_histogram(fill= "white", col= "black", binwidth = 0.1) +
  facet_wrap(~TREATMENT) +
  theme_bw()
# by sex and treatment
ggplot(DIG_prop_data, aes(Avg_DIG_prop)) +
  geom_histogram(fill= "white", col= "black", binwidth = 0.1) +
  facet_wrap(TREATMENT~SEX) +
  theme_bw()

# check for correlation of vars
cor(sapply(DIG_prop_data[, c("Avg_DIG_prop", "REC_AGE_D", "SEX", "TREATMENT",
                             "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                             "CONDITION", "CONDITION_PER")], as.numeric))

# plot a simple overview for treatment/age/sex
ggplot(DIG_prop_data, aes(REC_AGE_D, Avg_DIG_prop)) +
  geom_point() +
  facet_grid(TREATMENT~SEX) +
  theme_bw()

#check if there is an obvious treatment effect
plot(DIG_prop_data$Avg_DIG_prop~DIG_prop_data$TREATMENT)
fit <- lm(Avg_DIG_prop~TREATMENT * REC_AGE_D * SEX, DIG_prop_data) 
coef(fit)
S(fit)
rm(fit)

# overview plot 
plot_DIG_prop <- ggplot(DIG_prop_data, aes( REC_AGE_D, Avg_DIG_prop)) +
  geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
  geom_smooth(method = "lm", col= 1) +
  labs(y= "DIG call length") +
  facet_wrap(~ TREATMENT:SEX) +
  theme_bw()
annotate_figure(plot_DIG_prop, 
                top = text_grob("Avg_DIG_prop", 
                                color = "red", face = "bold", size = 14))

rm(plot_DIG_prop)

####OUTLIER CHECK ####
ggplot(DIG_prop_data) +
  aes(x = "", y = Avg_DIG_prop) +
  geom_boxplot(fill = "#0c4c8a") +
  scale_color_okabe_ito()+
  theme_clean()+
  facet_grid(~SEX:TREATMENT)
# get the outliers (including indices in df)
out <- boxplot.stats(DIG_prop_data$Avg_DIG_prop)$out
(out_ind <- which(DIG_prop_data$Avg_DIG_prop %in% c(out)))

# which datapoints seem to be outliers?
DIG_prop_data[out_ind, ]
DIG_prop_data[out_ind, ]$Avg_DIG_prop

# Rosner Test for outliers to double check
# k= number of suspected outliers
test <- rosnerTest(DIG_prop_data$Avg_DIG_prop,
                   k = length(out))
test 

#### BAYES MODELS ##############################################################
#### 1) Test random effects and nested level
B_DIG_prop_ran_00 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ 1 + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_DIG_prop_ran_00")
B_DIG_prop_ran_01 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ 1 + (1|LITTER_CODE), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_DIG_prop_ran_01")
B_DIG_prop_ran_02 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ 1 + (1|REC_GROUP_REF), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_DIG_prop_ran_02")
B_DIG_prop_ran_03 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ 1 + (1|LITTER_CODE/ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_DIG_prop_ran_03")
B_DIG_prop_ran_04 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_DIG_prop_ran_04")

loo(B_DIG_prop_ran_00, B_DIG_prop_ran_01, B_DIG_prop_ran_02, B_DIG_prop_ran_03, B_DIG_prop_ran_04, moment_match = F)

###### ICC ####
performance::variance_decomposition(B_DIG_prop_ran_00)
performance::r2_bayes(B_DIG_prop_ran_00) 

performance::variance_decomposition(B_DIG_prop_ran_01)
performance::r2_bayes(B_DIG_prop_ran_01)

performance::variance_decomposition(B_DIG_prop_ran_02)
performance::r2_bayes(B_DIG_prop_ran_02) 

performance::variance_decomposition(B_DIG_prop_ran_03)
performance::r2_bayes(B_DIG_prop_ran_03)

performance::variance_decomposition(B_DIG_prop_ran_04)
performance::r2_bayes(B_DIG_prop_ran_04)

rm(B_DIG_prop_ran_00, B_DIG_prop_ran_01, B_DIG_prop_ran_02, B_DIG_prop_ran_03, B_DIG_prop_ran_04)

#### 2) Define all (biologically plausible) models ####
B_DIG_prop_00 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ +1 + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_DIG_prop_00")
# only one predictor
B_DIG_prop_01 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_DIG_prop_01")
# additive predictors - no interactions
B_DIG_prop_02 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_DIG_prop_02")
B_DIG_prop_03 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + SEX + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_DIG_prop_03")
B_DIG_prop_04 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_DIG_prop_04")
B_DIG_prop_04_a <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT + REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                             save_pars = save_pars(all=T), cores=4, file = "B_DIG_prop_04_a")
B_DIG_prop_05 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER+  (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_DIG_prop_05")
# interactions
B_DIG_prop_06 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT *  REC_AGE_D + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, init_r=0.11, file = "B_DIG_prop_06")
B_DIG_prop_07 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT *  REC_AGE_D * SEX + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, init_r=0.11, file = "B_DIG_prop_07")
B_DIG_prop_08 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, init_r=0.11, file = "B_DIG_prop_08")
B_DIG_prop_08_a <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT * REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                             save_pars = save_pars(all=T), cores=4, init_r=0.005, file = "B_DIG_prop_08_a")
B_DIG_prop_09 <- brms::brm(formula =Sum_DIG | trials(Total_calls) ~ TREATMENT * REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|ID), data = DIG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, init_r=0.005, file = "B_DIG_prop_09")

#### CHECK MODELS ####
summary(B_DIG_prop_00) 
#posterior_summary(B_DIG_prop_00)
plot(B_DIG_prop_00) #ok
pp_check(B_DIG_prop_00, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_00, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_01) 
#posterior_summary(B_DIG_prop_01)
plot(B_DIG_prop_01) #ok
pp_check(B_DIG_prop_01, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_01, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_02) 
#posterior_summary(B_DIG_prop_02)
plot(B_DIG_prop_02) #ok
pp_check(B_DIG_prop_02, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_02, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_03) 
#posterior_summary(B_DIG_prop_03)
plot(B_DIG_prop_03) #ok
pp_check(B_DIG_prop_03, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_03, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_04) 
#posterior_summary(B_DIG_prop_04)
plot(B_DIG_prop_04) #ok
pp_check(B_DIG_prop_04, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_04, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_04_a) 
#posterior_summary(B_DIG_prop_04_a)
plot(B_DIG_prop_04_a) #ok
pp_check(B_DIG_prop_04_a, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_04_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_05) 
#posterior_summary(B_DIG_prop_05)
plot(B_DIG_prop_05) #ok
pp_check(B_DIG_prop_05, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_05, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_06) 
#posterior_summary(B_DIG_prop_06)
plot(B_DIG_prop_06) #ok
pp_check(B_DIG_prop_06, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_06, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_07) 
#posterior_summary(B_DIG_prop_07)
plot(B_DIG_prop_07) #ok
pp_check(B_DIG_prop_07, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_07, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_08) 
#posterior_summary(B_DIG_prop_08)
plot(B_DIG_prop_08) #ok
pp_check(B_DIG_prop_08, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_08, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_08_a) 
#posterior_summary(B_DIG_prop_08_a)
plot(B_DIG_prop_08_a) #ok
pp_check(B_DIG_prop_08_a, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_08_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_DIG_prop_09) 
#posterior_summary(B_DIG_prop_09)
plot(B_DIG_prop_09) #ok
pp_check(B_DIG_prop_09, ndraws = 100) #ok
model_loo <- loo(B_DIG_prop_09, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

rm(model_loo, k_rintercept, df)

#### MODEL COMPARISON ####
# leave-one-out-cross-validation (Vehtrai et al, 2017)
loo(B_DIG_prop_00, B_DIG_prop_01, B_DIG_prop_02, B_DIG_prop_03, B_DIG_prop_04, 
    B_DIG_prop_04_a, B_DIG_prop_05, B_DIG_prop_06, B_DIG_prop_07, B_DIG_prop_08, 
    B_DIG_prop_08_a, B_DIG_prop_09)

#### MOMENT MATCHING ####
loo(B_DIG_prop_00, B_DIG_prop_01, B_DIG_prop_02, B_DIG_prop_03, B_DIG_prop_04, 
    B_DIG_prop_04_a, B_DIG_prop_05, B_DIG_prop_06, B_DIG_prop_07, B_DIG_prop_08, 
    B_DIG_prop_08_a, B_DIG_prop_09, moment_match = T)

#### BEST MODEL ####
summary(B_DIG_prop_09)

mcmc_plot(B_DIG_prop_09, type = 'intervals')
mcmc_plot(B_DIG_prop_09, type = 'acf')

pd <- p_direction(B_DIG_prop_09)
plot(pd)

diagnostic_posterior(B_DIG_prop_09)
(cred_int <- ci(B_DIG_prop_09, ci=c(.89, .95), effects=c('fixed')))

# to odds ratio
exp(fixef(B_DIG_prop_09))

# to probability
inv_logit_scaled(fixef(B_DIG_prop_09))

# ROPE: different definition for logistic regression: -0.18, 0.18
rope_rang <- c((sd(DIG_prop_data$Avg_DIG_prop)*0.18)*-1, sd(DIG_prop_data$Avg_DIG_prop)*0.18)

describe_posterior(
  B_DIG_prop_09,
  effects = "fixed", #fixed vs all
  component = "all",
  rope_range = rope_rang,
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

equivalence_test(B_DIG_prop_09, range=rope_rang)

post_best <- describe_posterior(
  B_DIG_prop_09,
  effects = "fixed",
  component = "all",
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)
# convert to odds ratio!
post_best$Median <- exp(post_best$Median)
post_best$CI_high <- exp(post_best$CI_high)
post_best$CI_low <- exp(post_best$CI_low)
post_best <- post_best[-c(4,5)]
post_best

loo_R2(B_DIG_prop_09)
performance::variance_decomposition(B_DIG_prop_09)

################################################################################
## CALCULATE ESTIMATES FOR INTERACTIONS ####
helper_vars <- c(mean(DIG_prop_data$Helper_Pup_REC)-sd(DIG_prop_data$Helper_Pup_REC), 
                 mean(DIG_prop_data$Helper_Pup_REC), 
                 mean(DIG_prop_data$Helper_Pup_REC) + sd(DIG_prop_data$Helper_Pup_REC))

#treatment *age
(treat_age <- emtrends(B_DIG_prop_09, pairwise ~ TREATMENT, var="REC_AGE_D"))

pd(treat_age)

# age * helper
(age_helper <- emtrends(B_DIG_prop_09, specs = c('Helper_Pup_REC'), var = "REC_AGE_D",
                                 at = list(Helper_Pup_REC = helper_vars)))


pd(age_helper)

# treatment * helper * weight
(treat_helper_weight <- emtrends(B_DIG_prop_09, specs = c('TREATMENT', 'Helper_Pup_REC'), var = "WEIGHT_DIFF_PER",
                              at = list(Helper_Pup_REC = helper_vars)))

pd(treat_helper_weight)

### PLOTS: BEST MODEL ####
# plot all predictors
conditional_effects(B_DIG_prop_09)

#TREATMENT*AGE
min(DIG_prop_data$REC_AGE_D)#31
max(DIG_prop_data$REC_AGE_D)#130

age_dist <- B_DIG_prop_09 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_prop_data$TREATMENT),
                                    SEX = levels(DIG_prop_data$SEX),
                                    Helper_Pup_REC = mean(DIG_prop_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_prop_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5),
                                    Total_calls = 1), 
              re_formula = NULL, allow_new_levels=T)

age_dist$Avg_DIG_prop <- age_dist$.epred
# age & treatment
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_DIG_prop, color=TREATMENT, fill=TREATMENT)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = DIG_prop_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Digging call proportion") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  theme_clean()

#AGE*HELPER
helper_vars <- c(mean(DIG_prop_data$Helper_Pup_REC)-sd(DIG_prop_data$Helper_Pup_REC), 
                 mean(DIG_prop_data$Helper_Pup_REC), 
                 mean(DIG_prop_data$Helper_Pup_REC) + sd(DIG_prop_data$Helper_Pup_REC))

helper_dist <- B_DIG_prop_09 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_prop_data$TREATMENT),
                                    SEX = levels(DIG_prop_data$SEX),
                                    Helper_Pup_REC = helper_vars,
                                    WEIGHT_DIFF_PER = mean(DIG_prop_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5),
                                    Total_calls = 1), 
              re_formula = NULL, allow_new_levels=T)

helper_dist$Helper_Pup_REC <- as.factor(helper_dist$Helper_Pup_REC)
helper_dist$Avg_DIG_prop <- helper_dist$.epred

ggplot(helper_dist, aes(x = REC_AGE_D, y = Avg_DIG_prop, color=Helper_Pup_REC, fill=Helper_Pup_REC)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = DIG_prop_data, size = 2, color='black', fill='black',alpha=0.7, shape='+') +   # raw data
  scale_color_okabe_ito(name = "Helper/Pup ratio", labels=c('low', 'medium', 'high'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Helper/Pup ratio", labels=c('low', 'medium', 'high'))+
  labs(x = "Age (days)", y = "Digging call proportion") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  theme_clean()

#TREATMENT*HELPER*WEIGHT
helper_vars <- c(mean(DIG_prop_data$Helper_Pup_REC)-sd(DIG_prop_data$Helper_Pup_REC), 
                 mean(DIG_prop_data$Helper_Pup_REC), 
                 mean(DIG_prop_data$Helper_Pup_REC) + sd(DIG_prop_data$Helper_Pup_REC))

min(DIG_prop_data$WEIGHT_DIFF_PER)#-46
max(DIG_prop_data$WEIGHT_DIFF_PER)#48

# # labels for graphs
treatment <- c(
  "CTRL" = "Control",
  "FLUT" = "Flutamide"
)

treatment_dist <- B_DIG_prop_09 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_prop_data$TREATMENT),
                                    SEX = levels(DIG_prop_data$SEX),
                                    Helper_Pup_REC = helper_vars,
                                    WEIGHT_DIFF_PER = seq(-48, 48, by=2),
                                    REC_AGE_D = mean(DIG_prop_data$REC_AGE_D),
                                    Total_calls = 1), 
              re_formula = NULL, allow_new_levels=T)

treatment_dist$Helper_Pup_REC <- as.factor(treatment_dist$Helper_Pup_REC)
treatment_dist$Avg_DIG_prop <- treatment_dist$.epred

ggplot(treatment_dist, aes(x = WEIGHT_DIFF_PER, y = Avg_DIG_prop, color=Helper_Pup_REC, fill=Helper_Pup_REC)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = DIG_prop_data, size = 2, color='black', fill='black',alpha=0.7, shape='+') +   # raw data
  scale_color_okabe_ito(name = "Helper/Pup ratio", labels=c('low', 'medium', 'high'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Helper/Pup ratio", labels=c('low', 'medium', 'high'))+
  labs(x = "Weight offset (%)", y = "Digging call proportion") +
  facet_wrap(~TREATMENT, labeller = as_labeller(treatment))+
  theme_clean()

########################################################
# AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA)####
# no new data, just get the information of fitted values
post_avg <- posterior_average(B_DIG_prop_09, B_DIG_prop_05, weights = "stacking", 
                              missing = 0)

# posterior_summary(post_avg)
# ci(post_avg, ci=c(.89, .95)) 

# ROPE: different definition for logistic regression: -0.18, 0.18
rope_rang <- c((sd(DIG_prop_data$Avg_DIG_prop)*0.18)*-1, sd(DIG_prop_data$Avg_DIG_prop)*0.18)

describe_posterior(
  post_avg,
  effects = "fixed",
  component = "all",
  rope_range = rope_rang,
  centrality = "all",
  test = c("p_direction", "p_significance", "rope"),
  dispersion = TRUE
)

rm(post_avg)

### cleanup ####
rm(B_DIG_prop_00, B_DIG_prop_01, B_DIG_prop_02, B_DIG_prop_03, B_DIG_prop_04, 
   B_DIG_prop_04_a, B_DIG_prop_05, B_DIG_prop_06, B_DIG_prop_07, B_DIG_prop_08, 
   B_DIG_prop_08_a, B_DIG_prop_09, DIG_prop_data, DIG_data, full_data)