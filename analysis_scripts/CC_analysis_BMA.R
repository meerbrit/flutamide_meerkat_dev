#### Analysis script for Flutamide study: CC calls ###################
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
library(ggokabeito)   #  color palette
library(emmeans) # emtrends
library(ggnewscale) # change colour within graph
library(marginaleffects) # predictions

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
CC_data <- full_data[, c(34:48, 79, 105:113, 118, 120,  122, 125:126, 128, 147, 162, 166:169)]
# remove NAs
CC_data <- na.omit(CC_data) 
# how many observations of each sex in FLUT & CTRL?
table(CC_data$TREATMENT, CC_data$SEX) 

#create the different datasets for each response variable
CC_len_data <- subset(CC_data, CC_avg_Len >0 & CC_num >= MIN_BOUTS) # call length
CC_len_data$Avg_CC_mil <- (CC_len_data$CC_avg_Len)*1000 #convert avg_CC_LEN to milliseconds
CC_int_data <- subset(CC_data, CC_avg_Int >0 & CC_num >= MIN_BOUTS) # call interval length
CC_int_data$Avg_CC_mil <- (CC_int_data$CC_avg_Int)*1000 #convert avg_CC_INT to milliseconds
CC_rate_data <- subset(CC_data, CC_rate >0 & CC_num >= MIN_BOUTS) # call rate
CC_prop_data <- subset(CC_data, nBouts >= MIN_BOUTS) #call proportion
CC_prop_data$Total_calls <- CC_prop_data$Sum_BEG+
  CC_prop_data$Sum_DIG+CC_prop_data$Sum_CC+CC_prop_data$Sum_OTHER #save the number of total calls for calc later

# switch working directory to models folder
setwd("G:/My Drive/Uni/UZH/Projects/KMP_PUPS/Flutamide/scripts/CC_models/")

# ggplot theme 
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
corvif(sapply(CC_data[, 
                       c("REC_AGE_D", "SEX", "TREATMENT",
                         "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                         "CONDITION", "CONDITION_PER")], as.numeric))

corvif(sapply(CC_data[, 
                       c("REC_AGE_D", "SEX", "TREATMENT",
                         "Helper_Pup_REC", "WEIGHT_DIFF_PER")], as.numeric))

################################################################################
######################## Call length ###########################################
################################################################################

#### EXPLORE DATA ####
# how many observations of each sex in FLUT & CTRL?
table(CC_len_data$TREATMENT, CC_len_data$SEX) 
# unique individuals
sample_df <- CC_len_data[!duplicated(CC_len_data$ID),] 
table(sample_df$SEX)
table(sample_df$TREATMENT)

# test for normality of data
ggqqplot(CC_len_data$Avg_CC_mil) 
shapiro.test(CC_len_data$Avg_CC_mil) 

# use a histogram to check data distribution
ggplot(CC_len_data, aes(Avg_CC_mil)) +
  geom_histogram(fill= "white", col= "black", binwidth = 10) +
  facet_wrap(~TREATMENT) +
  theme_bw()
# by sex and treatment
ggplot(CC_len_data, aes(Avg_CC_mil)) +
  geom_histogram(fill= "white", col= "black", binwidth = 10) +
  facet_wrap(TREATMENT~SEX) +
  theme_bw()

# check for correlation of vars
cor(sapply(CC_len_data[, c("Avg_CC_mil", "REC_AGE_D", "SEX", "TREATMENT",
                            "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                            "CONDITION", "CONDITION_PER")], as.numeric))

# plot a simple overview for treatment/age/sex
ggplot(CC_len_data, aes(REC_AGE_D, Avg_CC_mil)) +
  geom_point() +
  facet_grid(TREATMENT~SEX) +
  theme_bw()

#check if there is an obvious treatment effect
plot(CC_len_data$Avg_CC_mil~CC_len_data$TREATMENT)
fit <- lm(Avg_CC_mil~TREATMENT * REC_AGE_D * SEX, CC_len_data) 
coef(fit)
S(fit)
rm(fit)

# overview plot 
plot_CC_len <- ggplot(CC_len_data, aes( REC_AGE_D, Avg_CC_mil)) +
  geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
  geom_smooth(method = "lm", col= 1) +
  labs(y= "CC call length") +
  facet_wrap(~ TREATMENT:SEX) +
  theme_bw()
annotate_figure(plot_CC_len, 
                top = text_grob("Avg_CC_mil", 
                                color = "red", face = "bold", size = 14))
rm(plot_CC_len)

####OUTLIER CHECK ####
ggplot(CC_len_data) +
  aes(x = "", y = Avg_CC_mil) +
  geom_boxplot(fill = "#0c4c8a") +
  scale_color_okabe_ito()+
  theme_clean()+
  facet_grid(~SEX:TREATMENT)
# get the outliers (including indices in df)
out <- boxplot.stats(CC_len_data$Avg_CC_mil)$out
(out_ind <- which(CC_len_data$Avg_CC_mil %in% c(out)))#40

# which datapoints seem to be outliers?
CC_len_data[out_ind, ]
CC_len_data[out_ind, ]$Avg_CC_mil

# Rosner Test for outliers to double check
# k= number of suspected outliers
test <- rosnerTest(CC_len_data$Avg_CC_mil,
                   k = length(out))
test

#### BAYES MODELS ##############################################################
#### 1) Test random effects and nested level
B_CC_len_ran_00 <- brms::brm(formula =Avg_CC_mil ~ 1 + (1|ID), data = CC_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_len_ran_00")
B_CC_len_ran_01 <- brms::brm(formula =Avg_CC_mil ~ 1 + (1|LITTER_CODE), data = CC_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_len_ran_01")
B_CC_len_ran_02 <- brms::brm(formula =Avg_CC_mil ~ 1 + (1|REC_GROUP_REF), data = CC_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_len_ran_02")
B_CC_len_ran_03 <- brms::brm(formula =Avg_CC_mil ~ 1 + (1|LITTER_CODE/ID), data = CC_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_len_ran_03")
B_CC_len_ran_04 <- brms::brm(formula =Avg_CC_mil ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = CC_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_len_ran_04")

loo(B_CC_len_ran_00, B_CC_len_ran_01, B_CC_len_ran_02, B_CC_len_ran_03, B_CC_len_ran_04, moment_match = T)

###### ICC ####
performance::variance_decomposition(B_CC_len_ran_00)
performance::r2_bayes(B_CC_len_ran_00) 

performance::variance_decomposition(B_CC_len_ran_01)
performance::r2_bayes(B_CC_len_ran_01)

performance::variance_decomposition(B_CC_len_ran_02)
performance::r2_bayes(B_CC_len_ran_02)

performance::variance_decomposition(B_CC_len_ran_03)
performance::r2_bayes(B_CC_len_ran_03)

performance::variance_decomposition(B_CC_len_ran_04)
performance::r2_bayes(B_CC_len_ran_04)

rm(B_CC_len_ran_00, B_CC_len_ran_01, B_CC_len_ran_02, B_CC_len_ran_03, B_CC_len_ran_04)

#### 2) Define all (biologically plausible) models ####
B_CC_len_00 <- brms::brm(formula =Avg_CC_mil ~ +1 + (1|ID), data = CC_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_len_00")
# only one predictor
B_CC_len_01 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT + (1|ID), data = CC_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_len_01")
# additive predictors - no interactions
B_CC_len_02 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT +  REC_AGE_D + (1|ID), data = CC_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_len_02")
B_CC_len_03 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT +  REC_AGE_D + SEX + (1|ID), data = CC_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_len_03")
B_CC_len_04 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|ID), data = CC_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_len_04")
B_CC_len_04_a <- brms::brm(formula =Avg_CC_mil ~ TREATMENT + REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|ID), data = CC_len_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                            save_pars = save_pars(all=T), cores=4, file = "B_CC_len_04_a")
B_CC_len_05 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER+  (1|ID), data = CC_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_len_05")
# interactions
B_CC_len_06 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT *  REC_AGE_D + (1|ID), data = CC_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_len_06")
B_CC_len_07 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT *  REC_AGE_D * SEX + (1|ID), data = CC_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_len_07")
B_CC_len_08 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|ID), data = CC_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_len_08")
B_CC_len_08_a <- brms::brm(formula =Avg_CC_mil ~ TREATMENT * REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|ID), data = CC_len_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                            save_pars = save_pars(all=T), cores=4, file = "B_CC_len_08_a")
B_CC_len_09 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT *REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|ID), data = CC_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_len_09")

#### CHECK MODELS ####
summary(B_CC_len_00) 
#posterior_summary(B_CC_len_00)
plot(B_CC_len_00) 
pp_check(B_CC_len_00, ndraws = 100) 
# Pareto-k diagnostic can be useful to identify problematic point(s)
model_loo <- loo(B_CC_len_00, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_01) 
#posterior_summary(B_CC_len_01)
plot(B_CC_len_01) 
pp_check(B_CC_len_01, ndraws = 100) 
model_loo <- loo(B_CC_len_01, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_02) 
#posterior_summary(B_CC_len_02)
plot(B_CC_len_02) 
pp_check(B_CC_len_02, ndraws = 100) 
model_loo <- loo(B_CC_len_02, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_03) 
#posterior_summary(B_CC_len_03)
plot(B_CC_len_03) 
pp_check(B_CC_len_03, ndraws = 100) 
model_loo <- loo(B_CC_len_03, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_04) 
#posterior_summary(B_CC_len_04)
plot(B_CC_len_04) 
pp_check(B_CC_len_04, ndraws = 100) 
model_loo <- loo(B_CC_len_04, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_04_a) 
#posterior_summary(B_CC_len_04_a)
plot(B_CC_len_04_a) 
pp_check(B_CC_len_04_a, ndraws = 100) 
model_loo <- loo(B_CC_len_04_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_05)
#posterior_summary(B_CC_len_05)
plot(B_CC_len_05) 
pp_check(B_CC_len_05, ndraws = 100) 
model_loo <- loo(B_CC_len_05, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_06) 
#posterior_summary(B_CC_len_06)
plot(B_CC_len_06) 
pp_check(B_CC_len_06, ndraws = 100) 
model_loo <- loo(B_CC_len_06, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_07) 
#posterior_summary(B_CC_len_07)
plot(B_CC_len_07) 
pp_check(B_CC_len_07, ndraws = 100) 
model_loo <- loo(B_CC_len_07, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_08) 
#posterior_summary(B_CC_len_08)
plot(B_CC_len_08) 
pp_check(B_CC_len_08, ndraws = 100) 
model_loo <- loo(B_CC_len_08, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_08_a) 
#posterior_summary(B_CC_len_08_a)
plot(B_CC_len_08_a) 
pp_check(B_CC_len_08_a, ndraws = 100) 
model_loo <- loo(B_CC_len_08_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_len_09) 
#posterior_summary(B_CC_len_09)
plot(B_CC_len_09) 
pp_check(B_CC_len_09, ndraws = 100) 
model_loo <- loo(B_CC_len_09, save_psis = TRUE, cores=2)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

rm(df, model_loo, k_rintercept)

#### MODEL COMPARISON ####
# leave-one-out-cross-validation (Vehtrai et al, 2017)
loo(B_CC_len_00, B_CC_len_01, B_CC_len_02, B_CC_len_03, B_CC_len_04, 
    B_CC_len_04_a, B_CC_len_05, B_CC_len_06, B_CC_len_07, B_CC_len_08, 
    B_CC_len_08_a, B_CC_len_09)

#### MOMENT MATCHING ####
loo(B_CC_len_00, B_CC_len_01, B_CC_len_02, B_CC_len_03, B_CC_len_04, 
    B_CC_len_04_a, B_CC_len_05, B_CC_len_06, B_CC_len_07, B_CC_len_08, 
    B_CC_len_08_a, B_CC_len_09, moment_match = T)

#### BEST MODEL ####
summary(B_CC_len_08)

mcmc_plot(B_CC_len_08, type = 'intervals')
mcmc_plot(B_CC_len_08, type = 'acf')
#mcmc_areas(B_CC_len_08)

pd <- p_direction(B_CC_len_08)
plot(pd)

diagnostic_posterior(B_CC_len_08)
(cred_int <- ci(B_CC_len_08, ci=c(.89, .95), effects=c('fixed')))

describe_posterior(
  B_CC_len_08,
  effects = "fixed", # fixed vs all
  component = "all",
  rope_range = rope_range(B_CC_len_08),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE)

equivalence_test(B_CC_len_08)

loo_R2(B_CC_len_08) 
performance::variance_decomposition(B_CC_len_08)

# plot all predictors
cond_eff <- conditional_effects(B_CC_len_08, re_formula = NULL, robust=TRUE) # include all random effects (vs NA), robust uses median
plot(cond_eff)

# AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA)####
# no new data, just get the information of fitted values
post_avg <- posterior_average(B_CC_len_08, B_CC_len_08_a, weights = "stacking", 
                              missing = 0)

#posterior_summary(post_avg)
# ci(post_avg, ci=c(.89, .95)) 

describe_posterior(
  post_avg,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_CC_len_08),
  test = c("p_direction", "p_significance", "rope"),
  centrality = 'all',
  dispersion = TRUE
)

rm(post_avg)

### cleanup ####
rm(B_CC_len_00, B_CC_len_01, B_CC_len_02, B_CC_len_03, B_CC_len_04, 
   B_CC_len_04_a, B_CC_len_05, B_CC_len_06, B_CC_len_07, B_CC_len_08, 
   B_CC_len_08_a, B_CC_len_09, CC_len_data)

################################################################################
######################## Call interval #########################################
################################################################################

# 1) explore data ####
# how many observations of each sex in FLUT & CTRL?
table(CC_int_data$TREATMENT, CC_int_data$SEX) 

# test for normality of data
ggqqplot(CC_int_data$Avg_CC_mil)
shapiro.test(CC_int_data$Avg_CC_mil)

# use a histogram to check data distribution
ggplot(CC_int_data, aes(Avg_CC_mil)) +
  geom_histogram(fill= "white", col= "black", binwidth = 300) +
  facet_wrap(~TREATMENT) +
  theme_bw()
# by sex and treatment
ggplot(CC_int_data, aes(Avg_CC_mil)) +
  geom_histogram(fill= "white", col= "black", binwidth = 300) +
  facet_wrap(TREATMENT~SEX) +
  theme_bw()

# check for correlation of vars
cor(sapply(CC_int_data[, c("Avg_CC_mil", "REC_AGE_D", "SEX", "TREATMENT",
                            "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                            "CONDITION", "CONDITION_PER")], as.numeric))

# plot a simple overview for treatment/age/sex
ggplot(CC_int_data, aes(REC_AGE_D, Avg_CC_mil)) +
  geom_point() +
  facet_grid(TREATMENT~SEX) +
  theme_bw()

#check if there is an obvious treatment effect
plot(CC_int_data$Avg_CC_mil~CC_int_data$TREATMENT)
fit <- lm(Avg_CC_mil~TREATMENT * REC_AGE_D * SEX, CC_int_data) 
coef(fit)
S(fit)
rm(fit)

# overview plot ########
plot_CC_int <- ggplot(CC_int_data, aes( REC_AGE_D, Avg_CC_mil)) +
  geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
  geom_smooth(method = "lm", col= 1) +
  labs(y= "CC call interval length") +
  facet_wrap(~ TREATMENT:SEX) +
  theme_bw()
annotate_figure(plot_CC_int, 
                top = text_grob("Avg_CC_mil", 
                                color = "red", face = "bold", size = 14))

rm(plot_CC_int)

####OUTLIER CHECK ####
ggplot(CC_int_data) +
  aes(x = "", y = Avg_CC_mil) +
  geom_boxplot(fill = "#0c4c8a") +
  scale_color_okabe_ito()+
  theme_clean()+
  facet_grid(~SEX:TREATMENT)
# get the outliers (including indices in df)
out <- boxplot.stats(CC_int_data$Avg_CC_mil)$out
(out_ind <- which(CC_int_data$Avg_CC_mil %in% c(out)))#4, 72

# which datapoints seem to be outliers?
CC_int_data[out_ind, ]
CC_int_data[out_ind, ]$Avg_CC_mil

# Rosner Test for outliers to double check
# k= number of suspected outliers
test <- rosnerTest(CC_int_data$Avg_CC_mil,
                   k = length(out))
test 

#### BAYES MODELS ##############################################################
#### 1) Test random effects and nested level
B_CC_int_ran_00 <- brms::brm(formula =Avg_CC_mil ~ 1 + (1|ID), data = CC_int_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_int_ran_00")
B_CC_int_ran_01 <- brms::brm(formula =Avg_CC_mil ~ 1 + (1|LITTER_CODE), data = CC_int_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_int_ran_01")
B_CC_int_ran_02 <- brms::brm(formula =Avg_CC_mil ~ 1 + (1|REC_GROUP_REF), data = CC_int_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_int_ran_02")
B_CC_int_ran_03 <- brms::brm(formula =Avg_CC_mil ~ 1 + (1|LITTER_CODE/ID), data = CC_int_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_int_ran_03")
B_CC_int_ran_04 <- brms::brm(formula =Avg_CC_mil ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = CC_int_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_int_ran_04")

loo(B_CC_int_ran_00, B_CC_int_ran_01, B_CC_int_ran_02, B_CC_int_ran_03, B_CC_int_ran_04, moment_match = F)

###### ICC ####
performance::variance_decomposition(B_CC_int_ran_00)
performance::r2_bayes(B_CC_int_ran_00) 

performance::variance_decomposition(B_CC_int_ran_01)
performance::r2_bayes(B_CC_int_ran_01)

performance::variance_decomposition(B_CC_int_ran_02)
performance::r2_bayes(B_CC_int_ran_02)

performance::variance_decomposition(B_CC_int_ran_03)
performance::r2_bayes(B_CC_int_ran_03)

performance::variance_decomposition(B_CC_int_ran_04)
performance::r2_bayes(B_CC_int_ran_04)

rm(B_CC_int_ran_00, B_CC_int_ran_01, B_CC_int_ran_02,
   B_CC_int_ran_03, B_CC_int_ran_04)

#### 2) Define all (biologically plausible) models ####
B_CC_int_00 <- brms::brm(formula =Avg_CC_mil ~ +1 + (1|ID), data = CC_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_int_00")
# only one predictor
B_CC_int_01 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT + (1|ID), data = CC_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_int_01")
# additive predictors - no interactions
B_CC_int_02 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT +  REC_AGE_D + (1|ID), data = CC_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_int_02")
B_CC_int_03 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT +  REC_AGE_D + SEX + (1|ID), data = CC_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_int_03")
B_CC_int_04 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|ID), data = CC_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_int_04")
B_CC_int_04_a <- brms::brm(formula =Avg_CC_mil ~ TREATMENT + REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|ID), data = CC_int_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, file = "B_CC_int_04_a")
B_CC_int_05 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER+  (1|ID), data = CC_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_int_05")
# interactions
B_CC_int_06 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT *  REC_AGE_D + (1|ID), data = CC_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_int_06")
B_CC_int_07 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT *  REC_AGE_D * SEX + (1|ID), data = CC_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_int_07")
B_CC_int_08 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|ID), data = CC_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_int_08")
B_CC_int_08_a <- brms::brm(formula =Avg_CC_mil ~ TREATMENT * REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|ID), data = CC_int_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                            save_pars = save_pars(all=T), cores=4, file = "B_CC_int_08_a")
B_CC_int_09 <- brms::brm(formula =Avg_CC_mil ~ TREATMENT *REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|ID), data = CC_int_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_int_09")

#### CHECK MODELS ####
summary(B_CC_int_00) 
#posterior_summary(B_CC_int_00)
plot(B_CC_int_00) 
pp_check(B_CC_int_00, ndraws = 100) 
model_loo <- loo(B_CC_int_00, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_01) 
#posterior_summary(B_CC_int_01)
plot(B_CC_int_01) 
pp_check(B_CC_int_01, ndraws = 100) #
model_loo <- loo(B_CC_int_01, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_02) 
#posterior_summary(B_CC_int_02)
plot(B_CC_int_02) 
pp_check(B_CC_int_02, ndraws = 100) #
model_loo <- loo(B_CC_int_02, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_03) 
#posterior_summary(B_CC_int_03)
plot(B_CC_int_03) 
pp_check(B_CC_int_03, ndraws = 100) #
model_loo <- loo(B_CC_int_03, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_04) 
#posterior_summary(B_CC_int_04)
plot(B_CC_int_04) 
pp_check(B_CC_int_04, ndraws = 100) #
model_loo <- loo(B_CC_int_04, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_04_a) 
#posterior_summary(B_CC_int_04_a)
plot(B_CC_int_04_a) 
pp_check(B_CC_int_04_a, ndraws = 100) #
model_loo <- loo(B_CC_int_04_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_05) 
#posterior_summary(B_CC_int_05)
plot(B_CC_int_05) 
pp_check(B_CC_int_05, ndraws = 100) #
model_loo <- loo(B_CC_int_05, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_06) 
#posterior_summary(B_CC_int_06)
plot(B_CC_int_06) 
pp_check(B_CC_int_06, ndraws = 100) #
model_loo <- loo(B_CC_int_06, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_07) 
#posterior_summary(B_CC_int_07)
plot(B_CC_int_07) 
pp_check(B_CC_int_07, ndraws = 100) #
model_loo <- loo(B_CC_int_07, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_08) 
#posterior_summary(B_CC_int_08)
plot(B_CC_int_08) 
pp_check(B_CC_int_08, ndraws = 100) #
model_loo <- loo(B_CC_int_08, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_08_a) 
#posterior_summary(B_CC_int_08_a)
plot(B_CC_int_08_a) 
pp_check(B_CC_int_08_a, ndraws = 100) #
model_loo <- loo(B_CC_int_08_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_int_09) 
#posterior_summary(B_CC_int_09)
plot(B_CC_int_09) 
pp_check(B_CC_int_09, ndraws = 100) #
model_loo <- loo(B_CC_int_09, save_psis = TRUE, cores=4)
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
loo(B_CC_int_00, B_CC_int_01, B_CC_int_02, B_CC_int_03, B_CC_int_04, 
    B_CC_int_04_a, B_CC_int_05, B_CC_int_06, B_CC_int_07, B_CC_int_08, 
    B_CC_int_08_a, B_CC_int_09)

#### MOMENT MATCHING ####
loo(B_CC_int_00, B_CC_int_01, B_CC_int_02, B_CC_int_03, B_CC_int_04, 
    B_CC_int_04_a, B_CC_int_05, B_CC_int_06, B_CC_int_07, B_CC_int_08, 
    B_CC_int_08_a, B_CC_int_09, moment_match = T)

#### BEST MODEL ####
summary(B_CC_int_02)

mcmc_plot(B_CC_int_02, type = 'intervals')
mcmc_plot(B_CC_int_02, type = 'acf')

pd <- p_direction(B_CC_int_02)
plot(pd)

diagnostic_posterior(B_CC_int_02)
(cred_int <- ci(B_CC_int_02, ci=c(.89, .95), effects=c('fixed')))

describe_posterior(
  B_CC_int_02,
  effects = "all", # all vs fixed
  component = "all",
  rope_range= rope_range(B_CC_int_02),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

equivalence_test(B_CC_int_02)

loo_R2(B_CC_int_02)
performance::variance_decomposition(B_CC_int_02)

####PLOTS: BEST MODEL ####
# plot all predictors
cond_eff <- conditional_effects(B_CC_int_02, re_formula = NULL, robust=TRUE)
plot(cond_eff)

# AGE
min(CC_int_data$REC_AGE_D)#31
max(CC_int_data$REC_AGE_D)#130

age_dist <- B_CC_int_02 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(CC_int_data$TREATMENT),
                                    SEX = levels(CC_int_data$SEX),
                                    Helper_Pup_REC = mean(CC_int_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(CC_int_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5)), 
              re_formula = NULL, allow_new_levels=T)

age_dist$Avg_CC_mil <- age_dist$.epred
#age only
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_CC_mil)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = CC_int_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Contact call interval length (ms)") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  guides(fill='none')+
  theme_clean()
# age & treatment
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_CC_mil, color=TREATMENT, fill=TREATMENT)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = CC_int_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Contact call interval length (ms)") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  theme_clean()

# AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA) ####
# no new data, just get the information of fitted values
post_avg <- posterior_average(B_CC_int_02, B_CC_int_03, B_CC_int_06,
                              B_CC_int_07, B_CC_int_05, B_CC_int_04, 
                              weights = "stacking", 
                              missing = 0)
posterior_summary(post_avg)

#ci(post_avg, ci=c(.89, .95)) 

describe_posterior(
  post_avg,
  effects = "fixed",
  component = "all",
  rope_range= rope_range(B_CC_int_02),
  test = c("p_direction", "p_significance", "rope"),
  #centrality = "all",
  dispersion = TRUE
)

rm(post_avg)

### cleanup #####
rm(B_CC_int_00, B_CC_int_01, B_CC_int_02, B_CC_int_03, B_CC_int_04, 
   B_CC_int_04_a, B_CC_int_05, B_CC_int_06, B_CC_int_07, B_CC_int_08, 
   B_CC_int_08_a, B_CC_int_09, CC_int_data)

################################################################################
######################## Call rate #############################################
################################################################################
# 1) explore data
# how many observations of each sex in FLUT & CTRL?
table(CC_rate_data$TREATMENT, CC_rate_data$SEX) #

# test for normality of data
ggqqplot(CC_rate_data$CC_rate)
shapiro.test(CC_rate_data$CC_rate) 

# use a histogram to check data distribution
ggplot(CC_rate_data, aes(CC_rate)) +
  geom_histogram(fill= "white", col= "black", binwidth = 0.1) +
  facet_wrap(~TREATMENT) +
  theme_bw()
# by sex and treatment
ggplot(CC_rate_data, aes(CC_rate)) +
  geom_histogram(fill= "white", col= "black", binwidth = 0.1) +
  facet_wrap(TREATMENT~SEX) +
  theme_bw()

# check for correlation of vars
cor(sapply(CC_rate_data[, c("CC_rate", "REC_AGE_D", "SEX", "TREATMENT",
                             "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                             "CONDITION", "CONDITION_PER")], as.numeric))

# plot a simple overview for treatment/age/sex
ggplot(CC_rate_data, aes(REC_AGE_D, CC_rate)) +
  geom_point() +
  facet_grid(TREATMENT~SEX) +
  theme_bw()

#check if there is an obvious treatment effect
plot(CC_rate_data$CC_rate~CC_rate_data$TREATMENT)
fit <- lm(CC_rate~TREATMENT * REC_AGE_D * SEX, CC_rate_data) 
coef(fit)
S(fit)
rm(fit)

# overview plot ########
plot_CC_rate <- ggplot(CC_rate_data, aes( REC_AGE_D, CC_rate)) +
  geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
  geom_smooth(method = "lm", col= 1) +
  labs(y= "CC call rate") +
  facet_wrap(~ TREATMENT:SEX) +
  theme_bw()
annotate_figure(plot_CC_rate, 
                top = text_grob("CC_rate", 
                                color = "red", face = "bold", size = 14))
rm(plot_CC_rate)

####OUTLIER CHECK ####
ggplot(CC_rate_data) +
  aes(x = "", y = CC_rate) +
  geom_boxplot(fill = "#0c4c8a") +
  scale_color_okabe_ito()+
  theme_clean()+
  facet_grid(~SEX:TREATMENT)
# get the outliers (including indices in df)
out <- boxplot.stats(CC_rate_data$CC_rate)$out
(out_ind <- which(CC_rate_data$CC_rate %in% c(out)))# no outlier

# which datapoints seem to be outliers?
CC_rate_data[out_ind, ]
CC_rate_data[out_ind, ]$CC_rate

# Rosner Test for outliers to double check
# k= number of suspected outliers
test <- rosnerTest(CC_rate_data$CC_rate,
                   k = length(out))
test  

#### BAYES MODELS ##############################################################
#### 1) Test random effects and nested level
B_CC_rat_ran_00 <- brms::brm(formula =CC_rate ~ 1 + (1|ID), data = CC_rate_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_rat_ran_00")
B_CC_rat_ran_01 <- brms::brm(formula =CC_rate ~ 1 + (1|LITTER_CODE), data = CC_rate_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_rat_ran_01")
B_CC_rat_ran_02 <- brms::brm(formula =CC_rate ~ 1 + (1|REC_GROUP_REF), data = CC_rate_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_rat_ran_02")
B_CC_rat_ran_03 <- brms::brm(formula =CC_rate ~ 1 + (1|LITTER_CODE/ID), data = CC_rate_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_rat_ran_03")
B_CC_rat_ran_04 <- brms::brm(formula =CC_rate ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = CC_rate_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file = "B_CC_rat_ran_04")

loo(B_CC_rat_ran_00, B_CC_rat_ran_01, B_CC_rat_ran_02, B_CC_rat_ran_03, B_CC_rat_ran_04, moment_match = F)

###### ICC ####
performance::variance_decomposition(B_CC_rat_ran_00)
performance::r2_bayes(B_CC_rat_ran_00) # marginal = fixed effects only

performance::variance_decomposition(B_CC_rat_ran_01)
performance::r2_bayes(B_CC_rat_ran_01)

performance::variance_decomposition(B_CC_rat_ran_02)
performance::r2_bayes(B_CC_rat_ran_02) 

performance::variance_decomposition(B_CC_rat_ran_03)
performance::r2_bayes(B_CC_rat_ran_03)

performance::variance_decomposition(B_CC_rat_ran_04)
performance::r2_bayes(B_CC_rat_ran_04)

rm(B_CC_rat_ran_00, B_CC_rat_ran_01, B_CC_rat_ran_02, B_CC_rat_ran_03, B_CC_rat_ran_04)

#### 2) Define all (biologically plausible) models ####
B_CC_rat_00 <- brms::brm(formula = CC_rate~ +1 + (1|ID), data = CC_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_00")
# only one predictor
B_CC_rat_01 <- brms::brm(formula =CC_rate ~ TREATMENT + (1|ID), data = CC_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_01")
# additive predictors - no interactions
B_CC_rat_02 <- brms::brm(formula =CC_rate ~ TREATMENT +  REC_AGE_D + (1|ID), data = CC_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_02")
B_CC_rat_03 <- brms::brm(formula =CC_rate ~ TREATMENT +  REC_AGE_D + SEX + (1|ID), data = CC_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_03")
B_CC_rat_04 <- brms::brm(formula =CC_rate ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|ID), data = CC_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_04")
B_CC_rat_04_a <- brms::brm(formula =CC_rate ~ TREATMENT + REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|ID), data = CC_rate_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_04_a")
B_CC_rat_05 <- brms::brm(formula =CC_rate ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER+  (1|ID), data = CC_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_05")
# interactions
B_CC_rat_06 <- brms::brm(formula =CC_rate ~ TREATMENT *  REC_AGE_D + (1|ID), data = CC_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_06")
B_CC_rat_07 <- brms::brm(formula =CC_rate ~ TREATMENT *  REC_AGE_D * SEX + (1|ID), data = CC_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_07")
B_CC_rat_08 <- brms::brm(formula =CC_rate ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|ID), data = CC_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_08")
B_CC_rat_08_a <- brms::brm(formula =CC_rate ~ TREATMENT * REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|ID), data = CC_rate_data, family = student(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                            save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_08_a")
B_CC_rat_09 <- brms::brm(formula =CC_rate ~ TREATMENT * REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|ID), data = CC_rate_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file = "B_CC_rat_09")

#### CHECK MODELS ####
summary(B_CC_rat_00) 
#posterior_summary(B_CC_rat_00)
plot(B_CC_rat_00) 
pp_check(B_CC_rat_00, ndraws = 100) 
model_loo <- loo(B_CC_rat_00, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_01) 
#posterior_summary(B_CC_rat_01)
plot(B_CC_rat_01) 
pp_check(B_CC_rat_01, ndraws = 100) 
model_loo <- loo(B_CC_rat_01, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_02) 
#posterior_summary(B_CC_rat_02)
plot(B_CC_rat_02) 
pp_check(B_CC_rat_02, ndraws = 100) 
model_loo <- loo(B_CC_rat_02, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_03) 
#posterior_summary(B_CC_rat_03)
plot(B_CC_rat_03) 
pp_check(B_CC_rat_03, ndraws = 100) 
model_loo <- loo(B_CC_rat_03, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_04) 
#posterior_summary(B_CC_rat_04)
plot(B_CC_rat_04) 
pp_check(B_CC_rat_04, ndraws = 100) 
model_loo <- loo(B_CC_rat_04, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_04_a) 
#posterior_summary(B_CC_rat_04_a)
plot(B_CC_rat_04_a) 
pp_check(B_CC_rat_04_a, ndraws = 100) 
model_loo <- loo(B_CC_rat_04_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_05) 
#posterior_summary(B_CC_rat_05)
plot(B_CC_rat_05) 
pp_check(B_CC_rat_05, ndraws = 100) 
model_loo <- loo(B_CC_rat_05, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_06) 
#posterior_summary(B_CC_rat_06)
plot(B_CC_rat_06) 
pp_check(B_CC_rat_06, ndraws = 100) 
model_loo <- loo(B_CC_rat_06, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_07) 
#posterior_summary(B_CC_rat_07)
plot(B_CC_rat_07) 
pp_check(B_CC_rat_07, ndraws = 100) 
model_loo <- loo(B_CC_rat_07, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_08) 
#posterior_summary(B_CC_rat_08)
plot(B_CC_rat_08) 
pp_check(B_CC_rat_08, ndraws = 100) 
model_loo <- loo(B_CC_rat_08, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_08_a) 
#posterior_summary(B_CC_rat_08_a)
plot(B_CC_rat_08_a) 
pp_check(B_CC_rat_08_a, ndraws = 100) 
model_loo <- loo(B_CC_rat_08_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_rat_09) 
#posterior_summary(B_CC_rat_09)
plot(B_CC_rat_09) 
pp_check(B_CC_rat_09, ndraws = 100) 
model_loo <- loo(B_CC_rat_09, save_psis = TRUE, cores=4)
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
loo(B_CC_rat_00, B_CC_rat_01, B_CC_rat_02, B_CC_rat_03, B_CC_rat_04, 
    B_CC_rat_04_a, B_CC_rat_05, B_CC_rat_06, B_CC_rat_07, B_CC_rat_08, 
    B_CC_rat_08_a, B_CC_rat_09)

#### MOMENT MATCHING ####
loo(B_CC_rat_00, B_CC_rat_01, B_CC_rat_02, B_CC_rat_03, B_CC_rat_04, 
    B_CC_rat_04_a, B_CC_rat_05, B_CC_rat_06, B_CC_rat_07, B_CC_rat_08, 
    B_CC_rat_08_a, B_CC_rat_09, moment_match = T)

#### BEST MODEL ####
summary(B_CC_rat_02)
mcmc_plot(B_CC_rat_02, type = 'intervals')
mcmc_plot(B_CC_rat_02, type = 'acf')

pd <- p_direction(B_CC_rat_02)
plot(pd)

diagnostic_posterior(B_CC_rat_02)
(cred_int <- ci(B_CC_rat_02, ci=c(.89, .95), effects=c('fixed')))

describe_posterior(
  B_CC_rat_02,
  effects = "fixed", # fixed vs all
  component = "all",
  rope_range = rope_range(B_CC_rat_02),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

equivalence_test(B_CC_rat_02)

loo_R2(B_CC_rat_02) 
performance::variance_decomposition(B_CC_rat_02)

####PLOTS: BEST MODEL ####
# plot all predictors
cond_eff <- conditional_effects(B_CC_rat_02, re_formula = NULL, robust=TRUE)
plot(cond_eff)

# AGE
min(CC_rate_data$REC_AGE_D)#31
max(CC_rate_data$REC_AGE_D)#130

age_dist <- B_CC_rat_02 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(CC_rate_data$TREATMENT),
                                    SEX = levels(CC_rate_data$SEX),
                                    Helper_Pup_REC = mean(CC_rate_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(CC_rate_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5)), 
              re_formula = NULL, allow_new_levels=T)

age_dist$CC_rate <- age_dist$.epred
#age only
ggplot(age_dist, aes(x = REC_AGE_D, y = CC_rate)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = CC_rate_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Contact call rate") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  guides(fill='none')+
  theme_clean()
# age & treatment
ggplot(age_dist, aes(x = REC_AGE_D, y = CC_rate, color=TREATMENT, fill=TREATMENT)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = CC_rate_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Contact call rate") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  theme_clean() 

# AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA) ####
# no new data, just get the information of fitted values
post_avg <- posterior_average(B_CC_rat_02, B_CC_rat_07, B_CC_rat_03, B_CC_rat_08,
                              B_CC_rat_04_a, B_CC_rat_04, B_CC_rat_06, B_CC_rat_05,
                              B_CC_rat_09,
                              weights = "stacking", 
                              missing = 0)

posterior_summary(post_avg)
# ci(post_avg, ci=c(.89, .95)) 

describe_posterior(
  post_avg,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_CC_rat_02),
  test = c("p_direction", "p_significance", "rope"),
#  centrality = 'all',
  dispersion = TRUE
)

rm(post_avg)

### cleanup ####
rm(B_CC_rat_00, B_CC_rat_01, B_CC_rat_02, B_CC_rat_03, B_CC_rat_04, 
   B_CC_rat_04_a, B_CC_rat_05, B_CC_rat_06, B_CC_rat_07, B_CC_rat_08, 
   B_CC_rat_08_a, B_CC_rat_09, CC_rate_data)

################################################################################
######################## Call proportion #######################################
################################################################################

#### EXPLORE DATA ####
# how many observations of each sex in FLUT & CTRL?
table(CC_prop_data$TREATMENT, CC_prop_data$SEX) #
sample_df <- CC_prop_data[!duplicated(CC_prop_data$ID),] 
table(sample_df$SEX)
table(sample_df$TREATMENT)

# test for normality of data
ggqqplot(CC_prop_data$Avg_CC_prop)
shapiro.test(CC_prop_data$Avg_CC_prop)

# use a histogram to check data distribution
ggplot(CC_prop_data, aes(Avg_CC_prop)) +
  geom_histogram(fill= "white", col= "black", binwidth = 0.1) +
  facet_wrap(~TREATMENT) +
  theme_bw()
# by sex and treatment
ggplot(CC_prop_data, aes(Avg_CC_prop)) +
  geom_histogram(fill= "white", col= "black", binwidth = 0.1) +
  facet_wrap(TREATMENT~SEX) +
  theme_bw()

# check for correlation of vars
cor(sapply(CC_prop_data[, c("Avg_CC_prop", "REC_AGE_D", "SEX", "TREATMENT",
                             "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                             "CONDITION", "CONDITION_PER")], as.numeric))

# plot a simple overview for treatment/age/sex
ggplot(CC_prop_data, aes(REC_AGE_D, Avg_CC_prop)) +
  geom_point() +
  facet_grid(TREATMENT~SEX) +
  theme_bw()

#check if there is an obvious treatment effect
plot(CC_prop_data$Avg_CC_prop~CC_prop_data$TREATMENT)
fit <- lm(Avg_CC_prop~TREATMENT * REC_AGE_D * SEX, CC_prop_data) 
coef(fit)
S(fit)
rm(fit)

# overview plot 
plot_CC_prop <- ggplot(CC_prop_data, aes( REC_AGE_D, Avg_CC_prop)) +
  geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
  geom_smooth(method = "lm", col= 1) +
  labs(y= "CC call length") +
  facet_wrap(~ TREATMENT:SEX) +
  theme_bw()
annotate_figure(plot_CC_prop, 
                top = text_grob("Avg_CC_prop", 
                                color = "red", face = "bold", size = 14))
rm(plot_CC_prop)

####OUTLIER CHECK ####
ggplot(CC_prop_data) +
  aes(x = "", y = Avg_CC_prop) +
  geom_boxplot(fill = "#0c4c8a") +
  scale_color_okabe_ito()+
  theme_clean()+
  facet_grid(~SEX:TREATMENT)
# get the outliers (including indices in df)
out <- boxplot.stats(CC_prop_data$Avg_CC_prop)$out
(out_ind <- which(CC_prop_data$Avg_CC_prop %in% c(out)))# no outliers

# which datapoints seem to be outliers?
CC_prop_data[out_ind, ]
CC_prop_data[out_ind, ]$Avg_CC_prop

# Rosner Test for outliers to double check
# k= number of suspected outliers
test <- rosnerTest(CC_prop_data$Avg_CC_prop,
                   k = length(out))
test

#### BAYES MODELS ##############################################################
#### 1) Test random effects and nested level
B_CC_prop_ran_00 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ 1 + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_CC_prop_ran_00")
B_CC_prop_ran_01 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ 1 + (1|LITTER_CODE), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_CC_prop_ran_01")
B_CC_prop_ran_02 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ 1 + (1|REC_GROUP_REF), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_CC_prop_ran_02")
B_CC_prop_ran_03 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ 1 + (1|LITTER_CODE/ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_CC_prop_ran_03")
B_CC_prop_ran_04 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_CC_prop_ran_04")

loo(B_CC_prop_ran_00, B_CC_prop_ran_01, B_CC_prop_ran_02, B_CC_prop_ran_03, B_CC_prop_ran_04, moment_match = F)

###### ICC #####################################################################
performance::variance_decomposition(B_CC_prop_ran_00)
performance::r2_bayes(B_CC_prop_ran_00) # marginal = fixed effects only

performance::variance_decomposition(B_CC_prop_ran_01)
performance::r2_bayes(B_CC_prop_ran_01)

performance::variance_decomposition(B_CC_prop_ran_02)
performance::r2_bayes(B_CC_prop_ran_02) 

performance::variance_decomposition(B_CC_prop_ran_03)
performance::r2_bayes(B_CC_prop_ran_03)

performance::variance_decomposition(B_CC_prop_ran_04)
performance::r2_bayes(B_CC_prop_ran_04)

rm(B_CC_prop_ran_00, B_CC_prop_ran_01, B_CC_prop_ran_02, B_CC_prop_ran_03, B_CC_prop_ran_04 )

#### 2) Define all (biologically plausible) models ####
B_CC_prop_00 <- brms::brm(formula = Sum_CC | trials(Total_calls)~ +1 + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_CC_prop_00")
# only one predictor
B_CC_prop_01 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_CC_prop_01")
# additive predictors - no interactions
B_CC_prop_02 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_CC_prop_02")
B_CC_prop_03 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + SEX + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_CC_prop_03")
B_CC_prop_04 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_CC_prop_04")
B_CC_prop_04_a <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT + REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                             save_pars = save_pars(all=T), cores=4, file = "B_CC_prop_04_a")
B_CC_prop_05 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER+  (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_CC_prop_05")
# interactions
B_CC_prop_06 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT *  REC_AGE_D + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_CC_prop_06")
B_CC_prop_07 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT *  REC_AGE_D * SEX + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, init_r=0.11, file = "B_CC_prop_07")
B_CC_prop_08 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, init_r=0.11, file = "B_CC_prop_08")
B_CC_prop_08_a <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT * REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                             save_pars = save_pars(all=T), cores=4, init_r=0.005, file = "B_CC_prop_08_a")
B_CC_prop_09 <- brms::brm(formula =Sum_CC | trials(Total_calls) ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|ID), data = CC_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, init_r=0.005, file = "B_CC_prop_09")
#### CHECK MODELS ####
summary(B_CC_prop_00) 
#posterior_summary(B_CC_prop_00)
plot(B_CC_prop_00) 
pp_check(B_CC_prop_00, ndraws = 100) 
model_loo <- loo(B_CC_prop_00, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_01) 
#posterior_summary(B_CC_prop_01)
plot(B_CC_prop_01) 
pp_check(B_CC_prop_01, ndraws = 100) 
model_loo <- loo(B_CC_prop_01, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_02) 
#posterior_summary(B_CC_prop_02)
plot(B_CC_prop_02) 
pp_check(B_CC_prop_02, ndraws = 100) 
model_loo <- loo(B_CC_prop_02, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_03) 
#posterior_summary(B_CC_prop_03)
plot(B_CC_prop_03) 
pp_check(B_CC_prop_03, ndraws = 100) 
model_loo <- loo(B_CC_prop_03, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_04) 
#posterior_summary(B_CC_prop_04)
plot(B_CC_prop_04) 
pp_check(B_CC_prop_04, ndraws = 100) 
model_loo <- loo(B_CC_prop_04, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_04_a) 
#posterior_summary(B_CC_prop_04_a)
plot(B_CC_prop_04_a) 
pp_check(B_CC_prop_04_a, ndraws = 100) 
model_loo <- loo(B_CC_prop_04_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_05) 
#posterior_summary(B_CC_prop_05)
plot(B_CC_prop_05) 
pp_check(B_CC_prop_05, ndraws = 100) 
model_loo <- loo(B_CC_prop_05, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_06) 
#posterior_summary(B_CC_prop_06)
plot(B_CC_prop_06) 
pp_check(B_CC_prop_06, ndraws = 100) 
model_loo <- loo(B_CC_prop_06, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_07) 
#posterior_summary(B_CC_prop_07)
plot(B_CC_prop_07) 
pp_check(B_CC_prop_07, ndraws = 100) 
model_loo <- loo(B_CC_prop_07, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_08) 
#posterior_summary(B_CC_prop_08)
plot(B_CC_prop_08) 
pp_check(B_CC_prop_08, ndraws = 100) 
model_loo <- loo(B_CC_prop_08, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_08_a) 
#posterior_summary(B_CC_prop_08_a)
plot(B_CC_prop_08_a) 
pp_check(B_CC_prop_08_a, ndraws = 100) 
model_loo <- loo(B_CC_prop_08_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

summary(B_CC_prop_09) 
#posterior_summary(B_CC_prop_09)
plot(B_CC_prop_09) 
pp_check(B_CC_prop_09, ndraws = 100) 
model_loo <- loo(B_CC_prop_09, save_psis = TRUE, cores=4)
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
loo(B_CC_prop_00, B_CC_prop_01, B_CC_prop_02, B_CC_prop_03, B_CC_prop_04, 
    B_CC_prop_04_a, B_CC_prop_05, B_CC_prop_06, B_CC_prop_07, B_CC_prop_08, 
    B_CC_prop_08_a, B_CC_prop_09)

#### MOMENT MATCHING ####
loo(B_CC_prop_00, B_CC_prop_01, B_CC_prop_02, B_CC_prop_03, B_CC_prop_04, 
    B_CC_prop_04_a, B_CC_prop_05, B_CC_prop_06, B_CC_prop_07, B_CC_prop_08, 
    B_CC_prop_08_a, B_CC_prop_09, moment_match = T)

#### BEST MODEL ####
summary(B_CC_prop_09)

mcmc_plot(B_CC_prop_09, type = 'intervals')
mcmc_plot(B_CC_prop_09, type= 'acf')

pd <- p_direction(B_CC_prop_09)
plot(pd)

diagnostic_posterior(B_CC_prop_09)
(cred_int <- ci(B_CC_prop_09, ci=c(.89, .95), effects=c('fixed')))

# convert from logit to proportional value again? same as probabilities: prob = odds / (1 + odds)
inv_logit_scaled(fixef(B_CC_prop_09)) 

exp(fixef(B_CC_prop_09)) # ODDS RATIO

# ROPE: different definition for logistic regression: -0.18, 0.18
rope_rang <- c((sd(CC_prop_data$Avg_CC_prop)*0.18)*-1, sd(CC_prop_data$Avg_CC_prop)*0.18)

describe_posterior(
  B_CC_prop_09,
  effects = "fixed", # fixed vs all
  component = "all",
  rope_range = rope_rang,
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

equivalence_test(B_CC_prop_09, range=rope_rang)

# convert coeff to odd ratio
post_best <- describe_posterior(
  B_CC_prop_09,
  effects = "fixed",
  component = "all",
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

post_best$Median <- exp(post_best$Median)
post_best$CI_high <- exp(post_best$CI_high)
post_best$CI_low <- exp(post_best$CI_low)
post_best <- post_best[-c(4,5)]
post_best

loo_R2(B_CC_prop_09) # loo adjusted R2
performance::variance_decomposition(B_CC_prop_09)

## CALCULATE ESTIMATES FOR INTERACTIONS ####

# treatment * sex * helper * weight
weight_vars <- c(mean(CC_prop_data$WEIGHT_DIFF_PER)-sd(CC_prop_data$WEIGHT_DIFF_PER), 
                 mean(CC_prop_data$WEIGHT_DIFF_PER), 
                 mean(CC_prop_data$WEIGHT_DIFF_PER) + sd(CC_prop_data$WEIGHT_DIFF_PER))

(treat_sex_helper_weight <- emtrends(B_CC_prop_09, specs = c('TREATMENT', 'SEX', 'WEIGHT_DIFF_PER'), var = "Helper_Pup_REC",
                              at = list(WEIGHT_DIFF_PER = weight_vars)))


pd(treat_sex_helper_weight)

#### PLOTS: BEST MODEL ####
#Plot all conditions
conditional_effects(B_CC_prop_09, re_formula = NULL, robust=TRUE)

# AGE
min(CC_prop_data$REC_AGE_D)#31
max(CC_prop_data$REC_AGE_D)#130

age_dist <- B_CC_prop_09 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(CC_prop_data$TREATMENT),
                                    SEX = levels(CC_prop_data$SEX),
                                    Helper_Pup_REC = mean(CC_prop_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(CC_prop_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5),
                                    Total_calls=1), 
              re_formula = NULL, allow_new_levels=T)

age_dist$Avg_CC_prop <- age_dist$.epred
#age only
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_CC_prop)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = CC_prop_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Contact call proportion") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  guides(fill='none')+
  theme_clean()
# age & treatment
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_CC_prop, color=TREATMENT, fill=TREATMENT)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = CC_prop_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Contact call proportion") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  theme_clean() 

# TREATMENT*SEX*HELPER*WEIGHT
min(CC_prop_data$Helper_Pup_REC)#0.33
max(CC_prop_data$Helper_Pup_REC)#16

weight_vars <- c(mean(CC_prop_data$WEIGHT_DIFF_PER)-sd(CC_prop_data$WEIGHT_DIFF_PER), 
                 mean(CC_prop_data$WEIGHT_DIFF_PER), 
                 mean(CC_prop_data$WEIGHT_DIFF_PER) + sd(CC_prop_data$WEIGHT_DIFF_PER))

# labels for graphs
treatment_sex <- c(
  "CTRL:F" = "Control female",
  "CTRL:M" = "Control male",
  "FLUT:F" = "Flutamide female",
  "FLUT:M" = "Flutamide male"
)

treat_dist <- B_CC_prop_09 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(CC_prop_data$TREATMENT),
                                    SEX = levels(CC_prop_data$SEX),
                                    Helper_Pup_REC = seq(floor(min(CC_prop_data$Helper_Pup_REC)), 
                                                         ceiling(max(CC_prop_data$Helper_Pup_REC)), by=1),
                                    WEIGHT_DIFF_PER = weight_vars,
                                    REC_AGE_D = mean(CC_prop_data$REC_AGE_D),
                                    Total_calls=1), 
              re_formula = NULL, allow_new_levels=T)

treat_dist$TREATMENT <- as.factor(treat_dist$TREATMENT)
treat_dist$SEX <- as.factor(treat_dist$SEX)
treat_dist$WEIGHT_DIFF_PER <- as.factor(treat_dist$WEIGHT_DIFF_PER)
treat_dist$Avg_CC_prop <- treat_dist$.epred

ggplot(treat_dist, aes(x = Helper_Pup_REC, y = Avg_CC_prop, color=WEIGHT_DIFF_PER, fill=WEIGHT_DIFF_PER) ) +
    stat_lineribbon(.width = c(.95))+
    geom_point(data = CC_prop_data, size = 2, color='black',fill='black', alpha=0.7, shape='+') +   # raw data
    scale_color_okabe_ito(name = "Weight condition", labels=c('Poor', 'Normal', 'Good'))+
    scale_fill_okabe_ito(alpha=0.2, name = "Weight condition", labels=c('Poor', 'Normal', 'Good'))+
    labs(x = "Helper/Pup ratio", y = "Contact call proportion") +
    theme_clean()+
    facet_wrap(~TREATMENT:SEX, labeller = as_labeller(treatment_sex))

### cleanup ####
rm(B_CC_prop_00, B_CC_prop_01, B_CC_prop_02, B_CC_prop_03, B_CC_prop_04, 
   B_CC_prop_04_a, B_CC_prop_05, B_CC_prop_06, B_CC_prop_07, B_CC_prop_08, 
   B_CC_prop_08_a, B_CC_prop_09, CC_prop_data, CC_data, full_data)
