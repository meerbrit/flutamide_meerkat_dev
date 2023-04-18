#### Analysis script for Flutamide study: BEG (repeat) calls ###################
###### Bayesian Multilevel models, best model & Bayesian Model Averaging (BMA)##
############## BWalkenhorst 2022 ###############################################

#### SETUP ####

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(EnvStats) #rosnerTest for outliers
library(ggplot2) 
library(ggpubr) 
library(car)
library(tidyverse)
library(tidybayes) 
library(brms)
library(bayestestR) #e.g. diagnostic_posterior
library(sjPlot) #e.g. plot_model
library(ggokabeito) # color palette
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
BEG_data <- full_data[, c(1:18, 79, 105:113, 118, 120, 122, 125:126, 
                          128, 147, 156, 162, 166:169)]
# remove NAs
BEG_data <- na.omit(BEG_data) 
# how many observations of each sex in FLUT & CTRL?
table(BEG_data$TREATMENT, BEG_data$SEX) 

#create the different datasets for each response variable
BEG_len_data <- subset(BEG_data, BEG_avg_Len >0 & BEG_num >= MIN_BOUTS) # call length
BEG_len_data$Avg_BEG_mil <- (BEG_len_data$BEG_avg_Len)*1000 #convert avg_BEG_LEN to milliseconds
BEG_int_data <- subset(BEG_data, BEG_avg_Int >0 & BEG_num >= MIN_BOUTS) # call interval length
BEG_int_data$Avg_BEG_mil <- (BEG_int_data$BEG_avg_Int)*1000 #convert avg_BEG_LEN to milliseconds
BEG_rate_data <- subset(BEG_data, BEG_rate >0 & BEG_num >= MIN_BOUTS) # call rate
BEG_prop_data <- subset(BEG_data, nBouts >= MIN_BOUTS) #call proportion
BEG_prop_data$Total_calls <- BEG_prop_data$Sum_BEG+
  BEG_prop_data$Sum_DIG+BEG_prop_data$Sum_CC+BEG_prop_data$Sum_OTHER #save the number of total calls for calc later

# switch working directory to models folder
setwd("G:/My Drive/Uni/UZH/Projects/KMP_PUPS/Flutamide/scripts/BEG_models/")

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
## check for collinearity to see which variables can be modeled together
corvif(sapply(BEG_data[, 
                       c("REC_AGE_D", "SEX", "TREATMENT",
                         "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                         "CONDITION", "CONDITION_PER")], as.numeric))

corvif(sapply(BEG_data[, 
                       c("REC_AGE_D", "SEX", "TREATMENT",
                         "Helper_Pup_REC", "WEIGHT_DIFF_PER")], as.numeric))

################################################################################
######################## Call length ###########################################
################################################################################

#### EXPLORE DATA ####
# how many observations of each sex in FLUT & CTRL?
table(BEG_len_data$TREATMENT, BEG_len_data$SEX) 


# unique individuals
sample_df <- BEG_len_data[!duplicated(BEG_len_data$ID),] 
# xx individuals
table(sample_df$SEX)

table(sample_df$TREATMENT)

# test for normality of data
# visually using qq-plots 
ggqqplot(BEG_len_data$Avg_BEG_mil)
# or using Shapiro-Wilk normality test (p > 0.05 = normal distribution)
shapiro.test(BEG_len_data$Avg_BEG_mil) 

# use a histogram to check data distribution
ggplot(BEG_len_data, aes(Avg_BEG_mil)) +
  geom_histogram(fill= "white", col= "black", binwidth = 10) +
  facet_wrap(~TREATMENT) +
  theme_bw()
# by sex and treatment
ggplot(BEG_len_data, aes(Avg_BEG_mil)) +
  geom_histogram(fill= "white", col= "black", binwidth = 10) +
  facet_wrap(TREATMENT~SEX) +
  theme_bw()

# check for correlation of vars
cor(sapply(BEG_len_data[, c("Avg_BEG_mil", "REC_AGE_D", "SEX", "TREATMENT",
                            "Helper_Pup_REC", "WEIGHT_DIFF_PER")], as.numeric))

# plot a simple overview for treatment/age/sex
ggplot(BEG_len_data, aes(REC_AGE_D, Avg_BEG_mil)) +
  geom_point() +
  facet_grid(TREATMENT~SEX) +
  theme_bw()

#check if there is an obvious treatment effect
plot(BEG_len_data$Avg_BEG_mil~BEG_len_data$TREATMENT)
fit <- lm(Avg_BEG_mil~TREATMENT * REC_AGE_D * SEX, BEG_len_data) 
coef(fit)
S(fit)
rm(fit)

# overview plot 
plot_beg_len <- ggplot(BEG_len_data, aes( REC_AGE_D, Avg_BEG_mil)) +
  geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
  geom_smooth(method = "lm", col= 1) +
  labs(y= "BEG call length") +
  facet_wrap(~ TREATMENT:SEX) +
  theme_bw()
annotate_figure(plot_beg_len, 
                top = text_grob("Avg_BEG_mil", 
                color = "red", face = "bold", size = 14))
rm(plot_beg_len)

####OUTLIER CHECK ####
ggplot(BEG_len_data) +
  aes(x = "", y = Avg_BEG_mil) +
  geom_boxplot(fill = "#0c4c8a") +
  scale_color_okabe_ito()+
  theme_clean()+
  facet_grid(~SEX:TREATMENT)
# get the outliers (including indices in df)
out <- boxplot.stats(BEG_len_data$Avg_BEG_mil)$out
(out_ind <- which(BEG_len_data$Avg_BEG_mil %in% c(out)))

# which datapoints seem to be outliers?
BEG_len_data[out_ind, ]
BEG_len_data[out_ind, ]$Avg_BEG_mil

# Rosner Test for outliers to double check
# k= number of suspected outliers
test <- rosnerTest(BEG_len_data$Avg_BEG_mil,
                   k = length(out))
test

#### BAYES MODELS ##############################################################

#### 1) Test random effects and nested level using intercept only models ####
B_BEG_len_ran_00 <- brms::brm(formula =Avg_BEG_mil ~ 1 + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=3, file="B_BEG_len_ran_00")
B_BEG_len_ran_01 <- brms::brm(formula =Avg_BEG_mil ~ 1 + (1|LITTER_CODE), data = BEG_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file="B_BEG_len_ran_01")
B_BEG_len_ran_02 <- brms::brm(formula =Avg_BEG_mil ~ 1 + (1|REC_GROUP_REF), data = BEG_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file="B_BEG_len_ran_02")
B_BEG_len_ran_03 <- brms::brm(formula =Avg_BEG_mil ~ 1 + (1|LITTER_CODE/ID), data = BEG_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                              save_pars = save_pars(all=T), cores=3, file="B_BEG_len_ran_03")
B_BEG_len_ran_04 <- brms::brm(formula =Avg_BEG_mil ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = BEG_len_data, family = student(link='identity'),
                              chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99), 
                              save_pars = save_pars(all=T), cores=3, file="B_BEG_len_ran_04")

loo(B_BEG_len_ran_00, B_BEG_len_ran_01, B_BEG_len_ran_02, B_BEG_len_ran_03, B_BEG_len_ran_04, moment_match = F)
  
###### ICC ####
performance::variance_decomposition(B_BEG_len_ran_00)
performance::r2_bayes(B_BEG_len_ran_00) 

performance::variance_decomposition(B_BEG_len_ran_01)
performance::r2_bayes(B_BEG_len_ran_01) 

performance::variance_decomposition(B_BEG_len_ran_02)
performance::r2_bayes(B_BEG_len_ran_02)

performance::variance_decomposition(B_BEG_len_ran_03)
performance::r2_bayes(B_BEG_len_ran_03)

performance::variance_decomposition(B_BEG_len_ran_04)
performance::r2_bayes(B_BEG_len_ran_04) 

rm(B_BEG_len_ran_00, B_BEG_len_ran_01, B_BEG_len_ran_02, B_BEG_len_ran_03, B_BEG_len_ran_04)

#### 2) Define all (biologically plausible) models ####
B_BEG_len_00 <- brms::brm(formula =Avg_BEG_mil ~ +1 + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_00")
#only one predictor
B_BEG_len_01 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_01")
# additive predictors - no interactions
B_BEG_len_02 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT +  REC_AGE_D + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_02")

B_BEG_len_03 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT +  REC_AGE_D + SEX + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_03")

B_BEG_len_04 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_04")

B_BEG_len_04_a <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT +  REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_04_a")

B_BEG_len_05 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_05")
# interactions
B_BEG_len_06 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT *  REC_AGE_D + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_06")

B_BEG_len_07 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT *  REC_AGE_D * SEX + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_07")

B_BEG_len_08 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_08")

B_BEG_len_08_a <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT *  REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_08_a")

B_BEG_len_09 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|ID), data = BEG_len_data, family = student(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                          save_pars = save_pars(all=T), cores=4, file="B_BEG_len_09")

#### CHECK MODELS ####
summary(B_BEG_len_00)
posterior_summary(B_BEG_len_00)
plot(B_BEG_len_00) 
pp_check(B_BEG_len_00, ndraws = 100) 
# Pareto-k diagnostic can be useful to identify problematic point(s)
model_loo <- loo(B_BEG_len_00, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_00, type = 'intervals', prob = 0.89, prob_outer=0.95)
mcmc_plot(B_BEG_len_00, type = 'acf') # autocorrelation

summary(B_BEG_len_01) 
posterior_summary(B_BEG_len_01)
plot(B_BEG_len_01) 
pp_check(B_BEG_len_01, ndraws = 100) 
model_loo <- loo(B_BEG_len_01, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_01, type = 'intervals')
mcmc_plot(B_BEG_len_01, type = 'acf') 

summary(B_BEG_len_02) 
posterior_summary(B_BEG_len_02)
plot(B_BEG_len_02) 
pp_check(B_BEG_len_02, ndraws = 100) 
model_loo <- loo(B_BEG_len_02, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_02, type = 'intervals')
mcmc_plot(B_BEG_len_02, type = 'acf') 

summary(B_BEG_len_03)
posterior_summary(B_BEG_len_03)
plot(B_BEG_len_03)
pp_check(B_BEG_len_03, ndraws = 100)
model_loo <- loo(B_BEG_len_03, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_03, type = 'intervals')
mcmc_plot(B_BEG_len_03, type = 'acf') 

summary(B_BEG_len_04)
posterior_summary(B_BEG_len_04)
plot(B_BEG_len_04) 
pp_check(B_BEG_len_04, ndraws = 100) 
model_loo <- loo(B_BEG_len_04, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_04, type = 'intervals')
mcmc_plot(B_BEG_len_04, type = 'acf') 

summary(B_BEG_len_04_a)
posterior_summary(B_BEG_len_04_a)
plot(B_BEG_len_04_a) 
pp_check(B_BEG_len_04_a, ndraws = 100) 
model_loo <- loo(B_BEG_len_04_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_04_a, type = 'intervals')
mcmc_plot(B_BEG_len_04_a, type = 'acf')

summary(B_BEG_len_05)
posterior_summary(B_BEG_len_05)
plot(B_BEG_len_05) 
pp_check(B_BEG_len_05, ndraws = 100) 
model_loo <- loo(B_BEG_len_05, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_05, type = 'intervals')
mcmc_plot(B_BEG_len_05, type = 'acf')

summary(B_BEG_len_06)
posterior_summary(B_BEG_len_06)
plot(B_BEG_len_06)
pp_check(B_BEG_len_06, ndraws = 100)
model_loo <- loo(B_BEG_len_06, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_06, type = 'intervals')
mcmc_plot(B_BEG_len_06, type = 'acf')

summary(B_BEG_len_07) 
posterior_summary(B_BEG_len_07)
plot(B_BEG_len_07) 
pp_check(B_BEG_len_07, ndraws = 100) 
model_loo <- loo(B_BEG_len_07, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_07, type = 'intervals')
mcmc_plot(B_BEG_len_07, type = 'acf')

summary(B_BEG_len_08) 
posterior_summary(B_BEG_len_08)
plot(B_BEG_len_08) 
pp_check(B_BEG_len_08, ndraws = 100)
model_loo <- loo(B_BEG_len_08, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_08, type = 'intervals')
mcmc_plot(B_BEG_len_08, type = 'acf')

summary(B_BEG_len_08_a) 
posterior_summary(B_BEG_len_08_a)
plot(B_BEG_len_08_a) 
pp_check(B_BEG_len_08_a, ndraws = 100) 
model_loo <- loo(B_BEG_len_08_a, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_08_a, type = 'intervals')
mcmc_plot(B_BEG_len_08_a, type = 'acf')

summary(B_BEG_len_09) 
posterior_summary(B_BEG_len_09)
plot(B_BEG_len_09) 
pp_check(B_BEG_len_09, ndraws = 100) 
model_loo <- loo(B_BEG_len_09, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(B_BEG_len_09, type = 'intervals')
mcmc_plot(B_BEG_len_09, type = 'acf')

rm(model_loo, k_rintercept, df)

#### MODEL COMPARISON ####
# leave-one-out-cross-validation (Vehtrai et al, 2017)
loo(B_BEG_len_00, B_BEG_len_01, B_BEG_len_02, B_BEG_len_03, B_BEG_len_04, 
    B_BEG_len_04_a, B_BEG_len_05, B_BEG_len_06, B_BEG_len_07, B_BEG_len_08, 
    B_BEG_len_08_a, B_BEG_len_09)

#If elpd difference < 4, the difference is small (Sivula, Magnusson and Vehtari, 2020)).
# ELPD = theoretical expected log pointwise predictive density

#### MOMENT MATCHING ####
# Moment matching LOO can be used to reduce the number of high Pareto k’s faster than by refitting all problematic cases.
# Paananen, T., Piironen, J., Buerkner, P.-C., Vehtari, A. (2020)
loo(B_BEG_len_00, B_BEG_len_01, B_BEG_len_02, B_BEG_len_03, B_BEG_len_04, 
    B_BEG_len_04_a, B_BEG_len_05, B_BEG_len_06, B_BEG_len_07, B_BEG_len_08, 
    B_BEG_len_08_a, B_BEG_len_09, moment_match = T)

#### BEST MODEL ####
summary(B_BEG_len_05)
mcmc_plot(B_BEG_len_05, type = 'intervals',  prob = 0.89, prob_outer=0.95)
mcmc_plot(B_BEG_len_05, type='acf')

diagnostic_posterior(B_BEG_len_05)
(cred_int <- ci(B_BEG_len_05, ci=c(.89, .95), effects=c('fixed')))

describe_posterior(
  B_BEG_len_05,
  effects = "fixed", # all for including random effects
  component = "all",
  rope_range = rope_range(B_BEG_len_05),
  test = c("p_direction", "p_significance", 'rope'),
  centrality = "all", 
  dispersion = TRUE
)

equivalence_test(B_BEG_len_05)

pd <- p_direction(B_BEG_len_05)
plot(pd)

loo_R2(B_BEG_len_05)

performance::variance_decomposition(B_BEG_len_05)

#### PLOTS: BEST MODEL ####
# plot all predictors
cond_eff <- conditional_effects(B_BEG_len_05, re_formula = NULL,  robust=TRUE) # include all random effects (vs NA), robust uses median
plot(cond_eff)

# TREATMENT:
treatment_dist <- B_BEG_len_05 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_len_data$TREATMENT),
                                    SEX = levels(BEG_len_data$SEX),
                                    Helper_Pup_REC = mean(BEG_len_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(BEG_len_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = mean(BEG_len_data$REC_AGE_D)), 
              re_formula = NULL, allow_new_levels=T)

treatment_dist$Avg_BEG_mil <- treatment_dist$.epred

ggplot(treatment_dist, aes(x = TREATMENT, y = Avg_BEG_mil, color= TREATMENT), fill=TREATMENT) +
  stat_boxplot(geom='errorbar', width=.8)+
  geom_boxplot(outlier.shape = NA)+
  scale_color_okabe_ito(name = "Treatment", labels=c("Control", "Flutamide"))+
  scale_fill_okabe_ito(alpha=0.3)+
  geom_point(data = BEG_len_data, size = 2, aes(color=TREATMENT), alpha=0.5) +   # raw data
  labs(x = "Treatment", y = "Begging call length (ms)") +
  scale_x_discrete(labels=c('Control', 'Flutamide'))+
  theme_clean()+
  guides(color='none', fill='none')

# AGE:
age_dist <- B_BEG_len_05 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_len_data$TREATMENT),
                                    SEX = levels(BEG_len_data$SEX),
                                    Helper_Pup_REC = mean(BEG_len_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(BEG_len_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5)), 
              re_formula = NULL, allow_new_levels=T)

age_dist$Avg_BEG_mil <- age_dist$.epred
#age only
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_BEG_mil)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = BEG_len_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Age (days)", y = "Begging call length (ms)") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
  guides(fill='none')+
  theme_clean()
# age & treatment
ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_BEG_mil, color=TREATMENT, fill=TREATMENT)) +  
    stat_lineribbon(.width=.95)+
    geom_point(data = BEG_len_data, size = 2, aes(color=TREATMENT)) +   # raw data
    scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
    scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
    labs(x = "Age (days)", y = "Begging call length (ms)") +
    scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
    theme_clean()

# crossing CIs:
#Helper
helper_dist <- B_BEG_len_05 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_len_data$TREATMENT),
                                    SEX = levels(BEG_len_data$SEX),
                                    Helper_Pup_REC = seq(round(min(BEG_prop_data$Helper_Pup_REC)), 
                                                         round(max(BEG_prop_data$Helper_Pup_REC)), by=1),
                                    WEIGHT_DIFF_PER = mean(BEG_len_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = mean(BEG_len_data$REC_AGE_D)), 
              re_formula = NULL, allow_new_levels=T)

helper_dist$Avg_BEG_mil <- helper_dist$.epred
#helper only
ggplot(helper_dist, aes(x = Helper_Pup_REC, y = Avg_BEG_mil)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = BEG_len_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Helper/Pup ratio", y = "Begging call length (ms)") +
  guides(fill='none')+
  theme_clean()
# helper & treatment
ggplot(helper_dist, aes(x = Helper_Pup_REC, y = Avg_BEG_mil, color=TREATMENT, fill=TREATMENT)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = BEG_len_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Helper/Pup ratio", y = "Begging call length (ms)") +
  theme_clean()

#Weight
min(BEG_len_data$WEIGHT_DIFF_PER)
max(BEG_len_data$WEIGHT_DIFF_PER)

weight_dist <- B_BEG_len_05 %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_len_data$TREATMENT),
                                    SEX = levels(BEG_len_data$SEX),
                                    Helper_Pup_REC = mean(BEG_len_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = seq(floor(min(BEG_len_data$WEIGHT_DIFF_PER)),
                                                          ceiling(max(BEG_len_data$WEIGHT_DIFF_PER)), by=2),
                                    REC_AGE_D = mean(BEG_len_data$REC_AGE_D)), 
              re_formula = NULL, allow_new_levels=T)

weight_dist$Avg_BEG_mil <- weight_dist$.epred
#weight only
ggplot(weight_dist, aes(x = WEIGHT_DIFF_PER, y = Avg_BEG_mil)) +
  stat_lineribbon(.width = c(.95))+
  geom_point(data = BEG_len_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Weight offset (%)", y = "Begging call length (ms)") +
  guides(fill='none')+
  theme_clean()
# weight & treatment
ggplot(weight_dist, aes(x = WEIGHT_DIFF_PER, y = Avg_BEG_mil, color=TREATMENT, fill=TREATMENT)) +  
  stat_lineribbon(.width=.95)+
  geom_point(data = BEG_len_data, size = 2, aes(color=TREATMENT)) +   # raw data
  scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
  labs(x = "Weight offset (%)", y = "Begging call length (ms)") +
  theme_clean()

# AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA) ####
# no new data, just get the information of fitted values
post_avg <- posterior_average(B_BEG_len_05, B_BEG_len_04_a, B_BEG_len_06, 
                              B_BEG_len_02, B_BEG_len_04, B_BEG_len_03, 
                              B_BEG_len_07,
                              weights = "stacking", 
                              missing = 0)

posterior_summary(post_avg)

ci(post_avg, ci=c(.89, .95)) 

describe_posterior(
  post_avg,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_BEG_len_05),
  test = c("p_direction", "p_significance", "rope"),
  dispersion = TRUE
)

rm(post_avg)

#cleanup ####
rm(B_BEG_len_00, B_BEG_len_01, B_BEG_len_02, B_BEG_len_03, B_BEG_len_04, 
   B_BEG_len_04_a, B_BEG_len_05, B_BEG_len_06, B_BEG_len_07, B_BEG_len_08, 
   B_BEG_len_08_a, B_BEG_len_09, BEG_len_data)

 ################################################################################
 ######################## Call interval #########################################
 ################################################################################
 # 1) explore data ####
 # how many observations of each sex in FLUT & CTRL?
 table(BEG_int_data$TREATMENT, BEG_int_data$SEX) 

 # test for normality of data
 ggqqplot(BEG_int_data$BEG_avg_Int)
 shapiro.test(BEG_int_data$BEG_avg_Int) 
 
 # use a histogram to check data distribution
 ggplot(BEG_int_data, aes(BEG_avg_Int)) +
   geom_histogram(fill= "white", col= "black", binwidth = 0.01) +
   facet_wrap(~TREATMENT) +
   theme_bw()
 # by sex and treatment
 ggplot(BEG_int_data, aes(BEG_avg_Int)) +
   geom_histogram(fill= "white", col= "black", binwidth = 0.01) +
   facet_wrap(TREATMENT~SEX) +
   theme_bw()
 
 # check for correlation of vars
 cor(sapply(BEG_int_data[, c("BEG_avg_Int", "REC_AGE_D", "SEX", "TREATMENT",
                             "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                             "CONDITION", "CONDITION_PER")], as.numeric))
 
 # plot a simple overview for treatment/age/sex
 ggplot(BEG_int_data, aes(REC_AGE_D, BEG_avg_Int)) +
   geom_point() +
   facet_grid(TREATMENT~SEX) +
   theme_bw()
 
 #check if there is an obvious treatment effect
 plot(BEG_int_data$BEG_avg_Int~BEG_int_data$TREATMENT)
 fit <- lm(BEG_avg_Int~TREATMENT * REC_AGE_D * SEX, BEG_int_data) 
 coef(fit)
 S(fit)
 rm(fit)
 
 # overview plot ########
 plot_beg_int <- ggplot(BEG_int_data, aes( REC_AGE_D, Avg_BEG_mil)) +
   geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
   geom_smooth(method = "lm", col= 1) +
   labs(y= "BEG call interval length") +
   facet_wrap(~ TREATMENT:SEX) +
   theme_bw()
 annotate_figure(plot_beg_int, 
                 top = text_grob("BEG_avg_Int", 
                                 color = "red", face = "bold", size = 14))
 
 rm(plot_beg_int)
 
 ####OUTLIER CHECK ####
 ggplot(BEG_int_data) +
   aes(x = "", y = Avg_BEG_mil) +
   geom_boxplot(fill = "#0c4c8a") +
   scale_color_okabe_ito()+
   theme_clean()+
   facet_grid(~SEX:TREATMENT)
 # get the outliers (including indices in df)
 out <- boxplot.stats(BEG_int_data$Avg_BEG_mil)$out
 (out_ind <- which(BEG_int_data$Avg_BEG_mil %in% c(out)))
 
 # which datapoints seem to be outliers?
 BEG_int_data[out_ind, ]
 BEG_int_data[out_ind, ]$Avg_BEG_mil
 
 # Rosner Test for outliers to double check
 # k= number of suspected outliers
 test <- rosnerTest(BEG_int_data$Avg_BEG_mil,
                    k = length(out))
 test 
 
 #### BAYES MODELS ##############################################################
 #### 1) Test random effects and nested level
 B_BEG_int_ran_00 <- brms::brm(formula =Avg_BEG_mil ~ 1 + (1|ID), data = BEG_int_data, family = student(link='identity'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_BEG_int_ran_00")
 B_BEG_int_ran_01 <- brms::brm(formula =Avg_BEG_mil ~ 1 + (1|LITTER_CODE), data = BEG_int_data, family = student(link='identity'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_BEG_int_ran_01")
 B_BEG_int_ran_02 <- brms::brm(formula =Avg_BEG_mil ~ 1 + (1|REC_GROUP_REF), data = BEG_int_data, family = student(link='identity'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_BEG_int_ran_02")
 B_BEG_int_ran_03 <- brms::brm(formula =Avg_BEG_mil ~ 1 + (1|LITTER_CODE/ID), data = BEG_int_data, family = student(link='identity'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_BEG_int_ran_03")
 B_BEG_int_ran_04 <- brms::brm(formula =Avg_BEG_mil ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = BEG_int_data, family = student(link='identity'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_BEG_int_ran_04")
 
 loo(B_BEG_int_ran_00, B_BEG_int_ran_01, B_BEG_int_ran_02, B_BEG_int_ran_03, B_BEG_int_ran_04, moment_match = F)

 ###### ICC #####################################################################
 performance::variance_decomposition(B_BEG_int_ran_00)
 performance::r2_bayes(B_BEG_int_ran_00) 

 performance::variance_decomposition(B_BEG_int_ran_01)
 performance::r2_bayes(B_BEG_int_ran_01) 
 
 performance::variance_decomposition(B_BEG_int_ran_02)
 performance::r2_bayes(B_BEG_int_ran_02)

 performance::variance_decomposition(B_BEG_int_ran_03)
 performance::r2_bayes(B_BEG_int_ran_03) 
 
 performance::variance_decomposition(B_BEG_int_ran_04)
 performance::r2_bayes(B_BEG_int_ran_04) 

 rm(B_BEG_int_ran_00, B_BEG_int_ran_01, B_BEG_int_ran_02, B_BEG_int_ran_03, B_BEG_int_ran_04)
 
 #### 2) Define all (biologically plausible) models ####
 B_BEG_int_00 <- brms::brm(formula =Avg_BEG_mil ~ +1 + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_00")
 # only one predictor
 B_BEG_int_01 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_01")
 # additive predictors - no interactions
 B_BEG_int_02 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT +  REC_AGE_D + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_02")
 B_BEG_int_03 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT +  REC_AGE_D + SEX + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_03")
 B_BEG_int_04 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_04")
 B_BEG_int_04_a <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT +  REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_04_a")
 B_BEG_int_05 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_05")
 # interactions
 B_BEG_int_06 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT *  REC_AGE_D + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_06")
 B_BEG_int_07 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT *  REC_AGE_D * SEX + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_07")
 B_BEG_int_08 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_08")
 B_BEG_int_08_a <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT *  REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_08_a")
 B_BEG_int_09 <- brms::brm(formula =Avg_BEG_mil ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|ID), data = BEG_int_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_int_09")
 
 #### CHECK MODELS ####
 summary(B_BEG_int_00) 
 #posterior_summary(B_BEG_int_00)
 plot(B_BEG_int_00) 
 pp_check(B_BEG_int_00, ndraws = 100) 
 model_loo <- loo(B_BEG_int_00, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 ranef(B_BEG_int_00)
 
 summary(B_BEG_int_01) 
 #posterior_summary(B_BEG_int_01)
 plot(B_BEG_int_01) 
 pp_check(B_BEG_int_01, ndraws = 100) #
 model_loo <- loo(B_BEG_int_01, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_int_02) 
 #posterior_summary(B_BEG_int_02)
 plot(B_BEG_int_02) 
 pp_check(B_BEG_int_02, ndraws = 100) #
 model_loo <- loo(B_BEG_int_02, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_int_03) 
 #posterior_summary(B_BEG_int_03)
 plot(B_BEG_int_03) 
 pp_check(B_BEG_int_03, ndraws = 100) #
 model_loo <- loo(B_BEG_int_03, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_int_04) 
 #posterior_summary(B_BEG_int_04)
 plot(B_BEG_int_04) 
 pp_check(B_BEG_int_04, ndraws = 100) #
 model_loo <- loo(B_BEG_int_04, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_int_04_a) 
 #posterior_summary(B_BEG_int_04_a)
 plot(B_BEG_int_04_a) 
 pp_check(B_BEG_int_04_a, ndraws = 100) #
 model_loo <- loo(B_BEG_int_04_a, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_int_05) 
 #posterior_summary(B_BEG_int_05)
 plot(B_BEG_int_05) 
 pp_check(B_BEG_int_05, ndraws = 100) #
 model_loo <- loo(B_BEG_int_05, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_int_06) 
 #posterior_summary(B_BEG_int_06)
 plot(B_BEG_int_06) 
 pp_check(B_BEG_int_06, ndraws = 100) #
 model_loo <- loo(B_BEG_int_06, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_int_07) 
 #posterior_summary(B_BEG_int_07)
 plot(B_BEG_int_07) 
 pp_check(B_BEG_int_07, ndraws = 100) #
 model_loo <- loo(B_BEG_int_07, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_int_08) 
 #posterior_summary(B_BEG_int_08)
 plot(B_BEG_int_08) 
 pp_check(B_BEG_int_08, ndraws = 100) #
 model_loo <- loo(B_BEG_int_08, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_int_08_a) 
 #posterior_summary(B_BEG_int_08_a)
 plot(B_BEG_int_08_a) 
 pp_check(B_BEG_int_08_a, ndraws = 100) #
 model_loo <- loo(B_BEG_int_08_a, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_int_09) 
 #posterior_summary(B_BEG_int_09)
 plot(B_BEG_int_09) 
 pp_check(B_BEG_int_09, ndraws = 100) #
 model_loo <- loo(B_BEG_int_09, save_psis = TRUE, cores=4)
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
 loo(B_BEG_int_00, B_BEG_int_01, B_BEG_int_02, B_BEG_int_03, B_BEG_int_04,
     B_BEG_int_04_a, B_BEG_int_05, B_BEG_int_06, B_BEG_int_07, B_BEG_int_08, 
     B_BEG_int_08_a, B_BEG_int_09) 

 #### MOMENT MATCHING ####
 loo(B_BEG_int_00, B_BEG_int_01, B_BEG_int_02, B_BEG_int_03, B_BEG_int_04,
     B_BEG_int_04_a, B_BEG_int_05, B_BEG_int_06, B_BEG_int_07, B_BEG_int_08, 
     B_BEG_int_08_a, B_BEG_int_09, moment_match = T)
 

 #### BEST MODEL ####
 summary(B_BEG_int_04_a)
 
 mcmc_plot(B_BEG_int_04_a, type = 'intervals', prob = 0.89, prob_outer=0.95)
 
 diagnostic_posterior(B_BEG_int_04_a)
 (cred_int <- ci(B_BEG_int_04_a, ci=c(.89, .95), effects=c('fixed')))
 
 describe_posterior(
   B_BEG_int_04_a,
   effects = "fixed", #all for random effect
   component = "all",
   rope_range = rope_range(B_BEG_int_04_a),
   test = c("p_direction", "p_significance", "rope"),
   centrality = "all",
   dispersion = TRUE
 )

 equivalence_test(B_BEG_int_04_a)
 
 pd <- p_direction(B_BEG_int_04_a)
 plot(pd)
 
 loo_R2(B_BEG_int_04_a)
 performance::variance_decomposition(B_BEG_int_04_a)

 #### PLOTS BEST MODEL ####
 # plot all predictors
 cond_eff <- conditional_effects(B_BEG_int_04_a, re_formula = NULL, robust=TRUE)
 plot(cond_eff)
 
 #WEIGHT
 weight_dist <- B_BEG_int_04_a %>% 
   epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_int_data$TREATMENT),
                                     SEX = levels(BEG_len_data$SEX),
                                     Helper_Pup_REC = mean(BEG_int_data$Helper_Pup_REC),
                                     WEIGHT_DIFF_PER =
                                       seq(floor(min(BEG_int_data$WEIGHT_DIFF_PER)),
                                           ceiling(max(BEG_int_data$WEIGHT_DIFF_PER)), by=2),
                                     REC_AGE_D = mean(BEG_int_data$REC_AGE_D)),
               re_formula = NULL, allow_new_levels=T)
 
 weight_dist$Avg_BEG_mil <- weight_dist$.epred
 #weight only
 ggplot(weight_dist, aes(x = WEIGHT_DIFF_PER, y = Avg_BEG_mil)) +
     stat_lineribbon(.width=.95)+
     scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
     scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
     geom_point(data = BEG_int_data, size = 2, aes(color=TREATMENT)) + # raw data
     labs(x = "Weight offset (%)",
          y = "Begging call interval length (ms)")+
     guides(fill='none')+
     theme_clean()
 # weight & treatment
 ggplot(weight_dist, aes(x = WEIGHT_DIFF_PER, y = Avg_BEG_mil, color=TREATMENT, fill=TREATMENT)) +  
   stat_lineribbon(.width=.95)+
   geom_point(data = BEG_int_data, size = 2, aes(color=TREATMENT)) +   # raw data
   scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
   scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
   labs(x = "Weight offset (%)",
        y = "Begging call interval length (ms)") +
   theme_clean()
 
 # CI crosses 0:
 # TREATMENT:
 treatment_dist <- B_BEG_int_04_a %>% 
   epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_int_data$TREATMENT),
                                     SEX = levels(BEG_int_data$SEX),
                                     Helper_Pup_REC = mean(BEG_int_data$Helper_Pup_REC),
                                     WEIGHT_DIFF_PER = mean(BEG_int_data$WEIGHT_DIFF_PER),
                                     REC_AGE_D = mean(BEG_int_data$REC_AGE_D)), 
               re_formula = NULL, allow_new_levels=T)
 
 treatment_dist$Avg_BEG_mil <- treatment_dist$.epred
 
 ggplot(treatment_dist, aes(x = TREATMENT, y = Avg_BEG_mil, color= TREATMENT), fill=TREATMENT) +
   stat_boxplot(geom='errorbar', width=.8)+
   geom_boxplot(outlier.shape = NA)+
   scale_color_okabe_ito(name = "Treatment", labels=c("Control", "Flutamide"))+
   scale_fill_okabe_ito(alpha=0.3)+
   geom_point(data = BEG_int_data, size = 2, aes(color=TREATMENT), alpha=0.5) +   # raw data
   labs(x = "Treatment", y = "Begging call interval length (ms)") +
   scale_x_discrete(labels=c('Control', 'Flutamide'))+
   theme_clean()+
   guides(color='none', fill='none')

 # AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA) ####
 # no new data, just get the information of fitted values
 post_avg <- posterior_average(B_BEG_int_04_a, B_BEG_int_05, B_BEG_int_08_a, 
                               weights = "stacking", 
                               missing = 0)
 
 posterior_summary(post_avg)


 ci(post_avg, ci=c(.89, .95)) 
 
 describe_posterior(
   post_avg,
   effects = "fixed",
   component = "all",
   rope_range= rope_range(B_BEG_int_05),  # can use any of the models
   test = c("p_direction", "p_significance", "rope"),
   dispersion = TRUE
 )
 
 rm(post_avg)
 
 ### cleanup ####
 rm(B_BEG_int_00, B_BEG_int_01, B_BEG_int_02, B_BEG_int_03, B_BEG_int_04, 
    B_BEG_int_04_a, B_BEG_int_05, B_BEG_int_06, B_BEG_int_07, B_BEG_int_08, 
    B_BEG_int_08_a, B_BEG_int_09, BEG_int_data, BEG_int_new, BEG_int_new_all)
 
 ################################################################################
 ######################## Call rate #############################################
 ################################################################################
 # 1) explore data
 # how many observations of each sex in FLUT & CTRL?
 table(BEG_rate_data$TREATMENT, BEG_rate_data$SEX) #

 # test for normality of data
 ggqqplot(BEG_rate_data$BEG_rate)
 shapiro.test(BEG_rate_data$BEG_rate) 
 
 # use a histogram to check data distribution
 ggplot(BEG_rate_data, aes(BEG_rate)) +
   geom_histogram(fill= "white", col= "black", binwidth = 0.01) +
   facet_wrap(~TREATMENT) +
   theme_bw()
 # by sex and treatment
 ggplot(BEG_rate_data, aes(BEG_rate)) +
   geom_histogram(fill= "white", col= "black", binwidth = 0.01) +
   facet_wrap(TREATMENT~SEX) +
   theme_bw()
 
 # check for correlation of vars
 cor(sapply(BEG_rate_data[, c("BEG_rate", "REC_AGE_D", "SEX", "TREATMENT",
                             "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                             "CONDITION", "CONDITION_PER")], as.numeric))
 

 # plot a simple overview for treatment/age/sex
 ggplot(BEG_rate_data, aes(REC_AGE_D, BEG_rate)) +
   geom_point() +
   facet_grid(TREATMENT~SEX) +
   theme_bw()
 
 #check if there is an obvious treatment effect
 plot(BEG_rate_data$BEG_rate~BEG_rate_data$TREATMENT)
 fit <- lm(BEG_rate~TREATMENT * REC_AGE_D * SEX, BEG_rate_data) 
 coef(fit)
 S(fit)
 
 # overview plot ########
 plot_beg_rate <- ggplot(BEG_rate_data, aes( REC_AGE_D, BEG_rate)) +
   geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
   geom_smooth(method = "lm", col= 1) +
   labs(y= "BEG call rate") +
   facet_wrap(~ TREATMENT:SEX) +
   theme_bw()
 annotate_figure(plot_beg_rate, 
                 top = text_grob("BEG_rate", 
                                 color = "red", face = "bold", size = 14))
 rm(plot_beg_rate)
 
 ####OUTLIER CHECK ####
 ggplot(BEG_rate_data) +
   aes(x = "", y = BEG_rate) +
   geom_boxplot(fill = "#0c4c8a") +
   scale_color_okabe_ito()+
   theme_clean()+
   facet_grid(~SEX:TREATMENT)
 # get the outliers (including indices in df)
 out <- boxplot.stats(BEG_rate_data$BEG_rate)$out
 (out_ind <- which(BEG_rate_data$BEG_rate %in% c(out)))
 
 # which datapoints seem to be outliers?
 BEG_rate_data[out_ind, ]
 BEG_rate_data[out_ind, ]$BEG_rate
 
 # Rosner Test for outliers to double check
 # k= number of suspected outliers
 test <- rosnerTest(BEG_rate_data$BEG_rate,
                    k = length(out))
 test 
 
 #### BAYES MODELS ##############################################################
 #### 1) Test random effects and nested level
 B_BEG_rat_ran_00 <- brms::brm(formula =BEG_rate ~ 1 + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_BEG_rat_ran_00")
 B_BEG_rat_ran_01 <- brms::brm(formula =BEG_rate ~ 1 + (1|LITTER_CODE), data = BEG_rate_data, family = student(link='identity'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_BEG_rat_ran_01")
 B_BEG_rat_ran_02 <- brms::brm(formula =BEG_rate ~ 1 + (1|REC_GROUP_REF), data = BEG_rate_data, family = student(link='identity'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_BEG_rat_ran_02")
 B_BEG_rat_ran_03 <- brms::brm(formula =BEG_rate ~ 1 + (1|LITTER_CODE/ID), data = BEG_rate_data, family = student(link='identity'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_BEG_rat_ran_03")
 B_BEG_rat_ran_04 <- brms::brm(formula =BEG_rate ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = BEG_rate_data, family = student(link='identity'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file = "B_BEG_rat_ran_04")
 
 loo(B_BEG_rat_ran_00, B_BEG_rat_ran_01, B_BEG_rat_ran_02, B_BEG_rat_ran_03, B_BEG_rat_ran_04, moment_match = F)
 
 ###### ICC ####
 performance::variance_decomposition(B_BEG_rat_ran_00)
 performance::r2_bayes(B_BEG_rat_ran_00) 

 
 performance::variance_decomposition(B_BEG_rat_ran_01)
 performance::r2_bayes(B_BEG_rat_ran_01) 

 performance::variance_decomposition(B_BEG_rat_ran_02)
 performance::r2_bayes(B_BEG_rat_ran_02)
 
 performance::variance_decomposition(B_BEG_rat_ran_03)
 performance::r2_bayes(B_BEG_rat_ran_03) 
 
 performance::variance_decomposition(B_BEG_rat_ran_04)
 performance::r2_bayes(B_BEG_rat_ran_04) 
 
 rm(B_BEG_rat_ran_00, B_BEG_rat_ran_01, B_BEG_rat_ran_02, B_BEG_rat_ran_03, B_BEG_rat_ran_04)
 
 #### 2) Define all (biologically plausible) models ####
 B_BEG_rat_00 <- brms::brm(formula = BEG_rate~ +1 + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_00")
 # only one predictor
 B_BEG_rat_01 <- brms::brm(formula =BEG_rate ~ TREATMENT + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_01")
 # additive predictors - no interactions
 B_BEG_rat_02 <- brms::brm(formula =BEG_rate ~ TREATMENT +  REC_AGE_D + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_02")
 B_BEG_rat_03 <- brms::brm(formula =BEG_rate ~ TREATMENT +  REC_AGE_D + SEX + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_03")
 B_BEG_rat_04 <- brms::brm(formula =BEG_rate ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_04")
 B_BEG_rat_04_a <- brms::brm(formula =BEG_rate ~ TREATMENT +  REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_04_a")
 B_BEG_rat_05 <- brms::brm(formula =BEG_rate ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER+  (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_05")
 # interactions
 B_BEG_rat_06 <- brms::brm(formula =BEG_rate ~ TREATMENT *  REC_AGE_D + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_06")
 B_BEG_rat_07 <- brms::brm(formula =BEG_rate ~ TREATMENT *  REC_AGE_D * SEX + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_07")
 B_BEG_rat_08 <- brms::brm(formula =BEG_rate ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_08")
 B_BEG_rat_08_a <- brms::brm(formula =BEG_rate ~ TREATMENT *  REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_08_a")
 B_BEG_rat_09 <- brms::brm(formula =BEG_rate ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|ID), data = BEG_rate_data, family = student(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 15),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_rat_09")
 
 #### CHECK MODELS ####
 summary(B_BEG_rat_00) 
 #posterior_summary(B_BEG_rat_00)
 plot(B_BEG_rat_00) 
 pp_check(B_BEG_rat_00, ndraws = 100) 
 model_loo <- loo(B_BEG_rat_00, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_01) 
 #posterior_summary(B_BEG_rat_01)
 plot(B_BEG_rat_01) 
 pp_check(B_BEG_rat_01, ndraws = 100) #
 model_loo <- loo(B_BEG_rat_01, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_02) 
 #posterior_summary(B_BEG_rat_02)
 plot(B_BEG_rat_02) 
 pp_check(B_BEG_rat_02, ndraws = 100) #
 model_loo <- loo(B_BEG_rat_02, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_03) 
 #posterior_summary(B_BEG_rat_03)
 plot(B_BEG_rat_03) 
 pp_check(B_BEG_rat_03, ndraws = 100) #
 model_loo <- loo(B_BEG_rat_03, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_04) 
 #posterior_summary(B_BEG_rat_04)
 plot(B_BEG_rat_04) 
 pp_check(B_BEG_rat_04, ndraws = 100) #
 model_loo <- loo(B_BEG_rat_04, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_04_a) 
 #posterior_summary(B_BEG_rat_04_a)
 plot(B_BEG_rat_04_a) 
 pp_check(B_BEG_rat_04_a, ndraws = 100) 
 model_loo <- loo(B_BEG_rat_04_a, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_05) 
 #posterior_summary(B_BEG_rat_05)
 plot(B_BEG_rat_05) 
 pp_check(B_BEG_rat_05, ndraws = 100) #
 model_loo <- loo(B_BEG_rat_05, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_06) 
 #posterior_summary(B_BEG_rat_06)
 plot(B_BEG_rat_06) 
 pp_check(B_BEG_rat_06, ndraws = 100) #
 model_loo <- loo(B_BEG_rat_06, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_07) 
 #posterior_summary(B_BEG_rat_07)
 plot(B_BEG_rat_07) 
 pp_check(B_BEG_rat_07, ndraws = 100) #
 model_loo <- loo(B_BEG_rat_07, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_08) 
 #posterior_summary(B_BEG_rat_08)
 plot(B_BEG_rat_08) 
 pp_check(B_BEG_rat_08, ndraws = 100) #
 model_loo <- loo(B_BEG_rat_08, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_08_a) 
 #posterior_summary(B_BEG_rat_08_a)
 plot(B_BEG_rat_08_a) 
 pp_check(B_BEG_rat_08_a, ndraws = 100) #
 model_loo <- loo(B_BEG_rat_08_a, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_rat_09) 
 #posterior_summary(B_BEG_rat_09)
 plot(B_BEG_rat_09) 
 pp_check(B_BEG_rat_09, ndraws = 100) #
 model_loo <- loo(B_BEG_rat_09, save_psis = TRUE, cores=4)
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
 loo(B_BEG_rat_00, B_BEG_rat_01, B_BEG_rat_02, B_BEG_rat_03, B_BEG_rat_04,
     B_BEG_rat_04_a, B_BEG_rat_05, B_BEG_rat_06, B_BEG_rat_07, B_BEG_rat_08, 
     B_BEG_rat_08_a, B_BEG_rat_09)

 #### MOMENT MATCHING ####
 loo(B_BEG_rat_00, B_BEG_rat_01, B_BEG_rat_02, B_BEG_rat_03, B_BEG_rat_04,
     B_BEG_rat_04_a, B_BEG_rat_05, B_BEG_rat_06, B_BEG_rat_07, B_BEG_rat_08, 
     B_BEG_rat_08_a, B_BEG_rat_09, moment_match = T)
 
 #### BEST MODEL ####
 summary(B_BEG_rat_05)
 
 mcmc_plot(B_BEG_rat_05, type = 'intervals', prob = 0.89, prob_outer=0.95)
 mcmc_plot(B_BEG_rat_05, type = 'acf')
 
 diagnostic_posterior(B_BEG_rat_05)
 (cred_int <- ci(B_BEG_rat_05, ci=c(.89, .95), effects=c('fixed')))

 
 describe_posterior(
   B_BEG_rat_05,
   effects = "fixed", #all
   component = "all",
   rope_range = rope_range(B_BEG_rat_05),
   test = c("p_direction", "p_significance", "rope"),
   centrality = "all",
   dispersion = TRUE
 )
 
equivalence_test(B_BEG_rat_05)
 
pd <- p_direction(B_BEG_rat_05)
plot(pd)
 
loo_R2(B_BEG_rat_05)
performance::variance_decomposition(B_BEG_rat_05)

#### PLOTS BEST MODEL ####
 # plot all predictors
 cond_eff <- conditional_effects(B_BEG_rat_05, re_formula = NULL, robust=TRUE) # include all random effects (vs NA), robust uses median
 plot(cond_eff)
 
 # TREATMENT:
 treatment_dist <- B_BEG_rat_05 %>% 
   epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_rate_data$TREATMENT),
                                     SEX = levels(BEG_rate_data$SEX),
                                     Helper_Pup_REC = mean(BEG_rate_data$Helper_Pup_REC),
                                     WEIGHT_DIFF_PER = mean(BEG_rate_data$WEIGHT_DIFF_PER),
                                     REC_AGE_D = mean(BEG_rate_data$REC_AGE_D)), 
               re_formula = NULL, allow_new_levels=T)
 
 treatment_dist$BEG_rate <- treatment_dist$.epred
 
 ggplot(treatment_dist, aes(x = TREATMENT, y = BEG_rate, color= TREATMENT), fill=TREATMENT) +
   stat_boxplot(geom='errorbar', width=.8)+
   geom_boxplot(outlier.shape = NA)+
   scale_color_okabe_ito(name = "Treatment", labels=c("Control", "Flutamide"))+
   scale_fill_okabe_ito(alpha=0.3)+
   geom_point(data = BEG_len_data, size = 2, aes(color=TREATMENT), alpha=0.5) +   # raw data
   labs(x = "Treatment", y = "Begging call rate") +
   scale_x_discrete(labels=c('Control', 'Flutamide'))+
   theme_clean()+
   guides(color='none', fill='none')
 
 #WEIGHT 
 weight_dist <- B_BEG_rat_05 %>% 
   epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_rate_data$TREATMENT),
                                     SEX = levels(BEG_rate_data$SEX),
                                     Helper_Pup_REC = mean(BEG_rate_data$Helper_Pup_REC),
                                     WEIGHT_DIFF_PER =
                                       seq(floor(min(BEG_rate_data$WEIGHT_DIFF_PER)),
                                           ceiling(max(BEG_rate_data$WEIGHT_DIFF_PER)), by=2),
                                     REC_AGE_D = mean(BEG_rate_data$REC_AGE_D)),
               re_formula = NULL, allow_new_levels=T)
 
 weight_dist$BEG_rate <- weight_dist$.epred
 #weight only
 ggplot(weight_dist, aes(x = WEIGHT_DIFF_PER, y = BEG_rate)) +
   stat_lineribbon(.width=.95)+
   scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
   scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
   geom_point(data = BEG_rate_data, size = 2, aes(color=TREATMENT)) + # raw data
   labs(x = "Weight offset (%)",
        y = "Begging call rate")+
   guides(fill='none')+
   theme_clean()
 # weight & treatment
 ggplot(weight_dist, aes(x = WEIGHT_DIFF_PER, y = BEG_rate, color=TREATMENT, fill=TREATMENT)) +  
   stat_lineribbon(.width=.95)+
   geom_point(data = BEG_rate_data, size = 2, aes(color=TREATMENT)) +   # raw data
   scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
   scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
   labs(x = "Weight offset (%)",
        y = "Begging call rate") +
   theme_clean()
 
 # crossing CIs  
 #Helper/Pup ratio
 min(BEG_rate_data$Helper_Pup_REC)#0.33
 max(BEG_rate_data$Helper_Pup_REC)#16
 
 helper_dist <- B_BEG_rat_05 %>% 
   epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_rate_data$TREATMENT),
                                     SEX = levels(BEG_rate_data$SEX),
                                     Helper_Pup_REC = seq(round(min(BEG_prop_data$Helper_Pup_REC)), 
                                                          round(max(BEG_prop_data$Helper_Pup_REC)), by=1),
                                     WEIGHT_DIFF_PER = mean(BEG_rate_data$WEIGHT_DIFF_PER),
                                     REC_AGE_D = mean(BEG_rate_data$REC_AGE_D)),
               re_formula = NULL, allow_new_levels=T)
 
 helper_dist$BEG_rate <- helper_dist$.epred
 #helper only
 ggplot(helper_dist, aes(x = Helper_Pup_REC, y = BEG_rate)) +
   stat_lineribbon(.width=.95)+
   scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
   scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
   geom_point(data = BEG_rate_data, size = 2, aes(color=TREATMENT)) + # raw data
   labs(x = "Helper/Pup ratio",
        y = "Begging call rate")+
   guides(fill='none')+
   theme_clean()
 # helper & treatment
 ggplot(helper_dist, aes(x = Helper_Pup_REC, y = BEG_rate, color=TREATMENT, fill=TREATMENT)) +  
   stat_lineribbon(.width=.95)+
   geom_point(data = BEG_rate_data, size = 2, aes(color=TREATMENT)) +   # raw data
   scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
   scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
   labs(x = "Helper/Pup ratio",
        y = "Begging call rate") +
   theme_clean()
 
 # AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA) ####
 # no new data, just get the information of fitted values
 post_avg <- posterior_average(B_BEG_rat_05, B_BEG_rat_04_a, B_BEG_rat_08_a, 
                               weights = "stacking", 
                               missing = 0)
 
 posterior_summary(post_avg)
 
 ci(post_avg, ci=c(.89, .95)) 
 
 describe_posterior(
   post_avg,
   effects = "fixed",
   component = "all",
  # centrality = "all",
   rope_range= rope_range(B_BEG_rat_05),  # can use any of the models
   test = c("p_direction", "p_significance", "rope"),
   dispersion = TRUE
 )
 
 rm(post_avg)
 
 # cleanup ####
 rm(B_BEG_rat_00, B_BEG_rat_01, B_BEG_rat_02, B_BEG_rat_03, B_BEG_rat_04, 
    B_BEG_rat_04_a, B_BEG_rat_05, B_BEG_rat_06, B_BEG_rat_07, B_BEG_rat_08, 
    B_BEG_rat_08_a, B_BEG_rat_09, BEG_rate_data)
 
################################################################################
 ######################## Call proportion #######################################
 ################################################################################
 
 #### EXPLORE DATA ####
 # how many observations of each sex in FLUT & CTRL?
 table(BEG_prop_data$TREATMENT, BEG_prop_data$SEX) #
 
 # unique individuals
 sample_df <- BEG_prop_data[!duplicated(BEG_prop_data$ID),] 
 table(sample_df$SEX)
 table(sample_df$TREATMENT)
 
 # test for normality of data
 ggqqplot(BEG_prop_data$Avg_BEG_prop)
 shapiro.test(BEG_prop_data$Avg_BEG_prop) 
 
 ggplot(BEG_prop_data, aes(Avg_BEG_prop)) +
   geom_histogram(fill= "white", col= "black", binwidth = 0.01) +
   facet_wrap(~TREATMENT) +
   theme_bw()
 # by sex and treatment
 ggplot(BEG_prop_data, aes(Avg_BEG_prop)) +
   geom_histogram(fill= "white", col= "black", binwidth = 0.01) +
   facet_wrap(TREATMENT~SEX) +
   theme_bw()
 
 # check for correlation of vars
 cor(sapply(BEG_prop_data[, c("Avg_BEG_prop", "REC_AGE_D", "SEX", "TREATMENT",
                             "Helper_Pup_REC", "WEIGHT_DIFF_PER", "WEIGHT_DIFF",
                             "CONDITION", "CONDITION_PER")], as.numeric))

 # plot a simple overview for treatment/age/sex
 ggplot(BEG_prop_data, aes(REC_AGE_D, Avg_BEG_prop)) +
   geom_point() +
   facet_grid(TREATMENT~SEX) +
   theme_bw()
 
 #check if there is an obvious treatment effect
 plot(BEG_prop_data$Avg_BEG_prop~BEG_prop_data$TREATMENT)
 fit <- lm(Avg_BEG_prop~TREATMENT * REC_AGE_D * SEX, BEG_prop_data) 
 coef(fit)
 S(fit)
 rm(fit)
 
 # overview plot 
 plot_beg_prop <- ggplot(BEG_prop_data, aes( REC_AGE_D, Avg_BEG_prop)) +
   geom_point(aes(col= TREATMENT:SEX), show.legend= T) +
   geom_smooth(method = "lm", col= 1) +
   labs(y= "BEG proportion") +
   facet_wrap(~ TREATMENT:SEX) +
   theme_bw()
 annotate_figure(plot_beg_prop, 
                 top = text_grob("Avg_BEG_prop", 
                                 color = "red", face = "bold", size = 14))
 
 rm(plot_beg_prop)

 BEG_prop_data %>% 
   summarise(mean = mean(Avg_BEG_prop),
             variance = var(Avg_BEG_prop)) 
 
 ####OUTLIER CHECK ####
 ggplot(BEG_prop_data) +
   aes(x = "", y = Avg_BEG_prop) +
   geom_boxplot(fill = "#0c4c8a") +
   scale_color_okabe_ito()+
   theme_clean()+
   facet_grid(~SEX:TREATMENT)
 # get the outliers (including indices in df)
 out <- boxplot.stats(BEG_prop_data$Avg_BEG_prop)$out
 (out_ind <- which(BEG_prop_data$Avg_BEG_prop %in% c(out)))
 
 # which datapoints seem to be outliers?
 BEG_prop_data[out_ind, ]
 BEG_prop_data[out_ind, ]$Avg_BEG_prop
 
 # Rosner Test for outliers to double check
 # k= number of suspected outliers
 test <- rosnerTest(BEG_prop_data$Avg_BEG_prop,
                    k = length(out))
 test 

 #### BAYES MODELS ##############################################################
 #### 1) Test random effects and nested level
 ## zero inflated beta binomial shows best fit
 B_BEG_prop_ran_00 <- brms::brm(formula = Sum_BEG | trials(Total_calls) ~ 1 + (1|ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file="B_BEG_prop_ran_00")
 B_BEG_prop_ran_01 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ 1 + (1|LITTER_CODE), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file="B_BEG_prop_ran_02")
 B_BEG_prop_ran_02 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ 1 + (1|REC_GROUP_REF), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file="B_BEG_prop_ran_02")
 B_BEG_prop_ran_03 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ 1 + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file="B_BEG_prop_ran_03")
 B_BEG_prop_ran_04 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ 1 + (1|REC_GROUP_REF/LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                               chains = 3, iter = 3000, warmup = 1000, seed = 23, control = list(max_treedepth = 15, adapt_delta=0.99),
                               save_pars = save_pars(all=T), cores=3, file="B_BEG_prop_ran_04")
 
 loo(B_BEG_prop_ran_00, B_BEG_prop_ran_01, B_BEG_prop_ran_02, B_BEG_prop_ran_03, B_BEG_prop_ran_04, moment_match = F)
 
 ###### ICC ####
 performance::variance_decomposition(B_BEG_prop_ran_00)
 performance::r2_bayes(B_BEG_prop_ran_00) 

 performance::variance_decomposition(B_BEG_prop_ran_01)
 performance::r2_bayes(B_BEG_prop_ran_01) 
 
 performance::variance_decomposition(B_BEG_prop_ran_02)
 performance::r2_bayes(B_BEG_prop_ran_02) 
 
 performance::variance_decomposition(B_BEG_prop_ran_03)
 performance::r2_bayes(B_BEG_prop_ran_03) 

 performance::variance_decomposition(B_BEG_prop_ran_04)
 performance::r2_bayes(B_BEG_prop_ran_04) 
 
 rm(B_BEG_prop_ran_00, B_BEG_prop_ran_01, B_BEG_prop_ran_02, B_BEG_prop_ran_03, B_BEG_prop_ran_04)
 
 #### 2) Define all (biologically plausible) models ####
 B_BEG_prop_00 <- brms::brm(formula = Sum_BEG | trials(Total_calls) ~ +1 + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_00")
 # only one predictor
 B_BEG_prop_01 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_01")
 # additive predictors - no interactions
 B_BEG_prop_02 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_02")
 B_BEG_prop_03 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + SEX + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_03")
 B_BEG_prop_04 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_04")
 B_BEG_prop_04_a <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + SEX + WEIGHT_DIFF_PER + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_04_a")
 B_BEG_prop_05 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT +  REC_AGE_D + SEX + Helper_Pup_REC + WEIGHT_DIFF_PER+  (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_05")
 # interactions
 B_BEG_prop_06 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT *  REC_AGE_D + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, init_r=0.1, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_06")
 B_BEG_prop_07 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT *  REC_AGE_D * SEX + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, init_r=0.1 ,control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_07")
 B_BEG_prop_08 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, init_r=0.11, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_08")
 B_BEG_prop_08_a <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT *  REC_AGE_D * SEX * WEIGHT_DIFF_PER + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, init_r=0.005, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_08_a")
 B_BEG_prop_09 <- brms::brm(formula =Sum_BEG | trials(Total_calls) ~ TREATMENT *  REC_AGE_D * SEX * Helper_Pup_REC * WEIGHT_DIFF_PER + (1|LITTER_CODE/ID), data = BEG_prop_data, family = zero_inflated_beta_binomial(link='logit'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, init_r=0.005, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, file = "B_BEG_prop_09")
 
 #### CHECK MODELS ####
 summary(B_BEG_prop_00) 
 #posterior_summary(B_BEG_prop_00)
 plot(B_BEG_prop_00) 
 pp_check(B_BEG_prop_00, ndraws = 100) 
 model_loo <- loo(B_BEG_prop_00, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_01) 
 #posterior_summary(B_BEG_prop_01)
 plot(B_BEG_prop_01) 
 pp_check(B_BEG_prop_01, ndraws = 100) 
 model_loo <- loo(B_BEG_prop_01, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_02) 
 #posterior_summary(B_BEG_prop_02)
 plot(B_BEG_prop_02) 
 pp_check(B_BEG_prop_02, ndraws = 100) 
 model_loo <- loo(B_BEG_prop_02, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_03) 
 #posterior_summary(B_BEG_prop_03)
 plot(B_BEG_prop_03) 
 pp_check(B_BEG_prop_03, ndraws = 100) 
 model_loo <- loo(B_BEG_prop_03, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_04) 
 #posterior_summary(B_BEG_prop_04)
 plot(B_BEG_prop_04) 
 pp_check(B_BEG_prop_04, ndraws = 100) 
 model_loo <- loo(B_BEG_prop_04, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_04_a) 
 #posterior_summary(B_BEG_prop_04_a)
 plot(B_BEG_prop_04_a) 
 pp_check(B_BEG_prop_04_a, ndraws = 100) 
 model_loo <- loo(B_BEG_prop_04_a, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_05) 
 #posterior_summary(B_BEG_prop_05)
 plot(B_BEG_prop_05) 
 pp_check(B_BEG_prop_05, ndraws = 100) 
 model_loo <- loo(B_BEG_prop_05, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_06) 
 #posterior_summary(B_BEG_prop_06)
 plot(B_BEG_prop_06) 
 pp_check(B_BEG_prop_06, ndraws = 100) 
 model_loo <- loo(B_BEG_prop_06, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_07) 
 #posterior_summary(B_BEG_prop_07)
 plot(B_BEG_prop_07) 
 pp_check(B_BEG_prop_07, ndraws = 100) 
 model_loo <- loo(B_BEG_prop_07, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_08) 
 #posterior_summary(B_BEG_prop_08)
 plot(B_BEG_prop_08) 
 pp_check(B_BEG_prop_08, ndraws = 100) #
 model_loo <- loo(B_BEG_prop_08, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_08_a) 
 #posterior_summary(B_BEG_prop_08_a)
 plot(B_BEG_prop_08_a) 
 pp_check(B_BEG_prop_08_a, ndraws = 100) #
 model_loo <- loo(B_BEG_prop_08_a, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 summary(B_BEG_prop_09) 
 #posterior_summary(B_BEG_prop_09)
 plot(B_BEG_prop_09) 
 pp_check(B_BEG_prop_09, ndraws = 100) 
 model_loo <- loo(B_BEG_prop_09, save_psis = TRUE, cores=4)
 k_rintercept <- model_loo$diagnostics$pareto_k
 df <- tibble(obs_idx = seq_along(k_rintercept), 
              khat = k_rintercept)
 ggplot(df, aes(x = obs_idx, y = khat)) +
   geom_point(alpha=0.5) +
   geom_hline(aes(yintercept=0)) +
   ylim(-1,1)
 
 rm(model_loo, k_rintercept, df)
 
 #### MODEL COMPARISON ####
 # leave-one-out-cross-validation (Vehtari et al, 2017)
 loo(B_BEG_prop_00, B_BEG_prop_01, B_BEG_prop_02, B_BEG_prop_03, B_BEG_prop_04,
     B_BEG_prop_04_a, B_BEG_prop_05, B_BEG_prop_06, B_BEG_prop_07, B_BEG_prop_08, 
     B_BEG_prop_08_a, B_BEG_prop_09)
 
 #### MOMENT MATCHING ####
 loo(B_BEG_prop_00, B_BEG_prop_01, B_BEG_prop_02, B_BEG_prop_03, B_BEG_prop_04,
     B_BEG_prop_04_a, B_BEG_prop_05, B_BEG_prop_06, B_BEG_prop_07, B_BEG_prop_08, 
     B_BEG_prop_08_a, B_BEG_prop_09, moment_match = T)
 
 #### BEST MODEL ####
 summary(B_BEG_prop_09)
 
 exp(fixef(B_BEG_prop_09)) # ODDS RATIO

 mcmc_plot(B_BEG_prop_09, type = 'intervals', prob = 0.89, prob_outer=0.95, )
 mcmc_plot(B_BEG_prop_09, type = 'acf')
 
 pd <- p_direction(B_BEG_prop_09)
 plot(pd)
 
 diagnostic_posterior(B_BEG_prop_09)
 (cred_int <- ci(B_BEG_prop_09, ci=c(.89, .95), effects=c('fixed')))
 
 # ROPE: different definition for logistic regression: -0.18, 0.18
 rope_rang <- c((sd(BEG_prop_data$Avg_BEG_prop)*0.18)*-1, sd(BEG_prop_data$Avg_BEG_prop)*0.18)
 
 describe_posterior(
   B_BEG_prop_09,
   effects = "fixed", #fixed vs all (for random effects)
   component = "all",
   rope_range = rope_rang,
   test = c("p_direction", "p_significance", "rope"),
   centrality = "all",
   dispersion = TRUE
 )

 equivalence_test(B_BEG_prop_09, range = rope_rang)
 
 post_best <- describe_posterior(
   B_BEG_prop_09,
   effects = "fixed",
   component = "all",
   test = c("p_direction", "p_significance"),
   centrality = "all",
   dispersion = TRUE
 )
# convert to odds!
 post_best$Median <- exp(post_best$Median)
 post_best$CI_high <- exp(post_best$CI_high)
 post_best$CI_low <- exp(post_best$CI_low)
 post_best <- post_best[-c(4,5)]
 post_best
 
 loo_R2(B_BEG_prop_09)
 performance::variance_decomposition(B_BEG_prop_09)
 
 ## CALCULATE ESTIMATES FOR INTERACTIONS #### 
 # TREATMENT*HELPER*WEIGHT
 weight_vars <- c(mean(BEG_prop_data$WEIGHT_DIFF_PER)-sd(BEG_prop_data$WEIGHT_DIFF_PER), 
                  mean(BEG_prop_data$WEIGHT_DIFF_PER), 
                  mean(BEG_prop_data$WEIGHT_DIFF_PER) + sd(BEG_prop_data$WEIGHT_DIFF_PER))
 
 (treat_helper_weight <- emtrends(B_BEG_prop_09, specs = c('TREATMENT', 'WEIGHT_DIFF_PER'), var = "Helper_Pup_REC",
                                              at = list(WEIGHT_DIFF_PER = weight_vars)))

 
 p_direction(treat_helper_weight)

 
 # SEX*HELPER*WEIGHT (using 3 categories for continuous vars)
 (sex_helper_weight <- emtrends(B_BEG_prop_09, specs = c('SEX', 'WEIGHT_DIFF_PER'), var = "Helper_Pup_REC",
                                  at = list(WEIGHT_DIFF_PER = weight_vars)))

 p_direction(sex_helper_weight)
 
  ### PLOTS: BEST MODEL ####
 # plot all predictors
 conditional_effects(B_BEG_prop_09, re_formula = NULL, robust=TRUE)
 
 # min(BEG_prop_data$REC_AGE_D)#31
 # max(BEG_prop_data$REC_AGE_D)#130
 
 #AGE
 age_dist <- B_BEG_prop_09 %>% 
   epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_prop_data$TREATMENT),
                                     SEX = levels(BEG_prop_data$SEX),
                                     Helper_Pup_REC = mean(BEG_prop_data$Helper_Pup_REC),
                                     WEIGHT_DIFF_PER = mean(BEG_prop_data$WEIGHT_DIFF_PER),
                                     REC_AGE_D = seq(30, 130, by=5),
                                     Total_calls = 1), 
               re_formula = NULL, allow_new_levels=T)
 
 age_dist$Avg_BEG_prop <- age_dist$.epred
 #age only
 ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_BEG_prop)) +
   stat_lineribbon(.width = c(.95))+
   geom_point(data = BEG_prop_data, size = 2, aes(color=TREATMENT)) +   # raw data
   scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
   scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
   labs(x = "Age (days)", y = "Begging call proportion") +
   scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
   guides(fill='none')+
   theme_clean()
 # age & treatment
 ggplot(age_dist, aes(x = REC_AGE_D, y = Avg_BEG_prop, color=TREATMENT, fill=TREATMENT)) +  
   stat_lineribbon(.width=.95)+
   geom_point(data = BEG_prop_data, size = 2, aes(color=TREATMENT)) +   # raw data
   scale_color_okabe_ito(name = "Treatment", labels=c('Control', 'Flutamide'))+
   scale_fill_okabe_ito(alpha=0.2, name = "Treatment", labels=c('Control', 'Flutamide'))+
   labs(x = "Age (days)", y = "Begging call proportion") +
   scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130))+
   theme_clean()
 
 #SEX * HELPER * WEIGHT
 weight_vars <- c(mean(BEG_prop_data$WEIGHT_DIFF_PER)-sd(BEG_prop_data$WEIGHT_DIFF_PER), 
                  mean(BEG_prop_data$WEIGHT_DIFF_PER), 
                  mean(BEG_prop_data$WEIGHT_DIFF_PER) + sd(BEG_prop_data$WEIGHT_DIFF_PER))
 
 sex_helper_weight_dist <- B_BEG_prop_09 %>% 
   epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_prop_data$TREATMENT),
                                     SEX = levels(BEG_prop_data$SEX),
                                     Helper_Pup_REC = seq(floor(min(BEG_prop_data$Helper_Pup_REC)), 
                                                          ceiling(max(BEG_prop_data$Helper_Pup_REC)), by=1),
                                     WEIGHT_DIFF_PER = weight_vars,
                                     REC_AGE_D = mean(BEG_prop_data$REC_AGE_D),
                                     Total_calls = 1), 
               re_formula = NULL, allow_new_levels=T)
 
 sex_helper_weight_dist$Avg_BEG_prop <- sex_helper_weight_dist$.epred
 sex_helper_weight_dist$WEIGHT_DIFF_PER <- as.factor(sex_helper_weight_dist$WEIGHT_DIFF_PER)
 
 sex.labs <- c("Female", "Male")
 names(sex.labs) <- c("F", "M")

 ggplot(sex_helper_weight_dist, aes(x = Helper_Pup_REC, y = Avg_BEG_prop, color=factor(WEIGHT_DIFF_PER), fill=factor(WEIGHT_DIFF_PER))) +
   stat_lineribbon(.width = c(.95))+
   geom_point(data = BEG_prop_data, size = 2, color='black', fill='black', alpha=0.7, shape='+') +
   scale_color_okabe_ito(name = "Weight condition", labels=c("Poor","Normal","Good"))+
   scale_fill_okabe_ito(alpha=0.2,name = "Weight condition", labels=c("Poor","Normal","Good"))+
   labs(x = "Helper/Pup ratio", y = "Begging call proportion") +
   guides(fill='none')+
   theme_clean()+
   facet_wrap(~SEX, labeller = labeller(SEX = sex.labs))

# TREATMENT * HELPER * WEIGHT
 weight_vars <- c(mean(BEG_prop_data$WEIGHT_DIFF_PER)-sd(BEG_prop_data$WEIGHT_DIFF_PER), 
                  mean(BEG_prop_data$WEIGHT_DIFF_PER), 
                  mean(BEG_prop_data$WEIGHT_DIFF_PER) + sd(BEG_prop_data$WEIGHT_DIFF_PER))
 
 treatment_helper_weight_dist <- B_BEG_prop_09 %>% 
   epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_prop_data$TREATMENT),
                                     SEX = levels(BEG_prop_data$SEX),
                                     Helper_Pup_REC = seq(floor(min(BEG_prop_data$Helper_Pup_REC)), 
                                                          ceiling(max(BEG_prop_data$Helper_Pup_REC)), by=1),
                                     WEIGHT_DIFF_PER = weight_vars,
                                     REC_AGE_D = mean(BEG_prop_data$REC_AGE_D),
                                     Total_calls = 1), 
               re_formula = NULL, allow_new_levels=T)
 
 treatment_helper_weight_dist$Avg_BEG_prop <- treatment_helper_weight_dist$.epred
 treatment_helper_weight_dist$WEIGHT_DIFF_PER <- as.factor(treatment_helper_weight_dist$WEIGHT_DIFF_PER)
 
 treatment.labs <- c("Control", "Flutamide")
 names(treatment.labs) <- c("CTRL", "FLUT")
 
 ggplot(treatment_helper_weight_dist, aes(x = Helper_Pup_REC, y = Avg_BEG_prop, color=factor(WEIGHT_DIFF_PER), fill=factor(WEIGHT_DIFF_PER))) +
   stat_lineribbon(.width = c(.95))+
   geom_point(data = BEG_prop_data, size = 2, color='black', fill='black', alpha=0.7, shape='+') +
   scale_color_okabe_ito(name = "Weight condition", labels=c("Poor","Normal","Good"))+
   scale_fill_okabe_ito(alpha=0.2,name = "Weight condition", labels=c("Poor","Normal","Good"))+
   labs(x = "Helper/Pup ratio", y = "Begging call proportion") +
   guides(fill='none')+
   theme_clean()+
   facet_wrap(~TREATMENT, labeller = labeller(TREATMENT = treatment.labs)) 
 
 ##CI crossing 0:
 # Treatment * sex
 treatment_dist <- B_BEG_prop_09 %>% 
   epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_prop_data$TREATMENT),
                                     SEX = levels(BEG_prop_data$SEX),
                                     Helper_Pup_REC = mean(BEG_prop_data$Helper_Pup_REC),
                                     WEIGHT_DIFF_PER = mean(BEG_prop_data$WEIGHT_DIFF_PER),
                                     REC_AGE_D = mean(BEG_prop_data$REC_AGE_D), 
                                     Total_calls = 1), 
               re_formula = NULL, allow_new_levels=T)
 
 treatment_dist$Avg_BEG_prop <- treatment_dist$.epred
 treatment_dist$TREATMENT <- as.factor(treatment_dist$TREATMENT)
 treatment_dist$SEX <- as.factor(treatment_dist$SEX)
 
 ggplot(treatment_dist, aes(x = TREATMENT, y = Avg_BEG_prop, color= SEX), fill=SEX) +
   stat_boxplot(geom='errorbar', width=.8)+
   geom_boxplot(outlier.shape = NA)+
   scale_color_okabe_ito(name = "Sex", labels=c("Female", "Male"))+
   scale_fill_okabe_ito(alpha=0.3)+
   #geom_point(data = BEG_len_data, size = 2, aes(color=TREATMENT:SEX), alpha=0.5) +   # raw data
   labs(x = "Treatment", y = "Begging call proportion") +
   scale_x_discrete(labels=c('Control', 'Flutamide'))+
   theme_clean()
 
 # AVERAGE BEST MODELS (BAYESIAN MODEL AVERAGE - BMA) ####
 # # no new data, just get the information of fitted values
 post_avg <- posterior_average(B_BEG_prop_09, B_BEG_prop_08_a, B_BEG_prop_06, weights = "stacking",
                               missing = 0)

 posterior_summary(post_avg)

 #ci(post_avg, ci=c(.89, .95))

 # ROPE: different definition for logistic regression: -0.18, 0.18
 rope_rang <- c((sd(BEG_prop_data$Avg_BEG_prop)*0.18)*-1, sd(BEG_prop_data$Avg_BEG_prop)*0.18)

 describe_posterior(
   post_avg,
   effects = "fixed",
   component = "all",
   rope_range = rope_rang,
   test = c("p_direction", "p_significance", "rope"),
   centrality = "all",
   dispersion = TRUE
 )

 rm(post_avg)

 # cleanup ####
 rm(B_BEG_prop_00, B_BEG_prop_01, B_BEG_prop_02, B_BEG_prop_03, B_BEG_prop_04,
    B_BEG_prop_04_a, B_BEG_prop_05, B_BEG_prop_06, B_BEG_prop_07, B_BEG_prop_08, 
    B_BEG_prop_08_a, B_BEG_prop_09, BEG_prop_data, BEG_data, full_data)