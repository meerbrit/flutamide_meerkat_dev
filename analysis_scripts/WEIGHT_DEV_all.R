###### Pup weight development ##################################################
################################################################################
# model weight development for pups over the first 370 days (long term data)
# to categorise condition = poor, normal, good
# or have a definitive number (point estimate), 
# or have the % off the expected weight
#### BWalkenhorst, 2022 ########################################################

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() 

# Load initial packages
library(readxl)
library(ggplot2)
library(writexl)
library(brms)
library(tidybayes)
library(tidyselect)
library(ggokabeito) 
library(ggmcmc)
library(lisa)
library(performance)

theme_clean <- function() {
  theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_text(face = "bold"))
}

# define the color palette
fk <- lisa_palette("FridaKahlo", n = 31, type = "continuous")

# custom trace plot
geom_trace <- function(subtitle = NULL, 
                       xlab = "iteration", 
                       xbreaks = 0:4 * 500) {
  list(
    annotate(geom = "rect", 
             xmin = 0, xmax = 1000, ymin = -Inf, ymax = Inf,
             fill = fk[16], alpha = 1/2, size = 0),
    geom_line(size = 1/3),
    scale_color_manual(values = fk[c(3, 8, 27, 31)]),
    scale_x_continuous(xlab, breaks = xbreaks, expand = c(0, 0)),
    labs(subtitle = subtitle),
    theme(panel.grid = element_blank())
  )
}

# set working directory
setwd("/Flutamide/")

# set seed to duplicate results
set.seed(42)

# load data
full_data <- read_excel(
  "sheets/ALL_Pup_weights.xlsx"
)
flut_data <- read_excel(
  "sheets/Flutamide_2ndGen_Metadata_BW.xlsx", 
  sheet = 'Flutamide_2ndGen_Metadata'
)

#remove data that does not make sense
full_data <- subset(full_data, AGE_D >= 0 & SEX != 'U')
# remove everything above 6 months
full_data <- subset(full_data, AGE_D <= 180)
# only use first 5 months to reduce dataset (runs forever otherwise!)
full_data <- subset(full_data, AGE_D <= 150)

# remove individuals from flutamide study
full_data <- subset(full_data, !(ID %in% flut_data$ID))
full_data <- na.omit(full_data)
full_data$ID <- as.factor(full_data$ID) 
full_data$SEX <- as.factor(full_data$SEX) 

table(full_data$SEX) 
# F      M  full
# 323429 376850 
# F      M  6 months
# 162768 191007
# F      M  5 months
# 133730 156695 

# unique individuals
sample_df <- full_data[!duplicated(full_data$ID),]
# 2714 individuals
table(sample_df$SEX)
# F    M 
# 1244 1470

# reset working directory
setwd("/scripts/WEIGHT_models/")

hist(full_data$AGE_D)

full_data %>% 
  ggplot(aes(x = Weight)) +
  geom_density(fill = fk[3], color = fk[3])

full_data %>% 
  ggplot(aes(x = AGE_D)) +
  geom_density(fill = fk[3], color = fk[3])

#plot the data
ggplot(full_data,
       aes(AGE_D, Weight, color = SEX))+
  scale_colour_manual(values=c("#ccf381", "#0d5b06", "#68bfdf", "#4831d4"))+
 # geom_jitter(data=full_data, aes(AGE_D, Weight),
#              width=0.1, height=0.03, alpha=0.5, size=1)+
  geom_smooth(se=TRUE, size=1.3, fill="gray87", method = "lm")+ #se True should already add CI
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_blank())+
  labs(x = "Age in days", y = "Weight")+
  theme(text = element_text(size=10))

ggplot(full_data, aes(x = AGE_D, y = Weight, color= SEX)) +  
  stat_lineribbon(.width = c(.95, .89, .50),  # regression line and CI
                  alpha = 0.3) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_okabe_ito()+
  labs(x = "Age in days", y = "Weight",
       fill = "CI: 50%, 89%, 95%") +
  scale_x_continuous(breaks=c(15, 30, 50, 60, 70, 90, 100, 120, 150))+
  theme_clean()

cor(sapply(full_data[, c( "Weight", "SEX", "AGE_D")], as.numeric))
#         Weight         SEX       AGE_D
# Weight 1.000000000 0.003932623 0.836194217
# SEX    0.003932623 1.000000000 0.003278847
# AGE_D  0.836194217 0.003278847 1.000000000

##############################################################################################################################
# initial values for better convergence
inits <- list(
  Intercept = 89,
  SEXM = -1
)
list_of_inits <- list(inits, inits, inits, inits)

#### MODELS

# these models can be found in Flutamide study/scripts/WEIGHT_models
# age_mod_00 <- brms::brm(formula =Weight ~ 1 + (1|ID), data = full_data, family = gaussian(link='identity'),
#                         chains = 4, iter = 20000, warmup = 5000, cores = 6, backend = "cmdstanr", threads = threading(3),
#                         init=list_of_inits, thin=10, file="dev_weight_00")
# 
# age_mod_01 <- brms::brm(formula =Weight ~ AGE_D + (1|ID), data = full_data, family = gaussian(link='identity'),
#                           chains = 4, iter = 20000, warmup = 5000, cores = 6, threads = threading(3), backend = "cmdstanr",
#                         init=list_of_inits, thin=10, file ="dev_weight_01")
# summary(age_mod_01)
# mcmc_plot(age_mod_01, type = 'acf') # high autocorrelation, thinning needed but ok

#this is the main one
age_mod_02 <- brms::brm(formula =Weight ~ AGE_D + SEX + (1|ID), data = full_data, family = gaussian(link='identity'),
                        chains = 4, iter = 20000, warmup = 5000, cores = 6, backend = "cmdstanr", threads = threading(3),
                        init=list_of_inits, thin=10, file="dev_weight_02")

###############################################################

# Model checks: ####

summary(age_mod_00) 
#posterior_summary(age_mod_00)
plot(age_mod_00) #
pp_check(age_mod_00, ndraws = 100) #
# Pareto-k diagnostic can be useful to identify problematic point(s)
model_loo <- brms::loo(age_mod_00, save_psis = TRUE, cores=2)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(age_mod_00, type = 'intervals')
mcmc_plot(age_mod_00, type = 'acf') 

summary(age_mod_01) 
#posterior_summary(age_mod_01)
plot(age_mod_01) 
pp_check(age_mod_01, ndraws = 100) 
model_loo <- brms::loo(age_mod_01, save_psis = TRUE, cores=2)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)
mcmc_plot(age_mod_01, type = 'intervals')
mcmc_plot(age_mod_01, type = 'acf')

summary(age_mod_02) 
#posterior_summary(age_mod_02)
plot(age_mod_02) 
pp_check(age_mod_02, ndraws = 10) 
# Pareto-k diagnostic can be useful to identify problematic point(s)
model_loo <- brms::loo(age_mod_02, save_psis = TRUE, cores=4)
k_rintercept <- model_loo$diagnostics$pareto_k
df <- tibble(obs_idx = seq_along(k_rintercept), 
             khat = k_rintercept)
ggplot(df, aes(x = obs_idx, y = khat)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1,1)

mcmc_plot(age_mod_02, type = 'intervals')
mcmc_plot(age_mod_02, type = 'acf') 
modeltranformed <- ggs(age_mod_02) 
ggplot(filter(modeltranformed, Parameter %in% c("b_Intercept", "b_AGE_D", "b_SEXM", "sigma")),
       aes(x   = Iteration,
           y   = value, 
           col = as.factor(Chain)))+
  geom_line()+
  facet_grid(Parameter ~ .,
             scale     = 'free_y',
             switch    = 'y')+
  labs(title = "Trace plots",
       col   = "Chains") +
  theme_minimal()
rm(modeltranformed)

pd <- p_direction(age_mod_02)
plot(pd)

hdi(age_mod_02)
# Highest Density Interval 
# 
# Parameter   |        95% HDI
# ----------------------------
# (Intercept) | [86.53, 91.80]
# AGE_D       | [ 2.32,  2.33] ###
# SEXM        | [-4.84,  1.99]

describe_posterior(
  age_mod_02,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(age_mod_02),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)
# Summary of Posterior Distribution 
# 
# Parameter   | Median |      MAD |  Mean |       SD |   MAP |         95% CI |     pd |   ps |            ROPE | % in ROPE |  Rhat |     ESS
# -------------------------------------------------------------------------------------------------------------------------------------------
# (Intercept) |  89.22 |     1.36 | 89.23 |     1.33 | 88.97 | [86.61, 91.90] |   100% | 1.00 | [-10.26, 10.26] |        0% | 1.014 |  368.00
# AGE_D       |   2.33 | 1.76e-03 |  2.33 | 1.75e-03 |  2.33 | [ 2.32,  2.33] |   100% | 0.00 | [-10.26, 10.26] |      100% | 1.000 | 5587.00
# SEXM        |  -1.43 |     1.77 | -1.42 |     1.76 | -1.44 | [-4.81,  2.05] | 78.83% | 0.00 | [-10.26, 10.26] |      100% | 1.007 |  341.00
# 
# # Fixed effects sigma
# 
# Parameter | Median |  MAD |  Mean |   SD |   MAP |         95% CI |   pd |   ps |            ROPE | % in ROPE |  Rhat |     ESS
# -------------------------------------------------------------------------------------------------------------------------------
# sigma     |  32.37 | 0.04 | 32.37 | 0.04 | 32.38 | [32.29, 32.46] | 100% | 1.00 | [-10.26, 10.26] |        0% | 1.000 | 6013.00

equivalence_test(age_mod_02)
# # Test for Practical Equivalence
# 
# ROPE: [-10.26 10.26]
# 
# Parameter |       H0 | inside ROPE |       95% HDI
# --------------------------------------------------
# Intercept | Rejected |      0.00 % | [86.61 91.90]
# AGE_D     | Accepted |    100.00 % | [ 2.32  2.33]
# SEXM      | Accepted |    100.00 % | [-4.81  2.05]

posterior_summary(age_mod_02)

performance::r2_loo(age_mod_02) # adjusted R2
performance::variance_decomposition(age_mod_02)

loo(age_mod_00, age_mod_01, age_mod_02)

#### PREDICTIONS ####
### general predictions ####
Weight ~ AGE_D + SEX + (1 | ID) 
prediction_data <- expand.grid(ID = seq(0, 10, 1),
                              SEX = factor(unique(full_data$SEX)),
                              AGE_D = seq(0, 200, 1))

general_weight_pred <- predict(age_mod_02, newdata = prediction_data, seed=23, allow_new_levels=T) 
saveRDS(general_weight_pred, "general_weight_pred.rds")
prediction_data$WEIGHT <- general_weight_pred[,1]

# plot predictions
ggplot(prediction_data, aes(x = AGE_D, y = WEIGHT)) +  
  stat_lineribbon(.width = c(.95, .89, .50),  # regression line and CI
                  alpha=0.3) +
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Age (days)", y = "Predicted weight (grams) \n",
       fill = "CI") +
  scale_x_continuous(breaks=c(0, 15, 30, 50, 60, 70, 90, 100, 120, 150, 180, 200))+
  scale_y_continuous(breaks=seq(0, 650, 50))+
  theme_clean()

ggplot(prediction_data,
       aes(AGE_D, WEIGHT, color=SEX))+
  geom_smooth(se=TRUE, size=1.3, fill="gray87", method = "lm")+ #se True should already add CI
  scale_x_continuous(breaks=c(0, 15, 30, 50, 60, 70, 90, 100, 120, 150, 180, 200))+
  labs(x = "Age (days)", y = "Weight (g)")+
  scale_color_okabe_ito()+
  theme_clean()

###### predictions using study data ####
copy_flut <- flut_data
copy_flut$AGE_D <- copy_flut$REC_AGE_D

# get the predictions from the model
weight_pred <- predict(age_mod_02, newdata = copy_flut, seed=23, allow_new_levels=T) 
saveRDS(weight_pred, "weight_pred.rds")
rm(copy_flut)

weight_pred <- readRDS("weight_pred.rds")

flut_data$WEIGHT_PRED <- weight_pred[,1]
flut_data$WEIGHT_CI.lo <- weight_pred[,3]
flut_data$WEIGHT_CI.hi <- weight_pred[,4]

#plot predictions and actual average weight (jitter)
ggplot(flut_data,
       aes(REC_AGE_D, WEIGHT_PRED, color = SEX))+
  scale_colour_manual(values=c("#ccf381", "#0d5b06", "#68bfdf", "#4831d4"))+
  geom_jitter(data=flut_data, aes(REC_AGE_D, AvgWeight),
              width=0.1, height=0.03, alpha=0.5, size=1)+
  geom_smooth(se=TRUE, size=1.3, fill="gray87", method = "lm")+ #se True should already add CI
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_blank())+
  scale_x_continuous(breaks=c(15, 30, 50, 60, 70, 90, 100, 120, 150))+
  labs(x = "Age in days", y = "Weight (g)")+
  theme(text = element_text(size=10))

ggplot(flut_data, aes(x = REC_AGE_D, y = WEIGHT_PRED, color= SEX)) +  
  stat_lineribbon(.width = c(.95, .89, .50),  # regression line and CI
                  alpha = 0.3,) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_okabe_ito()+
  geom_jitter(data=flut_data, aes(REC_AGE_D, AvgWeight),
              width=0.1, height=0.03, alpha=0.5, size=1)+
  labs(x = "Age in days", y = "Weight",
       fill = "CI: 50%, 89%, 95%") +
  scale_x_continuous(breaks=c(15, 30, 50, 60, 70, 90, 100, 120, 150))+
  theme_clean()


# determine the condition based on credible interval
determineConditionByCI <- function(x){
  condition <- 'NORMAL'
  if(as.numeric(x['AvgWeight']) < as.numeric(x['WEIGHT_CI.lo'])) {
    condition <- 'POOR'
  }else if (as.numeric(x['AvgWeight']) > as.numeric(x['WEIGHT_CI.hi'])) {
    condition <- 'GOOD' 
  } 
  return(condition)
}

# determine condition based on percentage off the predicted (point) weight 
determineConditionByPercent <- function(x){
  condition <- 'NORMAL'
  if(as.numeric(x['WEIGHT_DIFF_PER']) <= (-20)) {
    condition <- 'POOR'
  }else if (as.numeric(x['WEIGHT_DIFF_PER']) >= 20) {
    condition <- 'GOOD' 
  } 
  return(condition)
}

flut_data$CONDITION <- apply(flut_data, 1, determineConditionByCI)
table(flut_data$CONDITION, flut_data$TREATMENT) # all weights NORMAL 

# not taking CI into account and just using difference predicted vs actual??
flut_data$WEIGHT_DIFF <- flut_data$AvgWeight - flut_data$WEIGHT_PRED 

# or use percentage off from point estimate/prediction??
flut_data$WEIGHT_DIFF_PER <-(flut_data$AvgWeight /(flut_data$WEIGHT_PRED/100)) -100

# or use certain percentage to classify condition (0% off?)
flut_data$CONDITION_PER <- apply(flut_data, 1, determineConditionByPercent)
table(flut_data$CONDITION_PER, flut_data$TREATMENT) # all weights normal


# save the new data
filepath_xlsx <- ''
write_xlsx(flut_data, filepath_xlsx)