#### Analysis script for Flutamide study: Overall proportions / DEVELOPMENT #####
###### Bayesian Multilevel models, best model & Bayesian Model Averaging (BMA)##
############## BWAlkenhorst 2022 ###############################################

#### SETUP ####

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(ggplot2) 
library(ggokabeito)
library(tidyverse)
library(tidybayes)

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
BEG_data <- full_data[, c(1:18, 79, 105:113, 118, 120,  122,  125:126, 128, 147, 162, 166:169)]
# remove NAs
BEG_data <- na.omit(BEG_data) 
BEG_prop_data <- subset(BEG_data, nBouts >= MIN_BOUTS) #call proportion
BEG_prop_data$Total_calls <- BEG_prop_data$Sum_BEG+
  BEG_prop_data$Sum_DIG+BEG_prop_data$Sum_CC+BEG_prop_data$Sum_OTHER #save the number of total calls for calc later
BEG_prop_data$CallType <- as.factor('BEG')

DIG_data <- full_data[, c(19:33, 79, 105:113, 118, 120,  122, 125:126, 128, 147, 162, 166:169)]
# remove NAs
DIG_data <- na.omit(DIG_data) 
DIG_prop_data <- subset(DIG_data, nBouts >= MIN_BOUTS) #call proportion
DIG_prop_data$Total_calls <- DIG_prop_data$Sum_BEG+
  DIG_prop_data$Sum_DIG+DIG_prop_data$Sum_CC+DIG_prop_data$Sum_OTHER 
DIG_prop_data$CallType <- as.factor('DIG')

CC_data <- full_data[, c(34:48, 79, 105:113, 118, 120,  122, 125:126, 128, 147, 162, 166:169)]
# remove NAs
CC_data <- na.omit(CC_data) 
CC_prop_data <- subset(CC_data, nBouts >= MIN_BOUTS) #call proportion
CC_prop_data$Total_calls <- CC_prop_data$Sum_BEG+
  CC_prop_data$Sum_DIG+CC_prop_data$Sum_CC+CC_prop_data$Sum_OTHER
CC_prop_data$CallType <- as.factor('CC')

rm(BEG_data, DIG_data, CC_data, full_data, MIN_BOUTS)

# switch working directory to models folder
setwd("G:/My Drive/Uni/UZH/Projects/KMP_PUPS/Flutamide/scripts/")

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

#load best models for each call type ####
BEG_model <- readRDS('BEG_models/B_BEG_prop_09.rds') 
DIG_model <- readRDS('DIG_models/B_DIG_prop_09.rds') 
CC_model <- readRDS('CC_models/B_CC_prop_09.rds') 


#### ESTIMATES ####
#BEG
BEG_dist <- BEG_model %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_prop_data$TREATMENT),
                                    SEX = levels(BEG_prop_data$SEX),
                                    Helper_Pup_REC = mean(BEG_prop_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(BEG_prop_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5),
                                    Total_calls=1), 
              re_formula = NULL, allow_new_levels=T)

BEG_dist$CallType <- 'BEG'
BEG_dist$TREATMENT <- as.factor(BEG_dist$TREATMENT)
BEG_dist$SEX <- as.factor(BEG_dist$SEX)
BEG_dist$Call_prop <- BEG_dist$.epred

#DIG
DIG_dist <- DIG_model %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_prop_data$TREATMENT),
                                    SEX = levels(DIG_prop_data$SEX),
                                    Helper_Pup_REC = mean(DIG_prop_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_prop_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5),
                                    Total_calls=1), 
              re_formula = NULL, allow_new_levels=T)

DIG_dist$CallType <- 'DIG'
DIG_dist$TREATMENT <- as.factor(DIG_dist$TREATMENT)
DIG_dist$SEX <- as.factor(DIG_dist$SEX)
DIG_dist$Call_prop <- DIG_dist$.epred

#CC
CC_dist <- CC_model %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(CC_prop_data$TREATMENT),
                                    SEX = levels(CC_prop_data$SEX),
                                    Helper_Pup_REC = mean(CC_prop_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(CC_prop_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(30, 130, by=5),
                                    Total_calls=1), 
              re_formula = NULL, allow_new_levels=T)

CC_dist$CallType <- 'CC'
CC_dist$TREATMENT <- as.factor(CC_dist$TREATMENT)
CC_dist$SEX <- as.factor(CC_dist$SEX)
CC_dist$Call_prop <- CC_dist$.epred

# create new df combining all results ####
proportion_data <- do.call("rbind", list(BEG_dist, DIG_dist, CC_dist))
proportion_data$CallType <- as.factor(proportion_data$CallType)
proportion_data$TREATMENT <- as.factor(proportion_data$TREATMENT)

#### PLOT ESTIMATES ####

# call type and treatment: two plots
ggplot(proportion_data, aes(x = REC_AGE_D, y = Call_prop, colour=CallType, fill=CallType)) +  
  stat_lineribbon(.width = c(.95) ) +
  scale_color_okabe_ito(name = "Call types", breaks=c('BEG', 'DIG', 'CC'), labels=c('Begging call', 'Digging call','Contact call'))+
  scale_fill_okabe_ito(alpha=0.2, name = "Call types", breaks=c('BEG', 'DIG', 'CC'), labels=c('Begging call', 'Digging call','Contact call'))+
  geom_point(data = BEG_prop_data, size = 1, aes(y = Avg_BEG_prop)) +   # raw data
  geom_point(data = DIG_prop_data, size = 1, aes(y= Avg_DIG_prop)) +   # raw data
  geom_point(data = CC_prop_data, size = 1, aes(y = Avg_CC_prop)) +   # raw data
  labs(x = "Age (days)", y = "Call proportion") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130), guide = guide_axis(angle = 45))+
  theme_clean()+
  facet_wrap(~TREATMENT, labeller=as_labeller(c("CTRL" = "Control", "FLUT" = "Flutamide")))

## call type and treatment: single plot
ggplot(proportion_data, aes(x = REC_AGE_D, y = Call_prop, colour=CallType:TREATMENT, fill=CallType:TREATMENT, linetype=TREATMENT)) +
  stat_lineribbon(.width = c(.95)) +
  scale_color_okabe_ito(name = "Call types and treatment", labels=c('Begging call: Control', 'Begging call: Flutamide', 'Contact call: Control', 'Contact call: Flutamide','Digging call: Control', 'Digging call: Flutamide'))+
  scale_fill_okabe_ito(alpha=0.3, name = "Call types and treatment", labels=c('Begging call: Control', 'Begging call: Flutamide', 'Contact call: Control', 'Contact call: Flutamide','Digging call: Control', 'Digging call: Flutamide'))+
  labs(x = "Age (days)", y = "Call proportion") +
  scale_x_continuous(breaks=c(20, 30, 50, 60, 70, 90, 100, 120, 130), guide = guide_axis(angle = 45))+
  guides(linetype='none')+
  theme_clean()

#### PREDICTIONS (6MO)####
#BEG
BEG_dist_pred <- BEG_model %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(BEG_prop_data$TREATMENT),
                                    SEX = levels(BEG_prop_data$SEX),
                                    Helper_Pup_REC = mean(BEG_prop_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(BEG_prop_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(0, 180, by=5),
                                    Total_calls=1), 
              re_formula = NULL, allow_new_levels=T)

BEG_dist_pred$CallType <- 'BEG'
BEG_dist_pred$TREATMENT <- as.factor(BEG_dist_pred$TREATMENT)
BEG_dist_pred$SEX <- as.factor(BEG_dist_pred$SEX)
BEG_dist_pred$Call_prop <- BEG_dist_pred$.epred

#DIG
DIG_dist_pred <- DIG_model %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_prop_data$TREATMENT),
                                    SEX = levels(DIG_prop_data$SEX),
                                    Helper_Pup_REC = mean(DIG_prop_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(DIG_prop_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(0, 180, by=5),
                                    Total_calls=1), 
              re_formula = NULL, allow_new_levels=T)

DIG_dist_pred$CallType <- 'DIG'
DIG_dist_pred$TREATMENT <- as.factor(DIG_dist_pred$TREATMENT)
DIG_dist_pred$SEX <- as.factor(DIG_dist_pred$SEX)
DIG_dist_pred$Call_prop <- DIG_dist_pred$.epred

#CC
CC_dist_pred <- CC_model %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(CC_prop_data$TREATMENT),
                                    SEX = levels(CC_prop_data$SEX),
                                    Helper_Pup_REC = mean(CC_prop_data$Helper_Pup_REC),
                                    WEIGHT_DIFF_PER = mean(CC_prop_data$WEIGHT_DIFF_PER),
                                    REC_AGE_D = seq(0, 180, by=5),
                                    Total_calls=1), 
              re_formula = NULL, allow_new_levels=T)

CC_dist_pred$CallType <- 'CC'
CC_dist_pred$TREATMENT <- as.factor(CC_dist_pred$TREATMENT)
CC_dist_pred$SEX <- as.factor(CC_dist_pred$SEX)
CC_dist_pred$Call_prop <- CC_dist_pred$.epred

# create new df combining all results of predictions ####
proportion_data_pred <- do.call("rbind", list(BEG_dist_pred, DIG_dist_pred, CC_dist_pred))
proportion_data_pred$CallType <- as.factor(proportion_data_pred$CallType)
proportion_data_pred$TREATMENT <- as.factor(proportion_data_pred$TREATMENT)

#### PLOT PREDICTIONS ####

# call type and treatment: two plots
ggplot(proportion_data_pred, aes(x = REC_AGE_D, y = Call_prop, color=CallType, fill=CallType)) +  
  stat_lineribbon(.width = c(.95) ) +
  geom_point(data = BEG_prop_data, size = 1, aes(y = Avg_BEG_prop)) +   # raw data
  geom_point(data = DIG_prop_data, size = 1, aes(y= Avg_DIG_prop)) +   # raw data
  geom_point(data = CC_prop_data, size = 1, aes(y = Avg_CC_prop)) +   # raw data
 scale_color_okabe_ito(name = "Call types", breaks=c('BEG', 'DIG', 'CC'), labels=c('Begging call', 'Digging call','Contact call'))+
 scale_fill_okabe_ito(alpha=0.2, name = "Call types", breaks=c('BEG', 'DIG', 'CC'), labels=c('Begging call', 'Digging call','Contact call'))+
  labs(x = "Age (days)", y = "Predicted call proportion") +
  scale_x_continuous(breaks=c(0, 20, 30, 50, 60, 70, 90, 100, 120, 130, 150, 180, 200), guide = guide_axis(angle = 45))+
  theme_clean()+
  facet_wrap(~TREATMENT, labeller=as_labeller(c("CTRL" = "Control", "FLUT" = "Flutamide")))

## call type and treatment: single plot
ggplot(proportion_data_pred, aes(x = REC_AGE_D, y = Call_prop, colour=CallType:TREATMENT, fill=CallType:TREATMENT, linetype=TREATMENT)) +
  stat_lineribbon(.width = c(.95)) +
  scale_color_okabe_ito(name = "Call types and treatment", labels=c('Begging call: Control', 'Begging call: Flutamide', 'Contact call: Control', 'Contact call: Flutamide','Digging call: Control', 'Digging call: Flutamide'))+
  scale_fill_okabe_ito(alpha=0.3, name = "Call types and treatment", labels=c('Begging call: Control', 'Begging call: Flutamide', 'Contact call: Control', 'Contact call: Flutamide','Digging call: Control', 'Digging call: Flutamide'))+
  labs(x = "Age (days)", y = "Predicted call proportion") +
  scale_x_continuous(breaks=c(0, 20, 30, 50, 60, 70, 90, 100, 120, 130, 150, 180, 200), guide = guide_axis(angle = 45))+
  guides(linetype='none')+
  theme_clean()

#cleanup ####
rm(BEG_model, DIG_model, CC_model, BEG_dist, BEG_prop_data,
   DIG_dist, DIG_prop_data, CC_dist, CC_prop_data,
   proportion_data, BEG_dist_pred, DIG_dist_pred, CC_dist_pred,
   proportion_data_pred)




