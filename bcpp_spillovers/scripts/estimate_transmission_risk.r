# Author:    Lerato E. Magosi
# R version: 4.1.2 (2021-11-01)
# Platform:  x86_64-apple-darwin17.0 (64-bit)
# Date:      03092023

# ----------------------------------------------------------------------------------------



# Goal: Fit a negative binomial regression model to counts of HIV transmission events identified
#       between bcpp communities at baseline versus after baseline with: same community transmission, control community 
#       as the source of transmission and drive distance as predictors

# Background: nb.glm occasionally fails early when datasets have lots of zeros hence
#             previously fit the negative binomial models in stata to obtain alpha 
#             values that we could then use as the starting values for theta. 
#             Reminder! theta =  1 / alpha. 
#
#             For reference see stata model: estimate_transmission_risk.do                                  
#                                                        
 
# Note 1. Because there were no post baseline samples in "Digawana" community we shall exclude 
#     "Diawana" as a recipient community in the post baseline model


# Note 2. Acount for the sex ratio of males to females in the number of sampled participants. 
#
#     To account for the varying numbers of sampled males versus females in bcpp communities
#     and the heterogeneous prevalence of hiv by gender we fit the model using the distinct 
#     possible opposite-sex trm pairs in the sample (mf, fm) as the offset instead of all distinct 
#     possible trm pairs in the sample (mm, ff, mf, fm).


# Platform: local macOS



# Required libraries
# ---------------------

library(dplyr)         # for sorting and merging data.frames
library(tidyr)         # for combining or splitting fields
library(MASS)          # for predicting fitted values and fitting negative binomial models
library(psych)         # for descriptive statistics
library(arm)           # for extracting standard errors of regression coefficients



# Load functions 
# ---------------------



# Load datasets
# ---------------------

# set working dir
setwd("./bcpp_spillovers")
getwd()

# Load dataset of pairwise drive distances between bcpp communities
# Note: This dataset flags likely baseline and post-baseline transmission events
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline <- read.csv("data/sub_m_bcpp_pairwise_dist_same_comm_oc_gender_combined.csv")

# Let's view the data structure with glimpse() and str()
dplyr::glimpse(m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline)

utils::str(m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline)

# Community pairs with identified transmission events
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::count(num_linked_pairs_observed, name = "community_pairings")
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::count(num_linked_pairs_observed > 0, name = "community_pairings")

# Community pairs with identified transmission events at baseline
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::count(num_linked_pairs_observed_baseline, name = "community_pairings")
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::count(num_linked_pairs_observed_baseline > 0, name = "community_pairings")

# Community pairs with identified transmission events post baseline
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::count(num_linked_pairs_observed_post_baseline, name = "community_pairings")
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::count(num_linked_pairs_observed_post_baseline > 0, name = "community_pairings")

# ***

# Load dataset of pairwise drive distances between communities in the 2011 bw census
m_bw_2011_census_pairwise_dist_pop <- read.csv("data/sub_m_bw_2011_census_pairwise_dist_pop.csv")
dplyr::glimpse(m_bw_2011_census_pairwise_dist_pop)

# ***

# Load summary of number of participants sampled at baseline versus post baseline by gender
counts_bcpp_number_sampled_by_community_combined <- read.csv("data/counts_bcpp_number_sampled_by_community_gender_combined.csv")
counts_bcpp_number_sampled_by_community_combined



# Fit models
# ---------------------

# Set alpha values
# Note: alpha values obtained from stata versions of models: estimate_transmission_risk.do

# alpha values for linear models

# baseline model
alpha_m_nb_oc_same_comm_linear_baseline <- .3406354964777654

# baseline model that excludes "Digawana" as a recipient communty
alpha_m_nb_oc_same_comm_linear_baseline_2 <- .31495946297113

# post_baseline model
alpha_m_nb_oc_same_comm_linear_post_baseline <- .2297758067049759

# ***


# Linear models

# Baseline

# Fit model of num_linked_pairs_observed_baseline with 
#     origin_control_community, same community transmission and pairwise drive distance as covariates

# Note: Offset for model with baseline transmission pairs is the natural log of the product of the 
#       number sampled at baseline for communities i and j, where i is the source (origin) community 
#       and j the recipient (destination) community


# Accounting for the sex ratio, the offset for the baseline model becomes:
# 
#       number of distinct possible opposite-sex (male-female, female-male) hiv transmission pairs in the sample i.e.
#       (number males sampled at baseline for community i * number females sampled at baseline for community j) + (number females sampled at baseline for community i * number males sampled at baseline for community j)


nb_glm_offset_ln_sample_oc_same_comm_linear_baseline <- MASS::glm.nb(num_linked_pairs_observed_baseline ~ origin_control_community + same_community + curr_travel_dist_km + offset(ln_max_possible_opposite_sex_pairs_in_sample_baseline), data = m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline, link = log, init.theta = 1/alpha_m_nb_oc_same_comm_linear_baseline)

# Model summary
summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline)

summary_nb_glm_offset_ln_sample_oc_same_comm_linear_baseline <- data.frame(coef = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline)$coefficients[, "Estimate"], se = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline)$coefficients[, "Std. Error"]) %>% 
    dplyr::mutate(p_value = 2 * stats::pnorm(-abs(coef / se))) %>%     
    dplyr::mutate(lower = coef - (2 * se), upper = coef + (2 * se)) %>% 
    dplyr::mutate_at(.vars = c("coef", "se", "lower", "upper"), exp) %>% 
    dplyr::select(coef, se, lower, upper, p_value)
 
summary_nb_glm_offset_ln_sample_oc_same_comm_linear_baseline
    
# Sanity check!
summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline)

setequal(summary_nb_glm_offset_ln_sample_oc_same_comm_linear_baseline$p_value, summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline)$coefficients[, "Pr(>|z|)"])

# Combined
model_oc_baseline <- cbind(variable = c("Intercept", "origin_control_community", "same_community", "drive_distance_kilometers"),
      summary_nb_glm_offset_ln_sample_oc_same_comm_linear_baseline) 

model_oc_baseline

# ***


# Exclude Digawana as a recipient community for consistency with the post baseline model:

# Why? There were no post-baseline samples in Digawana
# Subset community pairs with maximum distinct possible pairs is zero
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::filter(max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb == 0) %>% dplyr::select(curr_orig_dest)

nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2 <- MASS::glm.nb(num_linked_pairs_observed_baseline ~ origin_control_community + same_community + curr_travel_dist_km + offset(ln_max_possible_opposite_sex_pairs_in_sample_baseline), data = m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::filter(!(dest == "Digawana")), link = log, init.theta = 1/alpha_m_nb_oc_same_comm_linear_baseline_2)

# Model summary
summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2)

summary_nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2 <- data.frame(coef = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2)$coefficients[, "Estimate"], se = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2)$coefficients[, "Std. Error"]) %>% 
    dplyr::mutate(p_value = 2 * stats::pnorm(-abs(coef / se))) %>%     
    dplyr::mutate(lower = coef - (2 * se), upper = coef + (2 * se)) %>% 
    dplyr::mutate_at(.vars = c("coef", "se", "lower", "upper"), exp) %>% 
    dplyr::select(coef, se, lower, upper, p_value)
 
summary_nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2
    
# Sanity check!
summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2)

setequal(summary_nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2$p_value, summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2)$coefficients[, "Pr(>|z|)"])

# Combined
model_oc_baseline_2 <- cbind(variable = c("Intercept", "origin_control_community", "same_community", "drive_distance_kilometers"),
      summary_nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2) 

model_oc_baseline_2

# ***


# Post baseline

# Fit model of num_linked_pairs_observed_post_baseline with 
#     origin_control_community, same community transmission and pairwise drive distance as covariates

# Note: Offset for model with post baseline transmission pairs is the natural log of the product of the 
#       total number sampled in community i and post baseline samples in community j where,
#       i is the source / origin community and j the recipient / destination community

# Accounting for the sex ratio, the offset for the post baseline model becomes:
# 
# number distinct possible opposite-sex (male-female, female-male) hiv transmission pairs in sample i.e.
# (total number males sampled for community i * number females sampled post baseline for community j) + (total number females sampled for community i * number males sampled post baseline for community j)

nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline <- MASS::glm.nb(num_linked_pairs_observed_post_baseline ~ origin_control_community + same_community + curr_travel_dist_km + offset(ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb), data = m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::filter(!(dest == "Digawana")), link = log, init.theta = 1/alpha_m_nb_oc_same_comm_linear_post_baseline)

# Model summary
summary(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)

summary_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline <- data.frame(coef = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)$coefficients[, "Estimate"], se = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)$coefficients[, "Std. Error"]) %>% 
    dplyr::mutate(p_value = 2 * stats::pnorm(-abs(coef / se))) %>%     
    dplyr::mutate(lower = coef - (2 * se), upper = coef + (2 * se)) %>% 
    dplyr::mutate_at(.vars = c("coef", "se", "lower", "upper"), exp) %>% 
    dplyr::select(coef, se, lower, upper, p_value)
 
summary_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline
    
# Sanity check!
summary(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)

setequal(summary_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline$p_value, summary(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)$coefficients[, "Pr(>|z|)"])

# Combined
model_oc_post_baseline <- cbind(variable = c("Intercept", "origin_control_community", "same_community", "drive_distance_kilometers"),
      summary_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline) 

model_oc_post_baseline

# *** * ***




# Assemble aic valies to compare models
# ---------------------

# Baseline
aic_oc_baseline <- c(aic = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline)$aic)         
aic_oc_baseline

aic_oc_baseline_2 <- c(aic = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2)$aic)         
aic_oc_baseline_2

# Post baseline
aic_oc_post_baseline <- c(aic = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)$aic)         
aic_oc_post_baseline


# Assemble deviance values
# ---------------------

# Baseline
deviance_oc_baseline <- c(deviance = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline)$deviance)
deviance_oc_baseline

deviance_oc_baseline_2 <- c(deviance = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2)$deviance)
deviance_oc_baseline_2

# Post baseline
deviance_oc_post_baseline <- c(deviance = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)$deviance)
deviance_oc_post_baseline


null.deviance_oc_baseline <- c(null.deviance = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline)$null.deviance)
null.deviance_oc_baseline

null.deviance_oc_baseline_2 <- c(null.deviance = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2)$null.deviance)
null.deviance_oc_baseline_2

null.deviance_oc_post_baseline <- c(null.deviance = summary(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)$null.deviance)
null.deviance_oc_post_baseline


model_oc_aic_baseline <- data.frame(model = c("model_oc"), aic_oc_baseline, deviance_oc_baseline, null.deviance_oc_baseline)
model_oc_aic_baseline

model_oc_aic_baseline_2 <- data.frame(model = c("model_oc"), aic_oc_baseline_2, deviance_oc_baseline_2, null.deviance_oc_baseline_2)
model_oc_aic_baseline_2

model_oc_aic_post_baseline <- data.frame(model = c("model_oc"), aic_oc_post_baseline, deviance_oc_post_baseline, null.deviance_oc_post_baseline)
model_oc_aic_post_baseline

# *** * ***


# Predict expected probability of genetic linkage between communities in the 2011 Botswana census
# ---------------------

# Notes:

# 1: To obtain the expected probability of genetic linkage the offset will be set to 1 (offset = 1, log_offset = log(1) = 0)

# 2: Comparatively, expected (mean) counts can be obtained by multiplying the expected probability of viral linkage by 
#        the exponent of the log_offset: pred * exp(ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb_orig). 
#        Note: the mean counts are real numbers and not integers

# bcpp
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline_log_offset_0 <- m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% 
    dplyr::rename(ln_max_possible_opposite_sex_pairs_in_sample_baseline_orig = ln_max_possible_opposite_sex_pairs_in_sample_baseline,
                  ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb_orig = ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb) %>% 
    dplyr::mutate(ln_max_possible_opposite_sex_pairs_in_sample_baseline = 0,
                  ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb = 0)

m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline_log_offset_0 %>% dplyr::glimpse()

m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline_log_offset_0 %>% utils::str()

# ***

# census
dplyr::glimpse(m_bw_2011_census_pairwise_dist_pop)

m_bw_2011_census_pairwise_dist_pop_postbaseline_log_offset_0 <- m_bw_2011_census_pairwise_dist_pop %>% 
    dplyr::mutate(ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb = 0)
    
m_bw_2011_census_pairwise_dist_pop_postbaseline_log_offset_0 %>% utils::str()

# *** * ***


# Linear model predictions

# Predictions for model: nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline


# bcpp: post baseline, exclude (dest == "Digawana")
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_excl_dest_digawana <- cbind(m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline_log_offset_0 %>% dplyr::filter(!(dest == "Digawana")), predict(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, newdata = m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline_log_offset_0 %>% dplyr::filter(!(dest == "Digawana")) %>% dplyr::select(origin_control_community, same_community, curr_travel_dist_km, ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb), type = "link", se.fit = TRUE))
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_excl_dest_digawana)

utils::str(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_excl_dest_digawana)

# bcpp: post baseline, (dest == "Digawana") only
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_dest_digawana_only <- cbind(m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline_log_offset_0 %>% dplyr::filter(dest == "Digawana"), predict(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, newdata = m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline_log_offset_0 %>% dplyr::filter(dest == "Digawana") %>% dplyr::select(origin_control_community, same_community, curr_travel_dist_km, ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb), type = "link", se.fit = TRUE))
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_dest_digawana_only)

utils::str(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_dest_digawana_only)

# combine datasets
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp <- rbind(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_excl_dest_digawana,
                                                                             pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_dest_digawana_only) %>%
                                                                       dplyr::arrange(dest)

utils::str(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp)

pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp %>% dplyr::select(curr_orig_dest)


# census
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline <- cbind(m_bw_2011_census_pairwise_dist_pop_postbaseline_log_offset_0, predict(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, newdata = m_bw_2011_census_pairwise_dist_pop_postbaseline_log_offset_0 %>% dplyr::select(origin_control_community, same_community, curr_travel_dist_km, ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb), type = "link", se.fit = TRUE))
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)

utils::str(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)

# *** * ***


# Get inverse of link function to convert fitted values to scale of response variable

# We access the inverse of the link function for the fitted negative binomial model object as: object$family$linkinv
# Reminder! The standard link in a negative binomial regression model is log(mu) = xTb hence the inverse is just the exponent: mu = exp(xTb)
ilink_nb_glm <- nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline$family$linkinv
ilink_nb_glm

# Template:

# Note: method 1 and 2 are equivalent

# method 1:
# nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline <- nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline %>% dplyr::mutate(pred = exp(fit), upper = exp(fit + (2 * se.fit)), lower = exp(fit - (2 * se.fit)))

# method 2:
# nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline <- nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline %>% dplyr::mutate(pred = ilink_nb_glm(fit), upper = ilink_nb_glm(fit + (2 * se.fit)), lower = ilink_nb_glm(fit - (2 * se.fit)))


# Linear model predictions converted to response scale

# bcpp
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp <- pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp %>% dplyr::mutate(pred = ilink_nb_glm(fit), upper = ilink_nb_glm(fit + (2 * se.fit)), lower = ilink_nb_glm(fit - (2 * se.fit)))
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp)

utils::str(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp)

# census
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline <- pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline %>% dplyr::mutate(pred = ilink_nb_glm(fit), upper = ilink_nb_glm(fit + (2 * se.fit)), lower = ilink_nb_glm(fit - (2 * se.fit)))
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)

utils::str(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)

# *** * ***




# Save workspace
# ---------------------

save.image("results/estimate_transmission_risk.rda")



# Display session information for record keeping
# ----------------------------------------------
sessionInfo()


# And, thats all folks!




