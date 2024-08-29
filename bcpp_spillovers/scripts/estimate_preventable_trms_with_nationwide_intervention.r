# Author:    Lerato E. Magosi
# R version: 4.1.2 (2021-11-01)
# Platform:  x86_64-apple-darwin17.0 (64-bit)
# Date:      16Mar2024



# Goal: Estimate the maximum proportion of hiv transmissions that could have been 
#         prevented were the intervention applied to all communities 
#         (intervention, control and non-trial) nationally.

# Important!
#       Because there were no post baseline samples collected in "Digawana" 
#       community we shall exclude "Digawana" as a recipient community 
#       from predictions of the risk of transmission.


# Background: Previously fit a negative binomial regression model to estimate the 
#             expected probability of transmission between a pair of individuals
#             randomly sampled from their communities. The model was fit to counts
#             of directed transmission pairs identified within and between bcpp 
#             communities post baseline.
#             For reference, see: bcpp_spillovers/scripts/estimate_transmission_risk.r
#            
#             To set up a counterfactual experiment that estimates the number of 
#             transmissions that could have been prevented were the intervention
#             applied to all communities nationally, we shall run two experiments,
#             first: we shall treat all communities as intervention communities and
#             second: treat all communities as control communities.



# Platform: local macOS


# Load functions --------------------- 

# Load functions to estimate degree of connectedness between communities in a cluster-randomized trial
source("./bcpp_spillovers/scripts/functions_estimate_d_degree_connectedness_between_locations.r")


# Load libraries ---------------------

library(dplyr)         # needed for sorting and merging data.frames
library(tidyr)         # needed for combining or splitting fields
library(MASS)          # for predicting fitted values and fitting negative binomial models
library(doParallel)    # for parallelizing computations
library(foreach)       # for parallelizing computations
library(psych)         # for five point data summaries
library(forcats)       # for manipulating factors
library(RColorBrewer)  # for colour palettes
library(ggplot2)       # for plotting
library(scales)
library(ggrepel)


# Set number of cores
registerDoParallel(cores = 4)


# Set paths
# ---------------------

# set path to project dir
setwd("./bcpp_spillovers")

# check current dir
getwd()


# Load datasets
# ---------------------

# Load hiv prevalence estimates from BAIS IV 2013 by district
hiv_prev_by_district_bais_13 <- read.csv("./data/bais_2013_hiv_prevalence_by_district.csv")

# View data structure
dplyr::glimpse(hiv_prev_by_district_bais_13)

# Load workspace that contains the predicted (expected) risk of hiv transmission between communities post baseline
load("./results/estimate_transmission_risk.rda")
ls()


# View model that shows the risk (expected probability) of hiv transmission between a pair 
#     of individuals randomly sampled from their respective communities post baseline

# Model description: Fit a negative binomial regression model of num_linked_pairs_observed_post_baseline 
#     with origin_control_community, same community transmission and pairwise drive distance as covariates
nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline

# Model summary
summary(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)

# View alpha value used to initialize theta in the neg binomial regression model
# Note. Alpha values obtained from: 
#    bcpp_spillovers/scripts/estimate_transmission_risk.r
alpha_m_nb_oc_same_comm_linear_post_baseline

print(alpha_m_nb_oc_same_comm_linear_post_baseline, digits = 16)

# ***

# View structure of model input dataset
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::glimpse()

m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% utils::str()

# Note: "Digawana" community was excluded from the post baseline model as a recipient community because there were no
#        post baseline samples in "Digawana"

# Subset community pairs with maximum distinct possible pairs: zero
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% 
    dplyr::filter(max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb == 0) %>% 
    dplyr::select(curr_orig_dest)

m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% 
    dplyr::filter(!(dest == "Digawana")) %>%
    utils::str()

# ****


# View predictions of the risk (expected probability) of hiv transmission between communities

# bcpp: input dataset
m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline_log_offset_0 %>% dplyr::glimpse()

m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline_log_offset_0 %>% utils::str()

# bcpp: predictions
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp)

pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp %>% utils::str()

# View risk prediction by intervention for same community transmission (control community source (origin) = 1, intervention community source (origin) = 0)
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp %>% dplyr::filter(same_community == 1) %>% dplyr::count(pred, origin_control_community)


# census: input dataset
m_bw_2011_census_pairwise_dist_pop_postbaseline_log_offset_0 %>% dplyr::glimpse()

m_bw_2011_census_pairwise_dist_pop_postbaseline_log_offset_0 %>% utils::str()

# census: predictions
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline)

pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline %>% utils::str()

# ****


# Expand census input dataset to predict the risk of transmission within census communities

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 <- m_bw_2011_census_pairwise_dist_pop_postbaseline_log_offset_0 %>% 
    # in preparation for setting up same community pairs select fields of interest
    dplyr::select(starts_with("orig"), ends_with("origin"),
                  ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb) %>% 
    # Extract census communities
    dplyr::distinct() %>% 
    # set-up origin and destiantion information for same commuity pairs
    dplyr::mutate(dest = orig, 
                  destination = origin, 
                  destination_district = origin_district, 
                  destination_intervention = origin_intervention, 
                  destination_lat = origin_lat, 
                  destination_lon = origin_lon, 
                  flag_control_destination = flag_control_origin, 
                  BTOTL_2011_destination = BTOTL_2011_origin, 
                  MTOTL_2011_destination = MTOTL_2011_origin, 
                  FTOTL_2011_destination = FTOTL_2011_origin,
                  curr_travel_dist = 0,
                  curr_travel_dist_km = 0, 
                  curr_travel_time = 0,  
                  curr_travel_time_h = 0,
                  ln_curr_travel_dist_km = 0,
                  same_community = 1) %>% 
    tidyr::unite(., col = curr_orig_dest, origin, destination, sep = "~", remove = F) %>%
    base::rbind(.,  m_bw_2011_census_pairwise_dist_pop_postbaseline_log_offset_0)
       

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% dplyr::glimpse()

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% utils::str()

# Sanity check! Expect 488 same community pairs
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% dplyr::count(same_community)


# Attach hiv prevalence information from bais 2013 dataset by district
# ---------------------
    
# Harmonize district names prior to merging
hiv_prev_by_district_bais_13_harmonize_district_names <- hiv_prev_by_district_bais_13 %>% 
    dplyr::mutate(district = gsub("Northeast", "North East", district),
			      district = gsub("Southeast", "South East", district),
			      district = gsub("Selebi Phikwe", "Selibe Phikwe", district),
			      district = gsub("Central Serowe", "Serowe Palapye", district),
			      district = gsub("Ngamiland South", "Ngamiland East", district),
			      district = gsub("Ngamiland North", "Ngamiland West", district))

hiv_prev_by_district_bais_13_harmonize_district_names %>% dplyr::glimpse()

hiv_prev_by_district_bais_13_harmonize_district_names

# census
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 <- m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>%
    # Harmonize district names  
    dplyr::mutate(origin_district = gsub("Ngamiland$", "Ngamiland East", origin_district),
                  destination_district = gsub("Ngamiland$", "Ngamiland East", destination_district)
                  ) %>%
    dplyr::inner_join(., 
                      hiv_prev_by_district_bais_13_harmonize_district_names %>% dplyr::select(district, ends_with("_prop")), 
                      by = c("origin_district" = "district")
                      ) %>%
    dplyr::rename(hivprev_m_prop_origin_district = hiv_prev_male_prop, 
                  hivprev_f_prop_origin_district = hiv_prev_female_prop, 
                  hivprev_total_prop_origin_district = hiv_prev_total_prop) %>%
    dplyr::inner_join(., 
                      hiv_prev_by_district_bais_13_harmonize_district_names %>% dplyr::select(district, ends_with("_prop")), 
                      by = c("destination_district" = "district")
                      ) %>%
    dplyr::rename(hivprev_m_prop_destination_district = hiv_prev_male_prop, 
                  hivprev_f_prop_destination_district = hiv_prev_female_prop, 
                  hivprev_total_prop_destination_district = hiv_prev_total_prop)

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% dplyr::glimpse()

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% utils::str()


# bcpp
m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 <- m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline_log_offset_0 %>%
    # Harmonize district names  
    dplyr::mutate(origin_district = gsub("Ngamiland$", "Ngamiland East", origin_district),
                  destination_district = gsub("Ngamiland$", "Ngamiland East", destination_district)
                  ) %>%
    dplyr::inner_join(., 
                      hiv_prev_by_district_bais_13_harmonize_district_names %>% dplyr::select(district, ends_with("_prop")), 
                      by = c("origin_district" = "district")
                      ) %>%
    dplyr::rename(hivprev_m_prop_origin_district = hiv_prev_male_prop, 
                  hivprev_f_prop_origin_district = hiv_prev_female_prop, 
                  hivprev_total_prop_origin_district = hiv_prev_total_prop) %>%
    dplyr::inner_join(., 
                      hiv_prev_by_district_bais_13_harmonize_district_names %>% dplyr::select(district, ends_with("_prop")), 
                      by = c("destination_district" = "district")
                      ) %>%
    dplyr::rename(hivprev_m_prop_destination_district = hiv_prev_male_prop, 
                  hivprev_f_prop_destination_district = hiv_prev_female_prop, 
                  hivprev_total_prop_destination_district = hiv_prev_total_prop)
                          
m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 %>% dplyr::glimpse()

m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 %>% utils::str()

# ****


# Attach number of people with HIV sampled in each bcpp community to the census dataset
# ---------------------

# Subset bcpp community names and locations
bcpp_comm_names_loc_n30 <- m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% 
    dplyr::filter(!(origin_intervention == "Outside")) %>% 
    dplyr::select(origin, orig) %>% 
    dplyr::distinct()

bcpp_comm_names_loc_n30


# Get number sampled in each bcpp community during baseline and post-baseline
bcpp_number_sampled_orig <- m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 %>% 
    dplyr::select(orig, contains("number_hosts_sampled_group_1")) %>%
    dplyr::distinct() %>%
    # Assemble bcpp locations and number sampled
    dplyr::inner_join(., bcpp_comm_names_loc_n30, by = "orig") %>%
    dplyr::select(orig, origin, everything())

bcpp_number_sampled_orig %>% dplyr::glimpse()
    
bcpp_number_sampled_orig

bcpp_number_sampled_dest <- m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 %>% 
    dplyr::select(dest, contains("number_hosts_sampled_group_2")) %>%
    dplyr::distinct() %>%
    # Assemble bcpp locations and number sampled
    dplyr::inner_join(., bcpp_comm_names_loc_n30, by = c("dest" = "orig")) %>%
    dplyr::rename(destination = origin) %>%
    dplyr::select(dest, destination, everything())

bcpp_number_sampled_dest %>% dplyr::glimpse()
    
bcpp_number_sampled_dest


m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 <- m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>%
    dplyr::left_join(., bcpp_number_sampled_dest %>% dplyr::select(!dest), by = "destination") %>%
    dplyr::rename_with(., ~ gsub("number_hosts_sampled_group_2", "destination_number_hosts_sampled_group_2", .x))
    
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% dplyr::glimpse()

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% utils::str()


m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 <- m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>%
    dplyr::left_join(., bcpp_number_sampled_orig %>% dplyr::select(!orig), by = "origin") %>%
    dplyr::rename_with(., ~ gsub("number_hosts_sampled_group_1", "origin_number_hosts_sampled_group_1", .x))

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% dplyr::glimpse()

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% utils::str()

# Sanity check!

# Expect recipient (destination) communities where the field: destination_number_hosts_sampled_group_2 has a NON-missing value to be bcpp communities 
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% 
    dplyr::filter(!is.na(destination_number_hosts_sampled_group_2)) %>% 
    dplyr::count(destination_intervention)

# Expect source (origin) communities where the field: destination_number_hosts_sampled_group_2 has a NON-missing value to include both bcpp communities (intervention, control) and non-trial communities (outside)
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% 
    dplyr::filter(!is.na(destination_number_hosts_sampled_group_2)) %>% 
    dplyr::count(origin_intervention)



# Counterfactual - Predict what the risk of hiv transmission to recipients in communities would have been
#     were the bcpp intervention applied to all communities nationally. Note that risk here describes the 
#     expected probability of genetic linkage between a pair of individuals randomly sampled from their 
#     respective communities 


# Treat all communities as control communities
#  ---------------------


# bcpp counterfactual

m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 %>% dplyr::glimpse()

m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 %>% utils::str()

# To setup the counterfactual experiment we shall assign all origin communities in community pairings control status
m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_ctrl <- m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 %>%
    dplyr::rename(origin_control_community_original = origin_control_community) %>%
    # Treat all communities as control communities
    dplyr::mutate(origin_control_community = 1) %>%
    # Because there were no post baseline samples in "Digawana" we shall exclude 
    # "Digawana" as a recipient community  
    dplyr::filter(!(dest == "Digawana"))

    
m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_ctrl %>% dplyr::glimpse()

m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_ctrl %>% utils::str()

# Sanity check! Expect: n = 870
m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_ctrl %>% dplyr::count(origin_control_community)


# census counterfactual

dplyr::glimpse(m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144)

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% utils::str()

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl <- m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>%
    dplyr::rename(origin_control_community_original = origin_control_community) %>%
    # Treat all communities as control communities
    dplyr::mutate(origin_control_community = 1) %>%
    # Because there were no post baseline samples in "Digawana" we shall exclude 
    # "Digawana" as a recipient community  
    dplyr::filter(!(dest == "Digawana"))
                     
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl %>% dplyr::glimpse()

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl %>% utils::str()

# Sanity check! Expect: n = 238,144 - 488 = 237,656
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl %>% dplyr::count(origin_control_community)

# ****


# Define function to get bootstrapped predictions of the expected probability of hiv transmission between communities

# Acknowledgements: https://www.r-bloggers.com/2015/12/prediction-intervals-for-poisson-regression/

get_negbin_prediction_expected_prob_hiv_trm_bootstrap <- function(model_in, alpha_model_in, model_data_in, prediction_data_in, n_bootstraps = 1000) {

	# Set seed for overall procedure
	set.seed(1691)

	# Get individual seed values for resampling the model data
	vec_seeds <- base::round(stats::runif(n = n_bootstraps, min = 1, max = 1000), digits = 0)

	# Get bootstrapped predictions of the expected response rate

	# Note: The individual iterations are run in parallel then assembled at
	#       the end with rbind

	boot_expected_response_rate <- foreach(i = 1:n_bootstraps, .combine = rbind) %dopar% {

		# Set seed for individual iteration
		set.seed(vec_seeds[i])

		# Get number of rows in model_data
		n_rows <- nrow(model_data_in)

		# Generate a sequence of the row numbers in model_data
		vec_row_numbers <- base::seq(n_rows)

		# Sample row numbers of rows to subset in model_data with replacement 
		subset_rows <- base::sample(x = vec_row_numbers, size = n_rows, replace = TRUE)

		# Subset rows in model_data to generate a bootstrap dataset
		bootstrap_data <- model_data_in[ subset_rows, ]

		# Fit model of num_linked_pairs_observed with 
		#     origin_control_community, same community transmission and pairwise drive distance as covariates

		# 'update' will update the model call stored in the model object and re-fit the model
		# Here the data object in the model call is updated before re-fitting the model
		updated_model <- stats::update(object = model_in, data = bootstrap_data, init.theta = 1/alpha_model_in, evaluate = TRUE)

		# Get predictions of response variable in the response scale
		bootstrap_pred <- stats::predict(updated_model, newdata = prediction_data_in, type = "response")

        bootstrap_pred
        
	}


    boot_expected_response_rate


}


# Call function:

# bcpp counterfactual:
res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_bcpp_post_baseline_cf_ctrl <- get_negbin_prediction_expected_prob_hiv_trm_bootstrap(model_in = nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, 
                                                                                                                                              alpha_model_in = alpha_m_nb_oc_same_comm_linear_post_baseline, 
                                                                                                                                              model_data_in = m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::filter(!(dest == "Digawana")), 
                                                                                                                                              prediction_data_in = m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_ctrl, 
                                                                                                                                              n_bootstraps = 1000)

# View results

# Note: The function returns a matrix of predicted response values in the response scale 
#       Response variable is the expected prob. of hiv trm between community pairs 

dplyr::glimpse(res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_bcpp_post_baseline_cf_ctrl)

res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_bcpp_post_baseline_cf_ctrl[1:2, 1:5]

# Transpose i.e. switch rows and cols
t(res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_bcpp_post_baseline_cf_ctrl[1:2, 1:5])


# census counterfactual:
res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_ctrl <- get_negbin_prediction_expected_prob_hiv_trm_bootstrap(model_in = nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, 
                                                                                                                                                alpha_model_in = alpha_m_nb_oc_same_comm_linear_post_baseline, 
                                                                                                                                                model_data_in = m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::filter(!(dest == "Digawana")), 
                                                                                                                                                prediction_data_in = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl, 
                                                                                                                                                n_bootstraps = 1000)

# View results
dplyr::glimpse(res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_ctrl)

res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_ctrl[1:2, 1:5]

# Transpose i.e. switch rows and cols
t(res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_ctrl[1:2, 1:5])



# Estimate hiv transmissions between census communities
#  ---------------------

# Define function:
estimate_hiv_trms_cs11_cf <- function(df_cs11_prediction_data_in, vec_cs11_bootstrap_pred_in, df_cs11_community_size_and_hivprev_in, verbose_output = TRUE) {

	# Attach predicted responses to prediction dataset
	# ---------------------

    # census	
	pred_census_cf <- base::cbind(df_cs11_prediction_data_in, vec_cs11_bootstrap_pred_in) %>%
	    dplyr::rename(pred_cf = vec_cs11_bootstrap_pred_in)
	
	if (verbose_output) {
	
		base::writeLines("\n\npred_census_cf %>% dplyr::glimpse(): ")
		pred_census_cf %>% dplyr::glimpse()
	
	} 


	# Estimate number of hivpos and number of transmissions between census communities        
	# ---------------------
	
	# If hiv prevalence info is already available in the dataset that was used to 
	#     predict the risk of hiv transmission
	if (is.null(df_cs11_community_size_and_hivprev_in)) {

		cs11_hiv_trms_cf <- pred_census_cf %>%
		    # Because we are interested in transmissions between census communities
		    #     let's exclude same community pairings
		    dplyr::filter(!(origin == destination)) %>%
			# Get number of people with hiv in each community
			# Note: we will multiply the estimated community size by the estimated hiv prevalence in each district
			dplyr::mutate(est_hivpos_origin = BTOTL_2011_origin * hivprev_total_prop_origin_district,
			              est_hivpos_f_origin = FTOTL_2011_origin * hivprev_f_prop_origin_district,
			              est_hivpos_m_origin = MTOTL_2011_origin * hivprev_m_prop_origin_district,
						  est_hivpos_destination = BTOTL_2011_destination * hivprev_total_prop_destination_district,
						  est_hivpos_f_destination = FTOTL_2011_destination * hivprev_f_prop_destination_district,
						  est_hivpos_m_destination = MTOTL_2011_destination * hivprev_m_prop_destination_district) %>%
			# Get distinct possible pairs in population between communities
			dplyr::mutate(max_possible_opposite_sex_pairs_in_cs11_population_original = ((est_hivpos_m_origin * est_hivpos_f_destination) + (est_hivpos_f_origin * est_hivpos_m_destination))) %>%
			# Tip! Because 5 * NA = NA and 5 + NA = NA community pairs with a non-trial recipient community will have: max_possible_opposite_sex_pairs_in_cs11_population = NA and est_hiv_trms_cf = NA
			dplyr::mutate(max_possible_opposite_sex_pairs_in_cs11_population = ((est_hivpos_m_origin * destination_number_hosts_sampled_group_2_f_post_baseline) + (est_hivpos_f_origin * destination_number_hosts_sampled_group_2_m_post_baseline))) %>%
			# Get estimated transmission events between origin-destination community pairs
			# Note: This is the product of the estimated probability of hiv transmission between
			#       two individuals randomly sampled from their respective communities and the maximum
			#       distinct possible pairs between those communities
			dplyr::mutate(est_hiv_trms_cf = max_possible_opposite_sex_pairs_in_cs11_population * pred_cf) %>%
			# Assign intervention status
			dplyr::mutate(
				flag_control_origin_num = case_when(
					flag_control_origin == "Control" ~ 0, 
					flag_control_origin == "Intervention" ~ 1
					),
				flag_control_destination_num = case_when(
					flag_control_destination == "Control" ~ 0, 
					flag_control_destination == "Intervention" ~ 1
					)      
				)
		
	} else {
	
		cs11_hiv_trms_cf <- df_cs11_community_size_and_hivprev_in %>%
			# Attach counterfactual predictions and intervals
			dplyr::inner_join(., pred_census_cf %>%
									 # Because we are interested in transmissions between census communities
									 #     let's exclude same community pairings
									 dplyr::filter(!(origin == destination)) %>%
									 dplyr::select(origin, destination, pred_cf),
							  by = c("origin", "destination")) %>%
			# Get estimated transmission events between origin-destination community pairs
			# Note: This is the product of the estimated probability of hiv transmission between
			#       two individuals randomly sampled from their respective communities and the maximum
			#       distinct possible pairs between those communities
			dplyr::mutate(est_hiv_trms_cf = max_possible_opposite_sex_pairs_in_cs11_population * pred_cf) %>%
			# Assign intervention status
			dplyr::mutate(
				flag_control_origin_num = case_when(
					flag_control_origin == "Control" ~ 0, 
					flag_control_origin == "Intervention" ~ 1
					),
				flag_control_destination_num = case_when(
					flag_control_destination == "Control" ~ 0, 
					flag_control_destination == "Intervention" ~ 1
					)      
				)
	
	}

    if (verbose_output) {
        
        base::writeLines("cs11_hiv_trms_cf %>% dplyr::glimpse(): ")
        cs11_hiv_trms_cf %>% dplyr::glimpse()
    
    }

   
    # Estimate degree of connectedness between census communities by genetic linkage
    # ---------------------

	cs11_get_d_genetic_cf <- get_degree_contact(method = "genetic",
												source_in = cs11_hiv_trms_cf$origin, 
												recipient_in = cs11_hiv_trms_cf$destination, 
												source_intervention_in = cs11_hiv_trms_cf$flag_control_origin_num, 
												recipient_intervention_in = cs11_hiv_trms_cf$flag_control_destination_num, 
												linked_pairs_in = cs11_hiv_trms_cf$est_hiv_trms_cf,
												verbose_output = TRUE)

    if (verbose_output) {
    
        base::writeLines("cs11_get_d_genetic_cf %>% dplyr::glimpse(): ")
        cs11_get_d_genetic_cf %>% dplyr::glimpse()
    
    }


	# Attach estimates of degree of connectedness between census communities
	# ---------------------

	cs11_hiv_trms_cf_d_genetic <- cs11_hiv_trms_cf %>%
		dplyr::select(origin, orig, origin_lon, origin_lat, origin_district, origin_intervention, flag_control_origin, flag_control_origin_num) %>% 
		dplyr::distinct() %>%
		# Attach estimated hiv transmissions originating from intervention communities
		dplyr::inner_join(., cs11_get_d_genetic_cf$d_directed_ref_intervention, by = c("origin" = "recipient")) %>%
		dplyr::rename(d_genetic_ref_intervention_cf = total_inflow) %>%
		# Attach estimated hiv transmissions originating from control communities
		dplyr::inner_join(., cs11_get_d_genetic_cf$d_directed_ref_control, by = c("origin" = "recipient")) %>%
		dplyr::rename(d_genetic_ref_control_cf = total_inflow)

    if (verbose_output) {
    
        base::writeLines("cs11_hiv_trms_cf_d_genetic %>% dplyr::glimpse(): ")
        cs11_hiv_trms_cf_d_genetic %>% dplyr::glimpse()
    
    }


	# Subset community pairs that contain bcpp communities
	#      Why? We are interested in the influence of non-trial communities on bcpp communities
	#      but only have seroconverter information for bcpp communities
	#
	#      Note: This represents the total inflow into bcpp communities from bcpp communities and non-trial communities
    # ---------------------
    
	sub_bcpp_cs11_hiv_trms_cf_d_genetic <- cs11_hiv_trms_cf_d_genetic %>% dplyr::filter(origin_intervention %in% c("Intervention", "Control")) 

    if (verbose_output) {
    
        base::writeLines("sub_bcpp_cs11_hiv_trms_cf_d_genetic %>% dplyr::glimpse(): ")
        sub_bcpp_cs11_hiv_trms_cf_d_genetic %>% dplyr::glimpse()
    
    }


	# Estimate degree of connectedness between bcpp communities by genetic linkage
	# ---------------------

	# Subset bcpp community pairs
	sub_bcpp_cs11_hiv_trms_cf <- cs11_hiv_trms_cf %>% dplyr::filter(origin_intervention %in% c("Intervention", "Control") & destination_intervention %in% c("Intervention", "Control"))
	
    if (verbose_output) {
    
        base::writeLines("sub_bcpp_cs11_hiv_trms_cf %>% dplyr::glimpse(): ")
        sub_bcpp_cs11_hiv_trms_cf %>% dplyr::glimpse()
    
    }

	bcpp_get_d_genetic_cf <- get_degree_contact(method = "genetic",
											    source_in = sub_bcpp_cs11_hiv_trms_cf$origin, 
											    recipient_in = sub_bcpp_cs11_hiv_trms_cf$destination, 
											    source_intervention_in = sub_bcpp_cs11_hiv_trms_cf$flag_control_origin_num, 
											    recipient_intervention_in = sub_bcpp_cs11_hiv_trms_cf$flag_control_destination_num, 
											    linked_pairs_in = sub_bcpp_cs11_hiv_trms_cf$est_hiv_trms_cf,
											    verbose_output = TRUE)

    if (verbose_output) {
    
        base::writeLines("bcpp_get_d_genetic_cf %>% dplyr::glimpse(): ")
        bcpp_get_d_genetic_cf %>% dplyr::glimpse()
    
    }


	# Attach estimates of degree of connectedness between bcpp communities
	# ---------------------
	
	bcpp_hiv_trms_cf_d_genetic <- cs11_hiv_trms_cf %>%
		dplyr::filter(origin_intervention %in% c("Intervention", "Control") & destination_intervention %in% c("Intervention", "Control")) %>%
		dplyr::select(origin, orig, origin_lon, origin_lat, origin_district, origin_intervention, flag_control_origin, flag_control_origin_num) %>% 
		dplyr::distinct() %>%
		# Attach estimated hiv transmissions originating from intervention communities
		dplyr::inner_join(., bcpp_get_d_genetic_cf$d_directed_ref_intervention, by = c("origin" = "recipient")) %>%
		dplyr::rename(d_genetic_ref_intervention_bcpp_cf = total_inflow) %>%
		# Attach estimated hiv transmissions originating from control communities
		dplyr::inner_join(., bcpp_get_d_genetic_cf$d_directed_ref_control, by = c("origin" = "recipient")) %>%
		dplyr::rename(d_genetic_ref_control_bcpp_cf = total_inflow)

    if (verbose_output) {
    
        base::writeLines("bcpp_hiv_trms_cf_d_genetic %>% dplyr::glimpse(): ")
        bcpp_hiv_trms_cf_d_genetic %>% dplyr::glimpse()
    
    }


    # List of items to return
    # ---------------------
    
    output <- list("bcpp_hiv_trms_cf_d_genetic" = bcpp_hiv_trms_cf_d_genetic,
                   "sub_bcpp_cs11_hiv_trms_cf_d_genetic" = sub_bcpp_cs11_hiv_trms_cf_d_genetic,
                   "cs11_hiv_trms_cf_d_genetic" = cs11_hiv_trms_cf_d_genetic)

    output

}


# Call function:

# View dataset with district hiv prevalence estimates
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl %>% dplyr::glimpse()

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl %>% utils::str()

# census counterfactual:
lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl <- base::apply(X = res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_ctrl, MARGIN = 1, FUN = estimate_hiv_trms_cs11_cf, df_cs11_prediction_data_in = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = F)

# view results
utils::str(lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl$result.1) 

# lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl[1]

lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl[[1]]$bcpp_hiv_trms_cf_d_genetic

lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl[[1]]$sub_bcpp_cs11_hiv_trms_cf_d_genetic



# Sanity check! Expect 29 rows with non-zero values (trial communities) and 458 with non-zero values (non-trial communities)
# Why? Because zero post-baseline samples were collected from non-trial communities

lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl[1]$result.1$cs11_hiv_trms_cf_d_genetic %>% 
    dplyr::filter(d_genetic_ref_control_cf > 0) %>% 
    dplyr::select() %>%
    base::nrow()

lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl[1]$result.1$cs11_hiv_trms_cf_d_genetic %>% 
    dplyr::filter(d_genetic_ref_control_cf == 0) %>% 
    dplyr::select(origin_intervention) %>% 
    base::nrow()

# sort results
lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl[1]$result.1$cs11_hiv_trms_cf_d_genetic %>% 
    dplyr::arrange(d_genetic_ref_control_cf, origin_intervention) %>% 
    dplyr::select(origin, origin_intervention, d_genetic_ref_intervention_cf, d_genetic_ref_control_cf) %>% 
    utils::head()

lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl[1]$result.1$cs11_hiv_trms_cf_d_genetic %>% 
    dplyr::arrange(d_genetic_ref_control_cf, origin_intervention) %>% 
    dplyr::select(origin, origin_intervention, d_genetic_ref_intervention_cf, d_genetic_ref_control_cf) %>% 
    utils::tail()

lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl[1]$result.1$cs11_hiv_trms_cf_d_genetic %>% 
    dplyr::arrange(-d_genetic_ref_control_cf, origin_intervention) %>% 
    dplyr::select(origin, origin_intervention, d_genetic_ref_intervention_cf, d_genetic_ref_control_cf)


# Estimate hiv transmissions within census communities
#  ---------------------

# Define function:
estimate_hiv_trms_within_comm_cs11_cf <- function(df_cs11_prediction_data_in, vec_cs11_bootstrap_pred_in, df_cs11_community_size_and_hivprev_in, verbose_output = TRUE) {

	# Attach predicted responses to prediction dataset
	# ---------------------

    # census	
	pred_census_cf <- base::cbind(df_cs11_prediction_data_in, vec_cs11_bootstrap_pred_in) %>%
	    dplyr::rename(pred_cf = vec_cs11_bootstrap_pred_in)
	
	if (verbose_output) {
	
		base::writeLines("\n\npred_census_cf %>% dplyr::glimpse(): ")
		pred_census_cf %>% dplyr::glimpse()
	
	} 


	# Estimate number of hivpos and number of transmissions within census communities        
	# ---------------------
	
	# If hiv prevalence info is already available in the dataset that was used to 
	#     predict the risk of hiv transmission
	if (is.null(df_cs11_community_size_and_hivprev_in)) {

		within_cs11_hiv_trms_cf <- pred_census_cf %>%
		    # Subset same community pairs
		    dplyr::filter(origin == destination) %>%
			# Get number of people with hiv in each community
			# Note: we will multiply the estimated community size by the estimated hiv prevalence in each district
			dplyr::mutate(est_hivpos_origin = BTOTL_2011_origin * hivprev_total_prop_origin_district,
			              est_hivpos_f_origin = FTOTL_2011_origin * hivprev_f_prop_origin_district,
			              est_hivpos_m_origin = MTOTL_2011_origin * hivprev_m_prop_origin_district,
						  est_hivpos_destination = BTOTL_2011_destination * hivprev_total_prop_destination_district,
						  est_hivpos_f_destination = FTOTL_2011_destination * hivprev_f_prop_destination_district,
						  est_hivpos_m_destination = MTOTL_2011_destination * hivprev_m_prop_destination_district) %>%
			# Get distinct possible pairs in population between communities
			# Note: This is within community transmission because origin and destination communities are the same
			dplyr::mutate(max_possible_opposite_sex_pairs_in_cs11_population_original = ((est_hivpos_m_origin * est_hivpos_f_destination) + (est_hivpos_f_origin * est_hivpos_m_destination))) %>%
			# Tip! Because 5 * NA = NA and 5 + NA = NA community pairs with a non-trial recipient community will have: max_possible_opposite_sex_pairs_in_cs11_population = NA and est_hiv_trms_cf = NA
			dplyr::mutate(max_possible_opposite_sex_pairs_in_cs11_population = ((est_hivpos_m_origin * destination_number_hosts_sampled_group_2_f_post_baseline) + (est_hivpos_f_origin * destination_number_hosts_sampled_group_2_m_post_baseline))) %>%
			# Get estimated transmission events between origin-destination community pairs
			# Note: This is the product of the estimated probability of hiv transmission between
			#       two individuals randomly sampled from their respective communities and the maximum
			#       distinct possible pairs between those communities
			dplyr::mutate(within_comm_est_hiv_trms_cf = max_possible_opposite_sex_pairs_in_cs11_population * pred_cf)
		
	} else {
	
		within_cs11_hiv_trms_cf <- df_cs11_community_size_and_hivprev_in %>%
			# Attach counterfactual predictions and intervals
			dplyr::inner_join(., pred_census_cf %>%
			                         # Subset same commuinty pairs
                                     dplyr::filter(origin == destination) %>%
									 dplyr::select(origin, destination, pred_cf),
							  by = c("origin", "destination")) %>%
			# Get estimated transmission events between origin-destination community pairs
			# Note: This is the product of the estimated probability of hiv transmission between
			#       two individuals randomly sampled from their respective communities and the maximum
			#       distinct possible pairs between those communities
			dplyr::mutate(within_comm_est_hiv_trms_cf = max_possible_opposite_sex_pairs_in_cs11_population * pred_cf)

	
	}

    if (verbose_output) {
        
        base::writeLines("within_cs11_hiv_trms_cf %>% dplyr::glimpse(): ")
        within_cs11_hiv_trms_cf %>% dplyr::glimpse()
    
    }


    within_cs11_hiv_trms_cf


}


# Call function:

# census counterfactual: 
lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl <- base::apply(X = res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_ctrl, MARGIN = 1, FUN = estimate_hiv_trms_within_comm_cs11_cf, df_cs11_prediction_data_in = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = F)

# view results
utils::str(lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl$result.1) 
       
# lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl[1]

lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl[[1]] %>% utils::head()

# view sorted results
lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl[1]$result.1 %>% 
    dplyr::arrange(within_comm_est_hiv_trms_cf) %>% dplyr::select(origin, origin_intervention, within_comm_est_hiv_trms_cf, max_possible_opposite_sex_pairs_in_cs11_population) %>% 
    utils::head()

lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl[1]$result.1 %>% 
    dplyr::arrange(within_comm_est_hiv_trms_cf) %>% 
    dplyr::select(origin, origin_intervention, within_comm_est_hiv_trms_cf, max_possible_opposite_sex_pairs_in_cs11_population) %>% 
    utils::tail()


# Sanity check! Compare with result below from function: estimate_hiv_trms_bcpp_cf()
lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl[1]$result.1 %>% 
    dplyr::arrange(within_comm_est_hiv_trms_cf) %>% 
    dplyr::select(origin, origin_intervention, within_comm_est_hiv_trms_cf, max_possible_opposite_sex_pairs_in_cs11_population) %>% 
    dplyr::filter(origin_intervention %in% c("Intervention", "Control"))



# Estimate hiv transmissions within bcpp communities
#  ---------------------

# Define function:
estimate_hiv_trms_bcpp_cf <- function(df_bcpp_prediction_data_in, vec_bcpp_bootstrap_pred_in, df_cs11_community_size_and_hivprev_in, verbose_output = TRUE) {

	# Attach predicted responses to prediction dataset
	# ---------------------

    # bcpp	
	pred_bcpp_cf <- base::cbind(df_bcpp_prediction_data_in, vec_bcpp_bootstrap_pred_in) %>%
	    dplyr::rename(pred_cf = vec_bcpp_bootstrap_pred_in)
	
	if (verbose_output) {
	
		base::writeLines("\n\npred_bcpp_cf %>% dplyr::glimpse(): ")
		pred_bcpp_cf %>% dplyr::glimpse()
	
	} 


	# Estimate number of transmissions within bcpp communities
	# ---------------------

	# If hiv prevalence info is already available in the dataset that was used to 
	#     predict the risk of hiv transmission
	if (is.null(df_cs11_community_size_and_hivprev_in)) {

		within_bcpp_hiv_trms_cf <- pred_bcpp_cf %>%
			dplyr::filter(orig == dest) %>%
			# Get number of people with hiv in each community
			# Note: we will multiply the estimated community size by the estimated hiv prevalence in each district
			dplyr::mutate(est_hivpos_origin = BTOTL_2011_origin * hivprev_total_prop_origin_district,
			              est_hivpos_f_origin = FTOTL_2011_origin * hivprev_f_prop_origin_district,
			              est_hivpos_m_origin = MTOTL_2011_origin * hivprev_m_prop_origin_district,
						  est_hivpos_destination = BTOTL_2011_destination * hivprev_total_prop_destination_district,
						  est_hivpos_f_destination = FTOTL_2011_destination * hivprev_f_prop_destination_district,
						  est_hivpos_m_destination = MTOTL_2011_destination * hivprev_m_prop_destination_district) %>%
			# Get distinct possible pairs in population between communities
			# Note: This is within community transmission because origin and destination communities are the same
			dplyr::mutate(max_possible_opposite_sex_pairs_in_cs11_population_original = ((est_hivpos_m_origin * est_hivpos_f_destination) + (est_hivpos_f_origin * est_hivpos_m_destination))) %>%
			# Tip! Because 5 * NA = NA and 5 + NA = NA community pairs with a non-trial recipient community will have: max_possible_opposite_sex_pairs_in_cs11_population = NA and est_hiv_trms_cf = NA
			dplyr::mutate(max_possible_opposite_sex_pairs_in_cs11_population = ((est_hivpos_m_origin * number_hosts_sampled_group_2_f_post_baseline) + (est_hivpos_f_origin * number_hosts_sampled_group_2_m_post_baseline))) %>%
			# Get estimated transmission events between origin-destination community pairs
			# Note: This is the product of the estimated probability of hiv transmission between
			#       two individuals randomly sampled from their respective communities and the maximum
			#       distinct possible pairs between those communities
			dplyr::mutate(within_comm_est_hiv_trms_cf = max_possible_opposite_sex_pairs_in_cs11_population * pred_cf)
	
	} else {

		within_bcpp_hiv_trms_cf <- pred_bcpp_cf %>%
			dplyr::filter(orig == dest) %>%
			dplyr::inner_join(., df_cs11_community_size_and_hivprev_in %>%
									 dplyr::filter(!(orig == "Otse" & est_hivpos_origin == 412.797)) %>%
									 dplyr::select(orig, hivprev_total_prop_origin_district) %>%
									 dplyr::distinct(),
							  by = "orig"
							  ) %>%	
			# Get number of people with hiv in each community
			# Note: we will multiply the estimated community size by the estimated hiv prevalence in each district
			dplyr::mutate(est_hivpos_origin = BTOTL_2011_origin * hivprev_total_prop_origin_district,
			              est_hivpos_f_origin = FTOTL_2011_origin * hivprev_f_prop_origin_district,
			              est_hivpos_m_origin = MTOTL_2011_origin * hivprev_m_prop_origin_district,
						  est_hivpos_destination = BTOTL_2011_destination * hivprev_total_prop_destination_district,
						  est_hivpos_f_destination = FTOTL_2011_destination * hivprev_f_prop_destination_district,
						  est_hivpos_m_destination = MTOTL_2011_destination * hivprev_m_prop_destination_district) %>%
			# Get distinct possible pairs in population between communities
			# Note: This is within community transmission because origin and destination communities are the same
			dplyr::mutate(max_possible_opposite_sex_pairs_in_cs11_population_original = ((est_hivpos_m_origin * est_hivpos_f_destination) + (est_hivpos_f_origin * est_hivpos_m_destination))) %>%
			# Tip! Because 5 * NA = NA and 5 + NA = NA community pairs with a non-trial recipient community will have: max_possible_opposite_sex_pairs_in_cs11_population = NA and est_hiv_trms_cf = NA
			dplyr::mutate(max_possible_opposite_sex_pairs_in_cs11_population = ((est_hivpos_m_origin * number_hosts_sampled_group_2_f_post_baseline) + (est_hivpos_f_origin * number_hosts_sampled_group_2_m_post_baseline))) %>%
			# Get estimated transmission events between origin-destination community pairs
			# Note: This is the product of the estimated probability of hiv transmission between
			#       two individuals randomly sampled from their respective communities and the maximum
			#       distinct possible pairs between those communities
			dplyr::mutate(within_comm_est_hiv_trms_cf = max_possible_opposite_sex_pairs_in_cs11_population * pred_cf)	
	
	}


	if (verbose_output) {
	
		base::writeLines("\n\nwithin_bcpp_hiv_trms_cf %>% dplyr::glimpse(): ")
		within_bcpp_hiv_trms_cf %>% dplyr::glimpse()
	
	} 


    within_bcpp_hiv_trms_cf

}


# Call function:

# bcpp counterfactual
lst_res_estimate_hiv_trms_bcpp_post_baseline_cf_ctrl <- base::apply(X = res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_bcpp_post_baseline_cf_ctrl, MARGIN = 1, FUN = estimate_hiv_trms_bcpp_cf, df_bcpp_prediction_data_in = m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_ctrl, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = F)

utils::str(lst_res_estimate_hiv_trms_bcpp_post_baseline_cf_ctrl$result.1) 
       
lst_res_estimate_hiv_trms_bcpp_post_baseline_cf_ctrl[1]$result.1 %>% 
    dplyr::arrange(within_comm_est_hiv_trms_cf) %>% 
    dplyr::select(orig, orig_intervention, within_comm_est_hiv_trms_cf, max_possible_opposite_sex_pairs_in_cs11_population)


# END: Treat all communities as control communities i.e. estimate trms in the absence of the bcpp intervention -----------------------


# Treat all communities as intervention communities
#  ---------------------


# bcpp counterfactual

m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 %>% dplyr::glimpse()

m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 %>% utils::str()

# To setup the counterfactual experiment we shall assign all origin communities in community pairings intervention status
m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_int <- m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0 %>%
    dplyr::rename(origin_control_community_original = origin_control_community) %>%
    # Treat all communities as intervention communities
    dplyr::mutate(origin_control_community = 0) %>%
    # Because there were no post baseline samples in "Digawana" we shall exclude 
    # "Digawana" as a recipient community  
    dplyr::filter(!(dest == "Digawana"))

    
m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_int %>% dplyr::glimpse()

m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_int %>% utils::str()

# Sanity check! Expect: n = 870
m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_int %>% dplyr::count(origin_control_community)


# census counterfactual

dplyr::glimpse(m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144)

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>% utils::str()

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int <- m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_n238144 %>%
    dplyr::rename(origin_control_community_original = origin_control_community) %>%
    # Treat all communities as intervention communities
    dplyr::mutate(origin_control_community = 0) %>%
    # Because there were no post baseline samples in "Digawana" we shall exclude 
    # "Digawana" as a recipient community  
    dplyr::filter(!(dest == "Digawana"))

                     
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int %>% dplyr::glimpse()

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int %>% utils::str()

# Sanity check! Expect: n = 238,144 - 488 = 237,656
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int %>% dplyr::count(origin_control_community)

# ****


# Get bootstrapped predictions of the expected probability of hiv transmission between communities

# Acknowledgements: https://www.r-bloggers.com/2015/12/prediction-intervals-for-poisson-regression/

# Call function: get_negbin_prediction_expected_prob_hiv_trm_bootstrap()

# bcpp counterfactual:
res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_bcpp_post_baseline_cf_int <- get_negbin_prediction_expected_prob_hiv_trm_bootstrap(model_in = nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, 
                                                                                                                                             alpha_model_in = alpha_m_nb_oc_same_comm_linear_post_baseline, 
                                                                                                                                             model_data_in = m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::filter(!(dest == "Digawana")), 
                                                                                                                                             prediction_data_in = m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_int, 
                                                                                                                                             n_bootstraps = 1000)

# View results

# Note: The function returns a matrix of predicted response values in the response scale 
#       Response variable is the expected prob. of hiv trm between community pairs 

dplyr::glimpse(res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_bcpp_post_baseline_cf_int)

res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_bcpp_post_baseline_cf_int[1:2, 1:5]

# Transpose i.e. switch rows and cols
t(res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_bcpp_post_baseline_cf_int[1:2, 1:5])


# census counterfactual:
res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_int <- get_negbin_prediction_expected_prob_hiv_trm_bootstrap(model_in = nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, 
                                                                                                                                               alpha_model_in = alpha_m_nb_oc_same_comm_linear_post_baseline, 
                                                                                                                                               model_data_in = m_bcpp_pairwise_dist_pop_flag_baseline_vs_postbaseline %>% dplyr::filter(!(dest == "Digawana")), 
                                                                                                                                               prediction_data_in = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int, 
                                                                                                                                               n_bootstraps = 1000)

# View results
dplyr::glimpse(res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_int)

res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_int[1:2, 1:5]

# Transpose i.e. switch rows and cols
t(res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_int[1:2, 1:5])



# Estimate hiv transmissions between census communities
#  ---------------------

# Call function: estimate_hiv_trms_cs11_cf()

# View dataset with district hiv prevalence estimates
m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int %>% dplyr::glimpse()

m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int %>% utils::str()

# census counterfactual:
lst_res_estimate_hiv_trms_cs11_post_baseline_cf_int <- base::apply(X = res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_int, MARGIN = 1, FUN = estimate_hiv_trms_cs11_cf, df_cs11_prediction_data_in = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = F)

# view results
utils::str(lst_res_estimate_hiv_trms_cs11_post_baseline_cf_int$result.1) 
       
#lst_res_estimate_hiv_trms_cs11_post_baseline_cf_int[1]

lst_res_estimate_hiv_trms_cs11_post_baseline_cf_int[[1]]$bcpp_hiv_trms_cf_d_genetic

lst_res_estimate_hiv_trms_cs11_post_baseline_cf_int[[1]]$sub_bcpp_cs11_hiv_trms_cf_d_genetic

# sort results
lst_res_estimate_hiv_trms_cs11_post_baseline_cf_int[1]$result.1$cs11_hiv_trms_cf_d_genetic %>% 
    dplyr::arrange(d_genetic_ref_control_cf, origin_intervention) %>% dplyr::select(origin, origin_intervention, d_genetic_ref_intervention_cf, d_genetic_ref_control_cf) %>% 
    utils::head()

lst_res_estimate_hiv_trms_cs11_post_baseline_cf_int[1]$result.1$cs11_hiv_trms_cf_d_genetic %>% 
    dplyr::arrange(d_genetic_ref_control_cf, origin_intervention) %>% dplyr::select(origin, origin_intervention, d_genetic_ref_intervention_cf, d_genetic_ref_control_cf) %>% 
    utils::tail()

lst_res_estimate_hiv_trms_cs11_post_baseline_cf_int[1]$result.1$cs11_hiv_trms_cf_d_genetic %>% 
    dplyr::arrange(d_genetic_ref_control_cf, origin_intervention) %>% 
    dplyr::select(origin, origin_intervention, d_genetic_ref_intervention_cf, d_genetic_ref_control_cf)


# Estimate hiv transmissions within census communities
#  ---------------------

# Call function: estimate_hiv_trms_within_comm_cs11_cf()

# census counterfactual: 
lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int <- base::apply(X = res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_census_post_baseline_cf_int, MARGIN = 1, FUN = estimate_hiv_trms_within_comm_cs11_cf, df_cs11_prediction_data_in = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = F)

# view results
utils::str(lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int$result.1) 
       
# lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int[1]

# view sorted results
lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int[1]$result.1 %>% 
    dplyr::arrange(within_comm_est_hiv_trms_cf) %>% 
    dplyr::select(origin, origin_intervention, within_comm_est_hiv_trms_cf, max_possible_opposite_sex_pairs_in_cs11_population) %>% 
    utils::head()

lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int[1]$result.1 %>% 
    dplyr::arrange(within_comm_est_hiv_trms_cf) %>% 
    dplyr::select(origin, origin_intervention, within_comm_est_hiv_trms_cf, max_possible_opposite_sex_pairs_in_cs11_population) %>% 
    utils::tail()


# Sanity check!

# Estimated trms within comms in the absence of the bcpp intervention
a <- lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl[1]$result.1 %>% 
    dplyr::arrange(within_comm_est_hiv_trms_cf) %>% 
    dplyr::select(origin, origin_intervention, within_comm_est_hiv_trms_cf, max_possible_opposite_sex_pairs_in_cs11_population) %>% 
    dplyr::select(origin, within_comm_est_hiv_trms_cf)

a %>% utils::tail()

a %>% utils::head()

# Estimated trms if all communities recieved the intervention
b <- lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int[1]$result.1 %>% 
    dplyr::arrange(within_comm_est_hiv_trms_cf) %>% 
    dplyr::select(origin, origin_intervention, within_comm_est_hiv_trms_cf, max_possible_opposite_sex_pairs_in_cs11_population) %>%  
    dplyr::select(origin, within_comm_est_hiv_trms_cf)
    
b %>% utils::tail()

b %>% utils::head()


# Estimate hiv transmissions within bcpp communities
#  ---------------------

# Call function: estimate_hiv_trms_bcpp_cf()

# bcpp counterfactual
lst_res_estimate_hiv_trms_bcpp_post_baseline_cf_int <- base::apply(X = res_get_negbin_prediction_expected_prob_hiv_trm_bootstrap_bcpp_post_baseline_cf_int, MARGIN = 1, FUN = estimate_hiv_trms_bcpp_cf, df_bcpp_prediction_data_in = m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_int, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = F)

utils::str(lst_res_estimate_hiv_trms_bcpp_post_baseline_cf_int$result.1) 
       
lst_res_estimate_hiv_trms_bcpp_post_baseline_cf_int[1]$result.1 %>% 
    dplyr::arrange(within_comm_est_hiv_trms_cf) %>% 
    dplyr::select(orig, orig_intervention, within_comm_est_hiv_trms_cf, max_possible_opposite_sex_pairs_in_cs11_population)

# END: Treat all communities as intervention communities i.e. estimate trms were the bcpp intervention applied to all communities -----------------------



# Estimate maximum proportion of preventable hiv transmissions
# ---------------------

# Define function:

get_max_prop_preventable_hiv_trms_cf <- function(within_cs11_hiv_trms_cf_ctrl_in, within_cs11_hiv_trms_cf_int_in, btwn_cs11_hiv_trms_cf_ctrl_in, btwn_cs11_hiv_trms_cf_int_in, verbose_output = TRUE) {

    # Assemble dataset
    # ---------------------
    
    # Estimates of degree of connectedness between census communities in the absence of the bcpp intervention
    cs11_hiv_trms_cf_d_genetic_ctrl <- btwn_cs11_hiv_trms_cf_ctrl_in$cs11_hiv_trms_cf_d_genetic
        
    # Estimates of degree of connectedness between census communities in the presence of the bcpp intervention
    cs11_hiv_trms_cf_d_genetic_int <- btwn_cs11_hiv_trms_cf_int_in$cs11_hiv_trms_cf_d_genetic

 
 	if (verbose_output) {
	
		base::writeLines("\n\ncs11_hiv_trms_cf_d_genetic_ctrl %>% dplyr::glimpse(): ")
		cs11_hiv_trms_cf_d_genetic_ctrl %>% dplyr::glimpse()

		base::writeLines("\n\ncs11_hiv_trms_cf_d_genetic_int %>% dplyr::glimpse(): ")
		cs11_hiv_trms_cf_d_genetic_int %>% dplyr::glimpse()

		base::writeLines("\n\nwithin_cs11_hiv_trms_cf_ctrl_in %>% dplyr::glimpse(): ")
        within_cs11_hiv_trms_cf_ctrl_in %>% utils::str()
        
		base::writeLines("\n\nwithin_cs11_hiv_trms_cf_int_in %>% dplyr::glimpse(): ")
	    within_cs11_hiv_trms_cf_int_in %>% utils::str()
	    
	} 
    
       
	# Estimate maximum number of preventable hiv infections
	# ---------------------

	# Assemble dataset
	cs11_hiv_trms_cf_preventable <- within_cs11_hiv_trms_cf_ctrl_in %>%
		                     dplyr::rename(pred_cf_ctrl = pred_cf, 
		                                   max_possible_opposite_sex_pairs_in_cs11_population_within_comm_ctrl = max_possible_opposite_sex_pairs_in_cs11_population,
		                                   within_comm_est_hiv_trms_cf_ctrl = within_comm_est_hiv_trms_cf) %>%
		# Attach estimates of degree of connectedness between census communities in the absence of the bcpp intervention
		dplyr::inner_join(., cs11_hiv_trms_cf_d_genetic_ctrl %>%
		                         dplyr::rename(d_genetic_ref_intervention_cf_ctrl = d_genetic_ref_intervention_cf,
		                                       d_genetic_ref_control_cf_ctrl = d_genetic_ref_control_cf) %>%
								 dplyr::select(origin, origin_district,
											   d_genetic_ref_intervention_cf_ctrl, 
											   d_genetic_ref_control_cf_ctrl),
						 by = c("origin", "origin_district")) %>%
		# Attach estimates of transmissions within census communities in the presence of the bcpp intervention
		dplyr::inner_join(., within_cs11_hiv_trms_cf_int_in %>%
		                     dplyr::rename(pred_cf_int = pred_cf, 
		                                   max_possible_opposite_sex_pairs_in_cs11_population_within_comm_int = max_possible_opposite_sex_pairs_in_cs11_population,
		                                   within_comm_est_hiv_trms_cf_int = within_comm_est_hiv_trms_cf) %>%
		                     dplyr::select(origin, origin_district, pred_cf_int, 
		                                   max_possible_opposite_sex_pairs_in_cs11_population_within_comm_int, within_comm_est_hiv_trms_cf_int),
						 by = c("origin", "origin_district")) %>%		
		# Attach estimates of degree of connectedness between census communities in the presence of the bcpp intervention
		dplyr::inner_join(., cs11_hiv_trms_cf_d_genetic_int %>%
		                         dplyr::rename(d_genetic_ref_intervention_cf_int = d_genetic_ref_intervention_cf,
		                                       d_genetic_ref_control_cf_int = d_genetic_ref_control_cf) %>%
							 dplyr::select(origin, origin_district, 
										   d_genetic_ref_intervention_cf_int, 
										   d_genetic_ref_control_cf_int),
						   by = c("origin", "origin_district")) %>%
		# Get estimates of total transmissions
        dplyr::mutate(total_est_hiv_trms_int = within_comm_est_hiv_trms_cf_int + d_genetic_ref_intervention_cf_int + d_genetic_ref_control_cf_int,
                      total_est_hiv_trms_ctrl = within_comm_est_hiv_trms_cf_ctrl + d_genetic_ref_intervention_cf_ctrl + d_genetic_ref_control_cf_ctrl,
                      diff_total_est_hiv_trms = total_est_hiv_trms_ctrl - total_est_hiv_trms_int) %>%
		# Get maximum proportion of preventable transmissions                
		dplyr::mutate(prop_preventable_hiv_trms_cf = diff_total_est_hiv_trms / total_est_hiv_trms_ctrl)


 	if (verbose_output) {
	
		base::writeLines("\n\ncs11_hiv_trms_cf_preventable %>% utils::str(): ")
		cs11_hiv_trms_cf_preventable %>% utils::str()
	
	} 
    
    
    cs11_hiv_trms_cf_preventable
    

}





# Call function:

# res_get_max_prop_preventable_hiv_trms_post_baseline_cf <- base::mapply(FUN = get_max_prop_preventable_hiv_trms_cf, within_cs11_hiv_trms_cf_ctrl_in = lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl[1:1], within_cs11_hiv_trms_cf_int_in = lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int[1:1], btwn_cs11_hiv_trms_cf_ctrl_in = lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl[1:1], btwn_cs11_hiv_trms_cf_int_in = lst_res_estimate_hiv_trms_cs11_post_baseline_cf_int[1:1])

res_get_max_prop_preventable_hiv_trms_post_baseline_cf <- base::mapply(FUN = get_max_prop_preventable_hiv_trms_cf, within_cs11_hiv_trms_cf_ctrl_in = lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl, within_cs11_hiv_trms_cf_int_in = lst_res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int, btwn_cs11_hiv_trms_cf_ctrl_in = lst_res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl, btwn_cs11_hiv_trms_cf_int_in = lst_res_estimate_hiv_trms_cs11_post_baseline_cf_int)


# Get 2.5% and 97.5% quantiles of the difference in hiv transmissions in the presence and absence of the bcpp intervention across all bootstraps
# ---------------------

get_quantiles_diff_hiv_trms_cf <- function(lst_dataset_in, significance_level_in = 0.05) {

	# Initialize data.frame()
	df <- lst_dataset_in[, 1]$diff_total_est_hiv_trms

    # Get number of bootstraps
    n_bootstraps <- ncol(lst_dataset_in)
    
	# Assemble diff_total_est_hiv_trms fields in each bootstrap into a single data.frame with cbind
	for (number in 2:n_bootstraps) { 

		df <- cbind(df, lst_dataset_in[, number]$diff_total_est_hiv_trms)

	}
    
    df %>% utils::str()
        
	lower_threshold <- significance_level_in / 2
	upper_threshold <- 1 - lower_threshold


	# Get 2.5% and 97.5% quantiles of difference in hiv transmissions in the presence and absence of the bcpp intervention across all bootstraps
	df_quantiles <- base::t( base::apply(X = df, MARGIN = 1, FUN = quantile, c(lower_threshold, upper_threshold), type = 7, digits = 7, na.rm = TRUE))
	
	
	df_quantiles


}

# Call function:

res_get_quantiles_diff_hiv_trms_post_baseline_cf <- get_quantiles_diff_hiv_trms_cf(lst_dataset_in = res_get_max_prop_preventable_hiv_trms_post_baseline_cf, significance_level_in = 0.05)


# View results

head(res_get_quantiles_diff_hiv_trms_post_baseline_cf)

res_get_quantiles_diff_hiv_trms_post_baseline_cf %>% 
    base::as.data.frame() %>% 
    dplyr::filter(!is.na(`2.5%`)) %>% 
    utils::head()
 
res_get_quantiles_diff_hiv_trms_post_baseline_cf %>% 
    base::as.data.frame() %>% 
    dplyr::filter(!is.na(`2.5%`))


 # Sanity check! 
 
# Expect: 29 rows 
res_get_quantiles_diff_hiv_trms_post_baseline_cf %>% 
    base::as.data.frame() %>% 
    dplyr::filter(!is.na(`2.5%`)) %>%
    base::nrow()

# Expect: 29
res_get_quantiles_diff_hiv_trms_post_baseline_cf %>% 
    base::as.data.frame() %>% 
    dplyr::filter(!is.na(`2.5%`) & !is.na(`97.5%`)) %>% 
    base::nrow()

# Expect: 458
res_get_quantiles_diff_hiv_trms_post_baseline_cf %>% 
    base::as.data.frame() %>% 
    dplyr::filter(is.na(`2.5%`)) %>% 
    base::nrow()
    


# Get 2.5% and 97.5% quantiles of the mean of the maximum proportion of preventable hiv transmissions across all bootstraps
# ---------------------

get_quantiles_mean_max_prop_preventable_hiv_trms_cf <- function(lst_dataset_in, significance_level_in = 0.05) {

	# Initialize data frames
    
    # Total inflow of trms from all communities in the presence of the bcpp intervention
    df_total_est_hiv_trms_int <- lst_dataset_in[,1]$total_est_hiv_trms_int
    
    # Total inflow of trms from all communities in the absence of the bcpp intervention
    df_total_est_hiv_trms_ctrl <- lst_dataset_in[,1]$total_est_hiv_trms_ctrl
    
    # Get number of bootstraps
    n_bootstraps <- ncol(lst_dataset_in)
    
	# Assemble prop_preventable_hiv_trms_cf fields in each bootstrap into a single data.frame with cbind
	for (number in 2:n_bootstraps) { 

		df_total_est_hiv_trms_int <- cbind(df_total_est_hiv_trms_int, lst_dataset_in[, number]$total_est_hiv_trms_int)

        df_total_est_hiv_trms_ctrl <- cbind(df_total_est_hiv_trms_ctrl, lst_dataset_in[, number]$total_est_hiv_trms_ctrl)
	}

	lower_threshold <- significance_level_in / 2
	upper_threshold <- 1 - lower_threshold


    # Get mean of maximum proportion of preventable hiv transmissions within bootstraps 
    # i.e. get mean of each column

    # Get sum of hiv trms in the presence of the bcpp intervention within bootstraps
    vec_sum_df_total_est_hiv_trms_int <- base::apply(X = df_total_est_hiv_trms_int, MARGIN = 2, FUN = sum, na.rm = TRUE)
    
    writeLines("vec_sum_df_total_est_hiv_trms_int: ")
    vec_sum_df_total_est_hiv_trms_int %>% dplyr::glimpse()
        
    # Get sum of hiv trms in the absence of the bcpp intervention within bootstraps
    vec_sum_df_total_est_hiv_trms_ctrl <- base::apply(X = df_total_est_hiv_trms_ctrl, MARGIN = 2, FUN = sum, na.rm = TRUE)
    
    writeLines("vec_sum_df_total_est_hiv_trms_ctrl: ")
    vec_sum_df_total_est_hiv_trms_ctrl %>% dplyr::glimpse()
    
    writeLines("vec_diff: ")
    vec_diff <- (vec_sum_df_total_est_hiv_trms_ctrl - vec_sum_df_total_est_hiv_trms_int)
    vec_diff %>% dplyr::glimpse() 
    
    writeLines("vec_mean: ")
    vec_mean <- (vec_sum_df_total_est_hiv_trms_ctrl - vec_sum_df_total_est_hiv_trms_int) / vec_sum_df_total_est_hiv_trms_ctrl
    vec_mean %>% dplyr::glimpse() 
    
        
	# Get 2.5% and 97.5% quantiles of mean maximum proportion of preventable hiv transmissions across all bootstraps
	vec_quantiles_mean <- stats::quantile(x = vec_mean, probs = c(lower_threshold, upper_threshold), type = 7, digits = 7)

	# Get 2.5% and 97.5% quantiles of difference in hiv transmissions across all bootstraps
	vec_quantiles_diff <- stats::quantile(x = vec_diff, probs = c(lower_threshold, upper_threshold), type = 7, digits = 7)
		
	output <- list("vec_quantiles_mean" = vec_quantiles_mean,
	               "vec_quantiles_diff" = vec_quantiles_diff)

    output


}


# Call function:

res_get_quantiles_mean_max_prop_preventable_hiv_trms_post_baseline_cf <- get_quantiles_mean_max_prop_preventable_hiv_trms_cf(lst_dataset_in = res_get_max_prop_preventable_hiv_trms_post_baseline_cf, significance_level_in = 0.05)

res_get_quantiles_mean_max_prop_preventable_hiv_trms_post_baseline_cf



# Get 2.5% and 97.5% quantiles of the mean of the maximum proportion of preventable hiv transmissions across all bootstraps
# ---------------------

get_quantiles_mean_max_prop_preventable_hiv_trms_bcpp_cf <- function(lst_dataset_in, significance_level_in = 0.05) {

	# Initialize data frames
    
    # Total inflow of trms from all communities in the presence of the bcpp intervention
    df_total_est_hiv_trms_int <- lst_dataset_in[,1]$total_est_hiv_trms_int
    
    # Total inflow of trms from all communities in the absence of the bcpp intervention
    df_total_est_hiv_trms_ctrl <- lst_dataset_in[,1]$total_est_hiv_trms_ctrl
    
    # Get community names
    df_cs11_communities <- lst_dataset_in[,1] %>%
        as.data.frame() %>% 
        dplyr::select(origin, orig, origin_intervention)
    
    # Get number of bootstraps
    n_bootstraps <- ncol(lst_dataset_in)
    
	# Assemble prop_preventable_hiv_trms_cf fields in each bootstrap into a single data.frame with cbind
	for (number in 2:n_bootstraps) { 

		df_total_est_hiv_trms_int <- cbind(df_total_est_hiv_trms_int, lst_dataset_in[, number]$total_est_hiv_trms_int)

        df_total_est_hiv_trms_ctrl <- cbind(df_total_est_hiv_trms_ctrl, lst_dataset_in[, number]$total_est_hiv_trms_ctrl)
	}

    # Subset transmissions into bcpp communities in the presence of the bcpp intervention
    df_total_est_hiv_trms_int_bcpp <- df_cs11_communities %>% 
        cbind(., df_total_est_hiv_trms_int) %>%
        dplyr::filter(origin_intervention %in% c("Intervention", "Control")) %>%
        dplyr::select(!c("origin", "orig", "origin_intervention"))
    
    writeLines("\n\ndf_total_est_hiv_trms_int_bcpp %>% dplyr::glimpse(): ")
    df_total_est_hiv_trms_int_bcpp %>% dplyr::glimpse()
    
    writeLines("\n\ndf_total_est_hiv_trms_int_bcpp %>% utils::str(): ")
    df_total_est_hiv_trms_int_bcpp %>% utils::str()
    

    # Subset transmissions into bcpp communities in the absence of the bcpp intervention
    df_total_est_hiv_trms_ctrl_bcpp <- df_cs11_communities %>% 
        cbind(., df_total_est_hiv_trms_ctrl) %>%
        dplyr::filter(origin_intervention %in% c("Intervention", "Control")) %>%
        dplyr::select(!c("origin", "orig", "origin_intervention"))
    
    writeLines("\n\ndf_total_est_hiv_trms_ctrl_bcpp %>% dplyr::glimpse(): ")
    df_total_est_hiv_trms_ctrl_bcpp %>% dplyr::glimpse()
    
    writeLines("\n\ndf_total_est_hiv_trms_ctrl_bcpp %>% utils::str(): ")
    df_total_est_hiv_trms_ctrl_bcpp %>% utils::str()

	lower_threshold <- significance_level_in / 2
	upper_threshold <- 1 - lower_threshold


    # Get mean of maximum proportion of preventable hiv transmissions within bootstraps 
    # i.e. get mean of each column

    # Get sum of hiv trms in the presence of the bcpp intervention within bootstraps
    vec_sum_df_total_est_hiv_trms_int <- base::apply(X = df_total_est_hiv_trms_int_bcpp, MARGIN = 2, FUN = sum, na.rm = TRUE)
    
    writeLines("vec_sum_df_total_est_hiv_trms_int: ")
    vec_sum_df_total_est_hiv_trms_int %>% dplyr::glimpse()
        
    # Get sum of hiv trms in the absence of the bcpp intervention within bootstraps
    vec_sum_df_total_est_hiv_trms_ctrl <- base::apply(X = df_total_est_hiv_trms_ctrl_bcpp, MARGIN = 2, FUN = sum, na.rm = TRUE)
    
    writeLines("vec_sum_df_total_est_hiv_trms_ctrl: ")
    vec_sum_df_total_est_hiv_trms_ctrl %>% dplyr::glimpse()
    
    writeLines("vec_diff: ")
    vec_diff <- (vec_sum_df_total_est_hiv_trms_ctrl - vec_sum_df_total_est_hiv_trms_int)
    vec_diff %>% dplyr::glimpse() 
    
    writeLines("vec_mean: ")
    vec_mean <- (vec_sum_df_total_est_hiv_trms_ctrl - vec_sum_df_total_est_hiv_trms_int) / vec_sum_df_total_est_hiv_trms_ctrl
    vec_mean %>% dplyr::glimpse() 
    
        
	# Get 2.5% and 97.5% quantiles of mean maximum proportion of preventable hiv transmissions across all bootstraps
	vec_quantiles_mean <- stats::quantile(x = vec_mean, probs = c(lower_threshold, upper_threshold), type = 7, digits = 7)

	# Get 2.5% and 97.5% quantiles of difference in hiv transmissions across all bootstraps
	vec_quantiles_diff <- stats::quantile(x = vec_diff, probs = c(lower_threshold, upper_threshold), type = 7, digits = 7)
		
	output <- list("vec_quantiles_mean" = vec_quantiles_mean,
	               "vec_quantiles_diff" = vec_quantiles_diff)

    output


}


# Call function:

res_get_quantiles_mean_max_prop_preventable_hiv_trms_post_baseline_bcpp_cf <- get_quantiles_mean_max_prop_preventable_hiv_trms_bcpp_cf(lst_dataset_in = res_get_max_prop_preventable_hiv_trms_post_baseline_cf, significance_level_in = 0.05)

res_get_quantiles_mean_max_prop_preventable_hiv_trms_post_baseline_bcpp_cf


# Important! Because samples were not collected in non-trial communities we expect the following two functions:
#            get_quantiles_mean_max_prop_preventable_hiv_trms_cf() and get_quantiles_mean_max_prop_preventable_hiv_trms_bcpp_cf() to give the same results.
#            Had there been sampling in non-trial communities the function: get_quantiles_mean_max_prop_preventable_hiv_trms_bcpp_cf() would report
#           results for the bcpp communities only and the function: get_quantiles_mean_max_prop_preventable_hiv_trms_cf() would report results for 
#           all 488 census communities.
           

res_get_quantiles_mean_max_prop_preventable_hiv_trms_post_baseline_cf

res_get_quantiles_mean_max_prop_preventable_hiv_trms_post_baseline_bcpp_cf

# ********** END: Bootstrap intervals *************


# Get point estimate of the maximum proportion of preventable hiv transmissions


# Counterfactual - Predict expected probability of hiv transmission between census communities on the response scale
# ---------------------

# Predictions for model: nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline

# bcpp: control
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_cf_ctrl <- cbind(m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_ctrl, stats::predict(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, newdata = m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_ctrl, type = "response", se.fit = TRUE))
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_cf_ctrl)

pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_cf_ctrl %>% utils::str()

# bcpp: intervention
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_cf_int <- cbind(m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_int, stats::predict(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, newdata = m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_int, type = "response", se.fit = TRUE))
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_cf_int)

pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_cf_int %>% utils::str()


# census: control
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_census_cf_ctrl <- cbind(m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl, stats::predict(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, newdata = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl, type = "response", se.fit = TRUE))
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_census_cf_ctrl)

pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_census_cf_ctrl %>% utils::str()

# census: intervention
pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_census_cf_int <- cbind(m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int, stats::predict(nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline, newdata = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int, type = "response", se.fit = TRUE))
dplyr::glimpse(pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_census_cf_int)

pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_census_cf_int %>% utils::str()

# ****


# Estimate hiv transmissions between census communities
#  ---------------------

# Call function: estimate_hiv_trms_cs11_cf()

# control
res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl_point_est <- estimate_hiv_trms_cs11_cf(df_cs11_prediction_data_in = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl, vec_cs11_bootstrap_pred_in = pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_census_cf_ctrl$fit, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = TRUE)

res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl_point_est %>% utils::str()

# intervention
res_estimate_hiv_trms_cs11_post_baseline_cf_int_point_est <- estimate_hiv_trms_cs11_cf(df_cs11_prediction_data_in = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int, vec_cs11_bootstrap_pred_in = pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_census_cf_int$fit, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = TRUE)

res_estimate_hiv_trms_cs11_post_baseline_cf_int_point_est %>% utils::str()


# Estimate hiv transmissions within census communities
#  ---------------------

# Call function: estimate_hiv_trms_within_comm_cs11_cf()

# control
res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl_point_est <- estimate_hiv_trms_within_comm_cs11_cf(df_cs11_prediction_data_in = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_ctrl, vec_cs11_bootstrap_pred_in = pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_census_cf_ctrl$fit, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = TRUE)

res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl_point_est %>% utils::str()

res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl_point_est %>% dplyr::filter(!(destination_intervention == "Outside")) %>% utils::str()

# intervention
res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int_point_est <- estimate_hiv_trms_within_comm_cs11_cf(df_cs11_prediction_data_in = m_bw_2011_census_pairwise_dist_pop_hiv_prev_postbaseline_log_offset_0_cf_int, vec_cs11_bootstrap_pred_in = pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_census_cf_int$fit, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = TRUE)

res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int_point_est %>% utils::str()

res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int_point_est %>% dplyr::filter(!(destination_intervention == "Outside")) %>% utils::str()


# Estimate hiv transmissions within bcpp communities
#  ---------------------

# Call function: estimate_hiv_trms_bcpp_cf()

# control
res_estimate_hiv_trms_bcpp_cf_point_est_ctrl <- estimate_hiv_trms_bcpp_cf(df_bcpp_prediction_data_in = m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_ctrl, vec_bcpp_bootstrap_pred_in = pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_cf_ctrl$fit, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = TRUE)

res_estimate_hiv_trms_bcpp_cf_point_est_ctrl %>% utils::str()

# intervention
res_estimate_hiv_trms_bcpp_cf_point_est_int <- estimate_hiv_trms_bcpp_cf(df_bcpp_prediction_data_in = m_bcpp_pairwise_dist_pop_hiv_prev_flag_baseline_vs_postbaseline_log_offset_0_cf_int, vec_bcpp_bootstrap_pred_in = pred_nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline_bcpp_cf_int$fit, df_cs11_community_size_and_hivprev_in = NULL, verbose_output = TRUE)

res_estimate_hiv_trms_bcpp_cf_point_est_int %>% utils::str()


# Estimate maximum proportion of preventable hiv transmissions
# ---------------------

# Call function: get_max_prop_preventable_hiv_trms_cf()
res_get_max_prop_preventable_hiv_trms_post_baseline_cf_point_est <- get_max_prop_preventable_hiv_trms_cf(within_cs11_hiv_trms_cf_ctrl_in = res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_ctrl_point_est, 
                                                                                                         within_cs11_hiv_trms_cf_int_in = res_estimate_hiv_trms_within_comm_cs11_post_baseline_cf_int_point_est, 
                                                                                                         btwn_cs11_hiv_trms_cf_ctrl_in = res_estimate_hiv_trms_cs11_post_baseline_cf_ctrl_point_est, 
                                                                                                         btwn_cs11_hiv_trms_cf_int_in = res_estimate_hiv_trms_cs11_post_baseline_cf_int_point_est)

res_get_max_prop_preventable_hiv_trms_post_baseline_cf_point_est %>% utils::str()

res_get_max_prop_preventable_hiv_trms_post_baseline_cf_point_est %>% 
    dplyr::filter(!(destination_intervention == "Outside")) %>% 
    utils::str()


mean_max_prop_preventable_hiv_trms_post_baseline_cf_point_est <- res_get_max_prop_preventable_hiv_trms_post_baseline_cf_point_est %>%
	# Get mean maximum proportion of preventable transmissions                
	dplyr::summarise(mean_prop_preventable_hiv_trms_cf = (sum(total_est_hiv_trms_ctrl) - sum(total_est_hiv_trms_int)) / sum(total_est_hiv_trms_ctrl))

mean_max_prop_preventable_hiv_trms_post_baseline_cf_point_est


sub_bcpp_mean_max_prop_preventable_hiv_trms_post_baseline_cf_point_est <- res_get_max_prop_preventable_hiv_trms_post_baseline_cf_point_est %>%
     # Subset bcpp communities
    dplyr::filter(origin_intervention %in% c("Intervention", "Control")) %>%
	# Get mean maximum proportion of preventable transmissions                
	dplyr::summarise(mean_prop_preventable_hiv_trms_cf = (sum(total_est_hiv_trms_ctrl) - sum(total_est_hiv_trms_int)) / sum(total_est_hiv_trms_ctrl))

sub_bcpp_mean_max_prop_preventable_hiv_trms_post_baseline_cf_point_est


# Assemble estimates of the maximum proportion of preventable hiv transmissions and bootstrap intervals
# ---------------------

# Individual community estimates
point_est_and_bootstrap_intervals_post_baseline_cf_cs11_n487 <- cbind(res_get_max_prop_preventable_hiv_trms_post_baseline_cf_point_est, res_get_quantiles_diff_hiv_trms_post_baseline_cf) 
point_est_and_bootstrap_intervals_post_baseline_cf_cs11_n487 %>% utils::head()

point_est_and_bootstrap_intervals_post_baseline_cf_bcpp_n29 <- point_est_and_bootstrap_intervals_post_baseline_cf_cs11_n487 %>% dplyr::filter(!(origin_intervention == "Outside"))
point_est_and_bootstrap_intervals_post_baseline_cf_bcpp_n29

sub_point_est_and_bootstrap_intervals_post_baseline_cf_cs11_n487 <- point_est_and_bootstrap_intervals_post_baseline_cf_cs11_n487 %>% 
    dplyr::select(orig, origin, total_est_hiv_trms_int, total_est_hiv_trms_ctrl, diff_total_est_hiv_trms, "2.5%", "97.5%", prop_preventable_hiv_trms_cf) %>%
    dplyr::rename(diff_total_est_hiv_trms_lower_0_025 = "2.5%", diff_total_est_hiv_trms_upper_0_975 = "97.5%")

sub_point_est_and_bootstrap_intervals_post_baseline_cf_cs11_n487 %>% utils::head()

sub_point_est_and_bootstrap_intervals_post_baseline_cf_bcpp_n29 <- point_est_and_bootstrap_intervals_post_baseline_cf_bcpp_n29 %>% 
    dplyr::filter(!(origin_intervention == "Outside")) %>%
    dplyr::select(orig, origin, total_est_hiv_trms_int, total_est_hiv_trms_ctrl, diff_total_est_hiv_trms, 
                  "2.5%", "97.5%", prop_preventable_hiv_trms_cf, 
                  destination_number_hosts_sampled_group_2_f_post_baseline, destination_number_hosts_sampled_group_2_m_post_baseline) %>%
    dplyr::rename(diff_total_est_hiv_trms_lower_0_025 = "2.5%", diff_total_est_hiv_trms_upper_0_975 = "97.5%")

sub_point_est_and_bootstrap_intervals_post_baseline_cf_bcpp_n29


# All community estimates

# all
mean_point_est_and_bootstrap_intervals_post_baseline_cf <- mean_max_prop_preventable_hiv_trms_post_baseline_cf_point_est %>%
    base::c(., res_get_quantiles_mean_max_prop_preventable_hiv_trms_post_baseline_cf$vec_quantiles_mean)
    
mean_point_est_and_bootstrap_intervals_post_baseline_cf

base::names(mean_point_est_and_bootstrap_intervals_post_baseline_cf) <- c("mean", "lower_0_025", "upper_0_975")

mean_point_est_and_bootstrap_intervals_post_baseline_cf    

# Note: NA is because there was no sampling during the bcpp trial in non-trial communities

# bcpp only
sub_bcpp_mean_point_est_and_bootstrap_intervals_post_baseline_cf <- sub_bcpp_mean_max_prop_preventable_hiv_trms_post_baseline_cf_point_est %>%
    base::c(., res_get_quantiles_mean_max_prop_preventable_hiv_trms_post_baseline_bcpp_cf$vec_quantiles_mean)
    
sub_bcpp_mean_point_est_and_bootstrap_intervals_post_baseline_cf

base::names(sub_bcpp_mean_point_est_and_bootstrap_intervals_post_baseline_cf) <- c("mean", "lower_0_025", "upper_0_975")

sub_bcpp_mean_point_est_and_bootstrap_intervals_post_baseline_cf    


# ********** END: Get point estimate of the maximum proportion of preventable hiv transmissions **********






# Save dataset
# -----------------

# Mean estimate of maximum proportion of preventable transmissions for bcpp communities only
filepath_maxprev <- "./results"
write.csv(sub_bcpp_mean_point_est_and_bootstrap_intervals_post_baseline_cf, file = file.path(filepath_maxprev, "preventable_trms_nationwide_intervention.csv"), row.names = F)




# Display session information for record keeping
# ----------------------------------------------
sessionInfo()


# And, thats all folks!
