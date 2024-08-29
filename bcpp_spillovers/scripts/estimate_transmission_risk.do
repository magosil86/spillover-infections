* Author:        Lerato E. Magosi
* Stata version: 13.1
* Platform:      x86_64-apple-darwin17.0 (64-bit)
* Date:          21082023

* ----------------------------------


* Goal: Fit negative binomial regression model on transmission flows between bcpp communities with: same community transmission, control community as the source of transmission and drive distance as predictors


* Note 1. Because there were no post baseline samples in "Digawana" we shall exclude Diawana as a 
*    recipient community in the post baseline model

* Note 2. To account for the sex ratio of males to females im the number of sampled participants 
*    we shall run models with the following offsets:

* Baseline model: 

* Distinct possible opposite-sex (male-female, female-male) hiv transmission pairs in sample is

* (number males sampled at baseline for community i * number females sampled at baseline for community j) + (number females sampled at baseline for community i * number males sampled at baseline for community j)


* Post baseline model: 

* Distinct possible opposite-sex (male-female, female-male) hiv transmission pairs in sample is

* (total number males sampled for community i * number females sampled post baseline for community j) + (total number females sampled for community i * number males sampled post baseline for community j)



* platform: local mac osx 10.15.7



set more off



* Set working dir
*---------------------------------------

cd ./bcpp_spillovers/


* Setting -up data input and output paths 
*---------------------------------------

* filepath to input file
* Load dataset of pairwise drive distances and transmission flows between bcpp communities
local path_2_input data/sub_m_bcpp_pairwise_dist_same_comm_oc_gender_combined.csv

* filepath to output file
local path_2_output results/


* Load data
*---------

insheet using "`path_2_input'", names clear


* Inspect data
*------------
d


* Generate category variable to flag community pairs where transmission events occurred during baseline
*------------

gen baseline = .

recode baseline (.=1) if (num_linked_pairs_observed_baseli > 0)
recode baseline (.=0)

*note: stata orders the categories for us
tab baseline 

*****


* Generate category variable to flag community pairs where transmission events occurred post baseline
*------------

gen post_baseline = .

recode post_baseline (.=1) if (num_linked_pairs_observed_post_b > 0)
recode post_baseline (.=0)

*note: stata orders the categories for us
tab post_baseline 

*****


* View number of participants sampled at baseline versus post baseline
duplicates examples orig number_hosts_sampled_group_1_fem number_hosts_sampled_group_1_mal number_hosts_sampled_group_1_f_b number_hosts_sampled_group_1_m_b number_hosts_sampled_group_1_f_p number_hosts_sampled_group_1_m_p


* Fit negative binomial regression model
*---------------------------------------

* Offsets are set as follows:

* Offset for model with baseline transmission pairs: product of number sampled at baseline for communities i and j

* Offset for model with post baseline transmission pairs: product of total number sampled in community i and post baseline samples in community j

*

* Accounting for the sex ratio, the offsets become:

* 

* Baseline model: Distinct possible opposite-sex (male-female, female-male) hiv transmission pairs in sample is

* (number males sampled at baseline for community i * number females sampled at baseline for community j) + (number females sampled at baseline for community i * number males sampled at baseline for community j)

*

* Post baseline model: Distinct possible opposite-sex (male-female, female-male) hiv transmission pairs in sample is

* (total number males sampled for community i * number females sampled post baseline for community j) + (total number females sampled for community i * number males sampled post baseline for community j)


* Linear models


* rename v43 to ln_max_poss_os_prs_samp_a i.e. ln_max_possible_opposite_sex_pairs_in_sample_all
rename v43 ln_max_poss_os_prs_samp_a
                
* fit model of num_linked_pairs_observed variable, employing intervention status (origin_control_community), same_community and pairwise drive distance as covariates
* offset: distinct possible opposite-sex hiv transmission pairs in sample
*         i.e. offset: (total number males sampled community i * total number females sampled community j) + (total number females sampled community i * total number males sampled community j)
nbreg num_linked_pairs_observed origin_control_community same_community curr_travel_dist_km, offset(ln_max_poss_os_prs_samp_a) vce(robust)

* retrieve stored regression parameters
* for reference see: 
mat li r(table)

test

* view stored results
eret li

* store model results
estimates store m_nb_oc_same_comm_linear_a

*** * ***


* baseline

* rename ln_max_possible_opposite_sex_pai to ln_max_poss_os_prs_samp_b i.e. ln_max_possible_opposite_sex_pairs_in_sample_baseline
rename ln_max_possible_opposite_sex_pai ln_max_poss_os_prs_samp_b

* fit model of num_linked_pairs_observed variable, employing intervention status (origin_control_community), same_community and pairwise drive distance as covariates
* offset: distinct possible opposite-sex hiv transmission pairs in sample
*         i.e. offset is: (number baseline males i * number baseline females j) + (number baseline females i * number baseline males j)
nbreg num_linked_pairs_observed_baseli origin_control_community same_community curr_travel_dist_km, offset(ln_max_poss_os_prs_samp_b) vce(robust)

* retrieve stored regression parameters
* for reference see: 
mat li r(table)

test

* view stored results
eret li

* store model results
estimates store m_nb_oc_same_comm_linear_b

***

* Exclude Digawana community for consistency with post baseline model:

* fit model of num_linked_pairs_observed variable, employing intervention status (origin_control_community), same_community and pairwise drive distance as covariates
* offset: distinct possible opposite-sex hiv transmission pairs in sample
*         i.e. offset is: (number baseline males i * number baseline females j) + (number baseline females i * number baseline males j)
nbreg num_linked_pairs_observed_baseli origin_control_community same_community curr_travel_dist_km if !(dest == "Digawana"), offset(ln_max_poss_os_prs_samp_b) vce(robust)

* retrieve stored regression parameters
* for reference see: 
mat li r(table)

test

* view stored results
eret li

* store model results
estimates store m_nb_oc_same_comm_linear_b2

*** * ***


* post baseline

* rename max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb to max_poss_os_prs_samp_pb2
rename v44 max_poss_os_prs_samp_pb2

* rename ln_max_possible_opposite_sex_pairs_in_sample_post_baseline_all_x_pb to ln_max_poss_os_prs_samp_pb2
rename v45 ln_max_poss_os_prs_samp_pb2

* fit model of num_linked_pairs_observed variable, employing intervention status (origin_control_community), same_community and pairwise drive distance as covariates
* offset: distinct possible opposite-sex hiv transmission pairs in sample
*         i.e. (total number males in community i * post-baseline females in j) + (total number females in community i * post-baseline males in j)
nbreg num_linked_pairs_observed_post_b origin_control_community same_community curr_travel_dist_km if !(dest == "Digawana"), offset(ln_max_poss_os_prs_samp_pb2) vce(robust)

* retrieve stored regression parameters
* for reference see: 
mat li r(table)

test

* view stored results
eret li

* store model results
estimates store m_nb_oc_same_comm_lin_pb2

*** * ***



* Compare models
*---------------------------------------

* Compare oc models
est stat m_nb_oc_same_comm_linear_a m_nb_oc_same_comm_linear_b m_nb_oc_same_comm_linear_b2 m_nb_oc_same_comm_lin_pb2


********** Saving the dataset ***********



*save the full dataset
saveold "`path_2_output'estimate_transmission_risk.dta", replace



********** end of Saving the dataset ***********



di "Stata version `c(version)'"

*And that's all folks


