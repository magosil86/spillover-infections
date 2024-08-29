# Author:    Lerato E. Magosi
# R version: 4.1.2 (2021-11-01)
# Platform:  x86_64-apple-darwin17.0 (64-bit)
# Date:      25Feb2022

# Goal: Estimate the degree of contact between clusters in a cluster-randomized trial

# Why? Degree of contact between clusters will be used to quantify spillover indirect effects

# Strategy: Degree of contact will be inferred from physical proximity or extent of genetic-linkage between clusters.
#           For example, consider a cluster-randomized trial where the randomization unit is communities i.e. a community-randomized trial.
#           For each location (or community), we will estimate the degree of contact that community has relative to other communities in the reference trial arm.
#           Either the intervention communities or control communities can be selected as the reference trial arm. The degree of contact is then
#           estimated with two methods:
#
# 1. Genetic linkage: For each location the extent of genetic linkage will be measured as the 
#                     number of directed transmission events into that location that originate
#                     from locations in the reference trial arm.

# 2. Distance radius or disc measure: For each location the disc measure will be computed as the
#                     number of locations in the reference trial arm that are within a specified
#                     distance radius of that location.
#
# Acknowledgements: section 2.3.1, Anaya-Izquierdo and Alexander, 2020; Biometrics (PMID: 32557560)



# Required libraries ---------------------------

library(dplyr)     # needed for sorting and merging data.frames
library(tidyr)     # needed for combining or splitting fields
library(stats)     # needed for computing Euclidean distances
library(utils)     # needed for the following functions: str, head, tail
library(sf)        # needed for computing great circle distances



# Function: get_degree_contact ---------------------------
#
# Wrapper function to estimate degree of contact between communities with
#     counts of directed transmission events or spatial point data (latitude and longitude values)
#
#
# parameters:
#
# method                    (character)    A string to indicate the method to use to estimate degree of contact between communities. Options: "genetic", "distance_radius".
# verbose_output            (boolean)      A value to produce detailed output (Defualt is FALSE).
#
# Additional required parameters if method to estimate degree of contact is: "genetic"
#
# source_in                 (character)    A vector of source location names. 
# recipient_in              (character)    A vector of recipient location names.
# source_intervention_in    (numeric)      A vector of trial arm assignments for source locations (1 = intervention, 0 = control).
# recipient_intervention_in (numeric)      A vector of trial arm assignments for recipient locations (1 = intervention, 0 = control).
# linked_pairs_in           (numeric)      A vector of counts of directed transmission events between source and recipient locations.
# 
# Additional required parameters if method to estimate degree of contact is: "distance_radius".
#           
# longitude_in              (numeric)      A vector of longitude values for locations.
# latitude_in               (numeric)      A vector of latitude values for locations.
# intervention_in           (numeric)      A vector of trial arm assignments for locations (1 = intervention, 0 = control).
# crs_in                    (numeric)      A value to indicate the coordinate reference system for supplied latitude and longitude coordinates.                  
# radius_in                 (numeric)      A value to indicate the distance radius within which to estimate the degree contact.
#
# returns: 
#
# The following list if method to estimate degree of contact is: "genetic"
#
# d_directed_ref_intervention                 # Degree of contact relative to intervention locations
# d_directed_ref_intervention_input_dataset   # Input dataset used to compute d_directed_ref_intervention
# d_directed_ref_control                      # Degree of contact relative to control locations
# d_directed_ref_control_input_dataset        # Input dataset used to compute d_directed_ref_control
#
# The following list if method to estimate degree of contact is: "distance_radius".
#
# d_disc_great_circle_ref_intervention        # Degree of contact relative to intervention locations based on great circle distances
# d_disc_great_circle_ref_control             # Degree of contact relative to control locations based on great circle distances

       

get_degree_contact <- function(method = c("genetic", "distance_radius"),
                              source_in = NULL, recipient_in = NULL, source_intervention_in = NULL, recipient_intervention_in = NULL, linked_pairs_in = NULL,
                              longitude_in = NULL, latitude_in = NULL, intervention_in = NULL, radius_in = 1, crs_in = NULL,
                              verbose_output = FALSE) {
                              
    chosen_method <- match.arg(arg = method, choices = c("genetic", "distance_radius"))                          

    if (chosen_method == "genetic") {
 
        writeLines("\n\nEstimating degree of contact with method: genetic")
           
        # Estimate degree of contact with directed transmission events
        output <- get_d_genetic_directed(source_in, recipient_in, source_intervention_in, recipient_intervention_in, linked_pairs_in, verbose_output)
 
     } else if (chosen_method == "distance_radius") {

        writeLines("\n\nEstimating degree of contact with method: distance_radius")    
        
        # Estimate degree of contact with distance radius
        
        # Great circle distance
        output <- get_d_disc_great_circle(longitude_in, latitude_in, intervention_in, radius_in, crs_in, verbose_output)
    
    } else {
    
        stop("Method to estimate degree of contact should be one of: 'genetic' or 'distance_radius'.")
    
    }
        
    output    

}



# Function: get_d_genetic_directed  ---------------------------
#
# Goal: Estimate degree of contact with directed transmission events
#
# parameters:
#
# source_in                 (character)    A vector of source location names. 
# recipient_in              (character)    A vector of recipient location names.
# source_intervention_in    (numeric)      A vector of trial arm assignments for source locations (1 = intervention, 0 = control).
# recipient_intervention_in (numeric)      A vector of trial arm assignments for recipient locations (1 = intervention, 0 = control).
# linked_pairs_in           (numeric)      A vector of counts of directed transmission events between source and recipient locations.
# verbose_output            (boolean)      A value to produce detailed output (Defualt is FALSE).
#
# returns: a list containing:
#
# d_directed_ref_intervention                 # Degree of contact relative to intervention locations
# d_directed_ref_intervention_input_dataset   # Input dataset used to compute d_directed_ref_intervention
# d_directed_ref_control                      # Degree of contact relative to control locations
# d_directed_ref_control_input_dataset        # Input dataset used to compute d_directed_ref_control
#
# ------------------------------------------------------------------------------------

get_d_genetic_directed <- function(source_in, recipient_in, 
								   source_intervention_in, recipient_intervention_in, 
								   linked_pairs_in, verbose_output = FALSE) {

	# Assemble datasets

    # check if the number of unique sources and recipients is the same
    unique_sources <- base::unique(source_in)
    unique_recipients <- base::unique(recipient_in) 	
    
	if (length(unique_sources) == length(unique_recipients)) {

		source <- base::factor(source_in) 
		recipient <- base::factor(recipient_in)
	
	} else {
	
		source <- base::as.character(source_in) 
		recipient <- base::as.character(recipient_in)
	
	}
							   
	source_intervention <- base::as.numeric(source_intervention_in)
	recipient_intervention <- base::as.numeric(recipient_intervention_in)
	linked_pairs <- base::as.numeric(linked_pairs_in)

	# Dataset to compute degree of contact with reference to intervention locations
	df_directed_ref_intervention <- base::data.frame(source, recipient, source_intervention, recipient_intervention, linked_pairs)

	if (verbose_output) {

		# View dataset structure
		base::writeLines("\nStructure of input dataset to compute degree of contact with reference to intervention locations: \n")
		utils::str(df_directed_ref_intervention)
	
		base::writeLines("\nFirst six rows of dataset: ")
		base::print(utils::head(df_directed_ref_intervention))
	
	}


	# Dataset to compute degree of contact with reference to control locations
	df_directed_ref_control <- df_directed_ref_intervention %>%
		# Swap trial arm labels such that intervention locations are labeled as control and control locations are labeled as intervention
		dplyr::rename(source_intervention_orig = source_intervention, recipient_intervention_orig = recipient_intervention) %>%
		dplyr::mutate(
		    source_intervention = dplyr::case_when(
			    # Swap trial arm assignments for source locations
			    source_intervention_orig == 0 ~ 1,
			    source_intervention_orig == 1 ~ 0),
			recipient_intervention = dplyr::case_when(
			    # Swap trial arm assignments for recipient locations
			    recipient_intervention_orig == 0 ~ 1, 
			    recipient_intervention_orig == 1 ~ 0)
			) 

	if (verbose_output) {

		# View dataset structure
		base::writeLines("\nStructure of input dataset to compute degree of contact with reference to control locations: \n")
		utils::str(df_directed_ref_control)

		base::writeLines("\nFirst six rows of dataset: ")
		base::print(utils::head(df_directed_ref_control))

	}
						   

	# -----------------------------------           


    
    # Define function: compute the degree of contact of each location relative to locations in the reference trial arm
    
    compute_d_genetic_directed <- function(df_dataset) {
        
        df_d_genetic_directed <- df_dataset %>%
            # Exclude within location transmission events if any and subset transmission events where the source is from an intervention location
            dplyr::filter(!(source == recipient) & (source_intervention == 1)) %>% 
            # Sort by recipient location
            dplyr::arrange(recipient) %>% 
            # Compute the total transmissions into each recipient location from intervention locations
            dplyr::group_by(recipient) %>% 
            dplyr::summarise(total_inflow = sum(linked_pairs, na.rm = TRUE))
                       
        df_d_genetic_directed
        
    }
    
    # Compute degree of contact relative to intervention locations
    d_intervention <- compute_d_genetic_directed(df_directed_ref_intervention)
    
    # Compute degree of contact relative to control locations
    d_control <- compute_d_genetic_directed(df_directed_ref_control)           
															
	# -----------------------------------           


	# List of items to return
	list("d_directed_ref_intervention" = d_intervention,                             # Degree of contact relative to intervention locations
		 "d_directed_ref_intervention_input_dataset" = df_directed_ref_intervention, # Input dataset used to compute d_directed_ref_intervention
		 "d_directed_ref_control" = d_control,                                       # Degree of contact relative to control locations
		 "d_directed_ref_control_input_dataset" = df_directed_ref_control)           # Input dataset used to compute d_directed_ref_control    
								   
}



# Function: get_d_disc_great_circle  ---------------------------
#
# Goal: Estimate degree of contact based on Great circle pairwise distances
#
# parameters:
#
# longitude_in              (numeric)      A vector of longitude values for locations.
# latitude_in               (numeric)      A vector of latitude values for locations.
# intervention_in           (numeric)      A vector of trial arm assignments for locations (1 = intervention, 0 = control).
# crs_in                    (numeric)      A value to indicate the coordinate reference system for supplied latitude and longitude coordinates.                  
# radius_in                 (numeric)      A value to indicate the distance radius within which to estimate the degree contact.
# verbose_output            (boolean)      A value to produce detailed output (Defualt is FALSE).
#
# returns: a list containing:
#
# d_disc_great_circle_ref_intervention        # Degree of contact relative to intervention locations based on Great circle distances
# d_disc_great_circle_ref_control             # Degree of contact relative to control locations based on Great circle distances
#
# ------------------------------------------------------------------------------------


get_d_disc_great_circle <- function(longitude_in, latitude_in, intervention_in, radius_in = 1, crs_in = 4326, verbose_output = FALSE) {


    # Assemble dataset

	longitude <- base::as.numeric(longitude_in)
	latitude <- base::as.numeric(latitude_in)
	intervention <- base::as.integer(intervention_in)
	
	radius <- base::as.numeric(radius_in)
	
	crs <- base::as.numeric(crs_in)

    # Generate indices
	index_no <- base::seq(longitude)

	# Dataset to compute degree of contact with reference to intervention locations    
    df_disc_ref_great_circle_intervention <- base::data.frame(longitude, latitude, index_no, intervention)

	# Display dataset structure 
	if (verbose_output) {

		base::writeLines("\nStructure of input dataset to compute degree of contact with reference to intervention locations: \n")
		utils::str(df_disc_ref_great_circle_intervention)

		base::writeLines("\n\n")
		base::writeLines("\nFirst six rows of dataset: ")
		base::print(head(df_disc_ref_great_circle_intervention))
	
	}
	
    
    # Dataset to compute degree of contact with reference to control locations    
    df_disc_ref_great_circle_control <- df_disc_ref_great_circle_intervention %>%
		# Swap trial arm labels such that intervention locations are labeled as control and control locations are labeled as intervention
        dplyr::rename(intervention_orig = intervention) %>%
		dplyr::mutate(intervention = dplyr::case_when(
		    intervention_orig == 0 ~ 1, 
		    intervention_orig == 1 ~ 0))

	# Display dataset structure 
	if (verbose_output) {

		base::writeLines("\nStructure of input dataset to compute degree of contact with reference to control locations: \n")
		utils::str(df_disc_ref_great_circle_control)

		base::writeLines("\n\n")
		base::writeLines("\nFirst six rows of dataset: ")
		base::print(head(df_disc_ref_great_circle_control))
	
	}

	            
	# -----------------------------------           

    
    # Compute the disc measure

    compute_d_disc_great_circle <- function(df_coords_current_location, df_coords_locations, radius, crs, verbose_output = FALSE) {
 
    
		# Get indices for individual intervention locations
		indices_intervention_locations <- base::which(df_coords_locations$intervention == 1)

		# Get index number of the location for which the degree of contact is to be computed
		idx_current_location <- df_coords_current_location["index_no"]

        
        # Compute pairwise great circle distances between locations

        df_coords_locations_sf <- sf::st_as_sf(df_coords_locations, coords = c("longitude", "latitude"), crs = sf::st_crs(crs))

        dist_mat <- sf::st_distance(df_coords_locations_sf, df_coords_locations_sf, by_element = FALSE)

		# Subset pairwise great circle distances between intervention locations and the location of interest
		sub_dist_mat <- dist_mat[idx_current_location, indices_intervention_locations]

		# -----------------------------------           


		# Compute the degree of contact within the specified radius
		
		# Reminder: This is computed as the number of intervention locations separated from 
		#           a location of interest by a distance that falls within a set radius.
	
		indices_intervention_locations_within_radius <- base::which(as.numeric(sub_dist_mat) <= radius)
	
		# Degree of contact for a control location
		if (df_coords_locations_sf$intervention[idx_current_location] == 0) {
	
			d_disc_great_circle <- base::length(indices_intervention_locations_within_radius)
	
		# Degree of contact for an intervention location
		} else if (df_coords_locations_sf$intervention[idx_current_location] == 1) {
	
			# Note: When a location of interest is in the intervention arm the number of 
			#       intervention locations within the specified radius will be inflated by 1
			#       because the distance of an intervention location to itself is zero.
			#       To address this inflation, we subtract one from the number of intervention 
			#       locations within the specified radius.

			d_disc_great_circle <- base::length(indices_intervention_locations_within_radius) - 1
	
		}

		# -----------------------------------           

		# Display results
		if (verbose_output) {

			base::writeLines("\n\n  Degree of contact for current spatial location: ")
			base::print(base::paste("idx_current_location: ", idx_current_location, "current_location_lon: ", df_coords_current_location["longitude"], "current_location_lat: ", df_coords_current_location["latitude"], "Degree of contact: ", d_disc_great_circle))    
	
		}

        # d_disc_great_circle
        
		output <- data.frame("idx" = idx_current_location, "longitude" = df_coords_current_location["longitude"], "latitude" = df_coords_current_location["latitude"], "degree_contact" =  d_disc_great_circle)
		
		output
		
	}


	# Compute degree of contact with intervention locations as the reference trial arm
	res_compute_d_disc_great_circle_ref_intervention <- base::apply(df_disc_ref_great_circle_intervention[, c("longitude", "latitude", "index_no")], 1, compute_d_disc_great_circle, df_disc_ref_great_circle_intervention[, c("longitude", "latitude", "index_no", "intervention")], radius = radius, crs = crs, verbose_output = verbose_output, simplify = TRUE) %>%
	    dplyr::bind_rows()


	# Compute degree of contact with control locations as the reference trial arm
	res_compute_d_disc_great_circle_ref_control <- base::apply(df_disc_ref_great_circle_control[, c("longitude", "latitude", "index_no")], 1, compute_d_disc_great_circle, df_disc_ref_great_circle_control[, c("longitude", "latitude", "index_no", "intervention")], radius = radius, crs = crs, verbose_output = verbose_output, simplify = TRUE) %>%
	    dplyr::bind_rows()


	# -----------------------------------           

	# List of items to return
	list("d_disc_great_circle_ref_intervention" = res_compute_d_disc_great_circle_ref_intervention,    # Degree of contact relative to intervention locations based on Great circle distances
		 "d_disc_great_circle_ref_control" = res_compute_d_disc_great_circle_ref_control)              # Degree of contact relative to control locations based on Great circle distances



}

