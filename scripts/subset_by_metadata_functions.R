#!/usr/bin/env Rscript

############################### Load R packages ################################
print('Loading R packages...')
shhh <- suppressPackageStartupMessages
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(tibble))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(Seurat))
shhh(library(SeuratObject))

######################## Functions used in subset_by_metadata.R ##################
# 1. Check covariates are in the donor/sample metadata files
# md_type <- names(metadata_list)[1]
# md_list <- metadata_list
# covs_list <- covs_levels.list
check_covariates <- function(md_type, md_list, covs_list){
  # print(md_type)
  cov_vars <- names(covs_list[[md_type]])
  metadata_vars <- colnames(md_list[[md_type]])
  if(!all(cov_vars%in%metadata_vars)){
    err_message <- paste0('Your ', md_type, ' metadata variables (', paste(cov_vars, collapse = ', '), ') should be in the metadata file.')
    stop(err_message)
  }else{
    print(paste0('Your ', md_type, ' metadata variables (', paste(cov_vars, collapse = ', '), ') are in the metadata file.'))
  }
  return(NULL)
}

# 2. Get selected donors
# donor_cov <- names(covs_levels.list$Donor)[1]
# donor_list <- covs_levels.list$Donor
# donor_df <- donor.df
get_donors <- function(donor_cov, donor_list, donor_df){
  # print(donor_cov)
  l <- donor_list[[donor_cov]]
  donor_df.l <- donor_df[donor_df[[donor_cov]]%in%l,]
  donors <- donor_df.l$Donor_Pool
  return(donors)
}

# 3. Get selected samples
# sample_cov <- names(covs_levels.list$Sample)[1]
# sample_list <- covs_levels.list$Sample
# sample_df <- sample.df
get_samples <- function(sample_cov, sample_list, sample_df){
  # print(sample_cov)
  l <- sample_list[[sample_cov]]
  sample_df.l <- sample_df[sample_df[[sample_cov]]%in%l,]
  sample_df.l$Donor_Pool <- paste(sample_df.l$Donor, sample_df.l$Pool, sep = ';;')
  samples <- sample_df.l$Donor_Pool
  return(samples)
}