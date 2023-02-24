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
shhh(library(collapse))

######################## Functions used in several scripts #####################
# 1. Convert covariates levels (if it is a number); accessory of extract_levels()
# df <- covs.list[[1]][[2]]
# l = level_var
convert_levels <- function(df, l){
  x <- df[[l]]
  x <- unlist(str_split(x,','))
  x_check <- all(grepl("[[:digit:]]",x))
  if(x_check){
    x <- as.numeric(x)
  }
  return(x)
}

# 2. Extract covariates levels
# df <- covs.df
# cnames = c('type', 'covariate')
# level_var = 'level'
extract_levels <- function(df, cnames = c('type', 'covariate'), level_var = 'level'){
  covs.list <- rsplit(df, df[,cnames])
  covs_converted.list <- lapply(covs.list, function(l) lapply(l, function(x) convert_levels(x, level_var)))
  return(covs_converted.list)
}

# 3. Get tags for the output directory; accessory of get_tag_type()
# i <- c
# i_list <- type_list
get_tag <- function(i, i_list){
  # print(i)
  l <- i_list[[i]]
  l_tag <- paste(l, collapse = '_')
  c_l_tag <- paste(i, l_tag, sep = '.')
  return(c_l_tag)
}

# 4. Get tags for the different covariate types for the output directory
# l <- covs_levels.list[[1]]
get_tag_type <- function(l){
  covariates <- names(l)
  covariates_list <- sapply(covariates, function(c) get_tag(c, l), simplify = FALSE)
  type_l <- unname(unlist(covariates_list))
  type_tag <- paste(type_l, collapse = '-')
  return(type_tag)
}