#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("-l", "--cell_level"), action="store", default="L1", type='character',
              help="Cell type levels: Low resolution (L1) or high resolution (L2)"),
  make_option(c("-c", "--cell_type"), action="store", default=NULL, type='character',
              help="Cell types in the low cell level resolution (l1) or in the high cell level resolution (l2)"),
  make_option(c("-m", "--donor_sample"), action="store", default='donor_sample.tab', type='character',
              help="Donor-sample file."),
  make_option(c("-i", "--in_dir"), action="store", default='inputs', type='character',
              help="Input directory: WG2 output directory."),
  make_option(c("-o","--out_dir"), action="store", default="subset_by_metadata", type='character',
              help="Output directory.")
)
opt = parse_args(OptionParser(option_list=option_list))

################################# Set Variables ################################
 # Set working directory
cwd <- getwd()
setwd(cwd)

# Loading functions
functions.fn <- 'scripts/subset_by_metadata_functions.R'
print(paste0('Loading functions from: ', functions.fn))
source(functions.fn)
cat('\n')

utils_functions.fn <- 'scripts/utils.R'
print(paste0('Loading functions from: ', utils_functions.fn))
source(utils_functions.fn)
cat('\n')

# Cell level
# opt$cell_level <- 'L2'

# Cell type
opt$cell_type <- 'B'
# opt$cell_type <- 'CD8_T'

# Output directory
# opt$out_dir <- 'subset_by_metadata'

# Input directory
# opt$in_dir <- 'inputs'
in.dir <- paste0(opt$in_dir, '/', opt$cell_level, '/')

## Covariates
# opt$donor_sample <- 'donor_sample.test.tab'
covs.fn <- opt$donor_sample

## Metadata
### donor
donor_md.fn <- paste0(in.dir, opt$cell_type, '.covariates.txt')

### sample
sample_md.fn <- paste0(opt$in_dir, '/', 'donor_pool_stim.txt')

# Print report
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Covariates file: ', covs.fn))
print(paste0('Input directory: ', opt$in_dir))
print(paste0('Output directory: ', opt$out_dir))
cat('\n')

################################ Read input data ###############################
# Read covariates file
print(paste0('Reading donor-sample file in: ',covs.fn))
covs.df <- read.table(covs.fn, header = T)
covs.list <- lapply(rsplit(covs.df, covs.df[,c('type','covariate')]), function(x) lapply(x, function(y) y$level))
covs_levels.list <- extract_levels(covs.df)

# Read metadata files
## donor
print(paste0('Reading donor metadata file in: ',donor_md.fn))
donor.df <- read.delim(donor_md.fn, check.names = FALSE)

## sample
print(paste0('Reading sample metadata file in: ',sample_md.fn))
sample.df <- read.delim(sample_md.fn)

## metadata list
metadata_list <- list(Donor = donor.df,
                      Sample = sample.df)
print('Checking covariates...')
check_covariates.res <- sapply(names(metadata_list), function(i) check_covariates(i, metadata_list, covs_levels.list), simplify = FALSE)
cat('\n')

################################ Pre-processing ###############################
# Output directory
## get output tag
type_tag.list <- lapply(covs_levels.list, function(i) get_tag_type(i))
out_sdir <- paste(unname(unlist(type_tag.list)), collapse = '/')

## create output directory
out.dir <- paste0(opt$out_dir, '/', out_sdir, '/')
out.sdir <- paste0(out.dir, opt$cell_level, '/')
if(!dir.exists(out.sdir)){dir.create(out.sdir, recursive = T)}

####################### Get the selected samples (Donor_Pool) ##################
# Get the set of selected Donors and Pool for specific traits
## donors
get_donors.list <- sapply(names(covs_levels.list$Donor), function(i) get_donors(i, covs_levels.list$Donor, donor.df), simplify = FALSE)
get_donors.list <- Reduce("intersect", get_donors.list)

## samples
get_samples.list <- sapply(names(covs_levels.list$Sample), function(i) get_samples(i, covs_levels.list$Sample, sample.df), simplify = FALSE)
get_samples.list <- Reduce("intersect", get_samples.list)

## combine donors and samples
donors_samples.list <- list(donors = get_donors.list, 
                            samples = get_samples.list)
donors_samples.combined <- Reduce("intersect", donors_samples.list)
donors_samples.tag <- paste(donors_samples.combined, collapse = ', ')

# Get the subset for the Donor_Pool_Stimulation file
print('Subsetting Donor_Pool_Stimulation sample file...')
sample_selected.df <- sample.df
sample_selected.df$Donor_Pool <- paste(sample_selected.df$Donor, sample_selected.df$Pool, sep = ';;')
sample_selected.df <- sample_selected.df[sample_selected.df$Donor_Pool%in%donors_samples.combined,]
check_samples <- all(donors_samples.combined%in%unique(sample_selected.df$Donor_Pool))
print(paste0('All samples: ', check_samples))
sample_selected.df <- sample_selected.df[,c('Donor', 'Pool', 'Stimulation')]
sample_selected.fn <- paste0(out.dir, 'donor_pool_stim.txt')
print(paste0('Saving new file in: ', sample_selected.fn))
write.table(sample_selected.df, sample_selected.fn, sep = '\t', row.names = FALSE, quote = FALSE)
cat('\n')

######################### Subset sc/pseudobulk-DEA inputs #####################
print('Subsetting sc/pseudobulk-DEA input files...')
print(paste0('# of samples: ', length(donors_samples.combined)))
print(paste0('Sample IDs: ', donors_samples.tag))
cat('\n')

# Covariates metadata
print(paste0('Subsetting Metadata file in: ', donor_md.fn))
donor_selected.df <- droplevels(donor.df[donor.df$Donor_Pool%in%donors_samples.combined,])
check_samples <- all(unique(donor_selected.df$Donor_Pool)%in%donors_samples.combined)
print(paste0('All samples: ', check_samples))
donor_selected.fn <- paste0(out.sdir, opt$cell_type, '.covariates.txt')
print(paste0('Saving new file in: ', donor_selected.fn))
write.table(donor_selected.df, donor_selected.fn, sep = '\t', row.names = FALSE, quote = FALSE)
cat('\n')

# Pseudobulk gene expression matrix (mean+PFlogPF)
mean_mat.fn <- paste0(in.dir, '/', opt$cell_type, '.Exp.txt')
print(paste0('Reading and Subsetting Pseudobulk gene expression (mean+PFlogPF) file in: ',mean_mat.fn))
system.time(mean_mat <- read.delim(mean_mat.fn, check.names = FALSE))
samples.idx <- which(colnames(mean_mat)%in%donors_samples.combined)
mean_mat_selected <- mean_mat[,samples.idx]
mean_mat_selected <- cbind(mean_mat[,1],mean_mat_selected)
colnames(mean_mat_selected)[1] <- ''
samples <- colnames(mean_mat_selected)[-which(colnames(mean_mat_selected)=='')]
check_samples <- all(samples%in%donors_samples.combined)
print(paste0('All samples: ', check_samples))
mean_mat_selected.fn <- paste0(out.sdir, opt$cell_type, '.Exp.txt')
print(paste0('Saving new file in: ', mean_mat_selected.fn))
write.table(mean_mat_selected, mean_mat_selected.fn, sep = '\t', row.names = FALSE, quote = FALSE)
cat('\n')

# Single-cell gene expression matrix
so.fn <- paste0(in.dir, '/', opt$cell_type, '.Qced.Normalized.SCs.Rds')
print(paste0('Single-cell gene expression file in: ',so.fn))
system.time(so <- readRDS(so.fn))
cells_selected <- rownames(so@meta.data[so@meta.data$Donor_Pool%in%donors_samples.combined,])
so_selected <- so[,cells_selected]
check_samples <- all(unique(so@meta.data$Donor_Pool)%in%donors_samples.combined)
print(paste0('All samples: ', check_samples))
so_selected.fn <- paste0(out.sdir, opt$cell_type,'.Qced.Normalized.SCs.Rds')
print(paste0('Saving new file in: ', so_selected.fn))
saveRDS(so_selected, so_selected.fn)

# qtlInput
qtlInput.fn <- paste0(in.dir, '/', opt$cell_type, '.qtlInput.txt')
print(paste0('Reading and Subsetting qtlInput file in: ',qtlInput.fn))
system.time(qtlInput <- read.delim(qtlInput.fn, check.names = FALSE))
samples.idx <- which(colnames(qtlInput)%in%donors_samples.combined)
qtlInput_selected <- qtlInput[,samples.idx]
qtlInput_selected <- cbind(qtlInput[,1],qtlInput_selected)
colnames(qtlInput_selected)[1] <- ''
samples <- colnames(qtlInput_selected)[-which(colnames(qtlInput_selected)=='')]
check_samples <- all(samples%in%donors_samples.combined)
print(paste0('All samples: ', check_samples))
qtlInput_selected.fn <- paste0(out.sdir, opt$cell_type, '.qtlInput.txt')
print(paste0('Saving new file in: ', qtlInput_selected.fn))
write.table(qtlInput_selected, qtlInput_selected.fn, sep = '\t', row.names = FALSE, quote = FALSE)
cat('\n')

# qtlInput Pcs
qtlInput_Pcs.fn <- paste0(in.dir, '/', opt$cell_type, '.qtlInput.Pcs.txt')
print(paste0('Reading and Subsetting qtlInput Pcs file in: ',qtlInput_Pcs.fn))
system.time(qtlInput_Pcs <- read.delim(qtlInput_Pcs.fn, check.names = FALSE))
samples.idx <- which(qtlInput_Pcs[[1]]%in%donors_samples.combined)
qtlInput_Pcs_selected <- qtlInput_Pcs[samples.idx,]
samples <- qtlInput_Pcs_selected[[1]]
check_samples <- all(samples%in%donors_samples.combined)
print(paste0('All samples: ', check_samples))
qtlInput_Pcs_selected.fn <- paste0(out.sdir, opt$cell_type, '.qtlInput.Pcs.txt')
print(paste0('Saving new file in: ', qtlInput_Pcs_selected.fn))
write.table(qtlInput_Pcs_selected, qtlInput_Pcs_selected.fn, sep = '\t', row.names = FALSE, quote = FALSE)
cat('\n')