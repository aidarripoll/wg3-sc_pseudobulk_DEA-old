#!/usr/bin/env Rscript

############################### Load R packages ################################
print('Loading R packages...')
shhh <- suppressPackageStartupMessages
shhh(library(edgeR))
shhh(library(limma))
shhh(library(variancePartition))
shhh(library(BiocParallel))
shhh(library(data.table))
shhh(library(lme4))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(tibble))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))

######################## Functions used in pseudobulkDEA_limmadream.R ##################
# 1. Modify metadata
# metadata <- metadata_df
# phenotype <- opt$phenotype
# covs_df <- covs.df
modify_metadata_pseudo <- function(metadata, phenotype, covs_df){
  # Create random factors variables in the metadata
  random_vars <- covs_df[covs_df$type=='random',]$covariate
  split_len <- length(unlist(str_split(random_vars, '_')))
  random_names <- str_split_fixed(random_vars, '_', split_len)[1,]
  metadata[,random_names] <- str_split_fixed(metadata[[random_vars]], ';;', split_len)[,c(1:split_len)]
  aggregates <- random_names[1]
  vars_model <- covs.df[covs.df$type=='fixed',]$covariate
  vars_model <- c(vars_model, random_names[-1])
  md_vars <- c(aggregates, vars_model)
  
  # Phenotype modifications
  if(phenotype=='SEX'){
    print(paste0('Relevel ', phenotype, '...'))
    metadata[[phenotype]] <-  ifelse(metadata[[phenotype]]==1,'M','F')
  }
  
  if(phenotype=='age_cat'){
    print(paste0('Creating a new metadata variable: ', phenotype, '...'))
    metadata[[phenotype]] <- ifelse(metadata$age<=40, 'Y',
                                    ifelse(metadata$age>=60, 'O', 'M')) 
    metadata <- metadata[metadata[[phenotype]]%in%c('Y','O'),] #only keep Y and O samples (remove if we want to consider all samples)
  }
  
  if(phenotype=='age_cat_all'){
    print(paste0('Creating a new metadata variable: ', phenotype, '...'))
    metadata[[phenotype]] <- ifelse(metadata$age<=40, 'Y', 'O')
  }
  
  if(phenotype=='age_squared'){
    print(paste0('Creating a new metadata variable: ', phenotype, '...'))
    metadata[[phenotype]] <- metadata$age^2
  }
  
  if(phenotype%in%c('SEX','age_cat', 'age_cat_all')){
    print(paste0('Order ', phenotype, ' factor levels...'))
    Group_order <- c('M','F')
    if(phenotype%in%c('age_cat', 'age_cat_all')){
      Group_order <- c('Y','O')
      # # add if we want to consider all samples
      # if(phenotype%in%c('age_cat')){Group_order <- c('Y','M','O')}
      # if(phenotype%in%c('age_cat_all')){Group_order <- c('Y','O')}
    }
    metadata[[phenotype]] <- factor(metadata[[phenotype]],
                                    levels = Group_order)
    phenotype_table <- sort(table(metadata[[phenotype]]),decreasing=TRUE)
    print(phenotype_table)
  }else{
    print(summary(metadata[[phenotype]]))
  }
  
  # Declare factors
  factor_vars <- covs_df[covs_df$class=='factor' & covs_df$type=='fixed',]$covariate
  factor_vars <- c(factor_vars, random_names)
  if(phenotype%in%c('age_cat','age_cat_all','age_squared')){
    selected_vars <- c(factor_vars, phenotype)
    if(phenotype%in%c('age_cat','age_cat_all')){
      factor_vars <- c(factor_vars, phenotype)
    }
  }else{
    selected_vars <- factor_vars
  }
  
  # Pick selected columns
  selected_vars <- union(md_vars, selected_vars)
  aggregate_metadata <- droplevels(unique(metadata[,selected_vars]))
  
  # Make factors
  aggregate_metadata %>% 
    mutate_at(all_of(factor_vars), as.factor) %>%
    as.data.frame() -> aggregate_metadata
  aggregate_metadata %>% 
    mutate_at(all_of(factor_vars), droplevels) %>%
    as.data.frame() -> aggregate_metadata
  
  # Reformat
  selected_vars <- selected_vars[!selected_vars%in%aggregates]
  rownames(aggregate_metadata) <- as.character(unique(aggregate_metadata[[aggregates]]))
  aggregate_metadata <- aggregate_metadata[,selected_vars]
  
  return(aggregate_metadata)
}

# 2. Modify gene expression matrix 
# mean_mat <- mean_pseudo.mat
modify_ge <- function(mean_mat){
  # colnames
  donors <- str_split_fixed(colnames(mean_mat), ';;', 2)[,1]
  colnames(mean_mat) <- donors
  
  # rownames
  rownames(mean_mat) <- mean_mat[,1]
  mean_mat <- mean_mat[,-1]
  
  return(mean_mat)
}

# 3. Subset gene expression matrix (if needed --> if opt$phenotype==age_cat, keep only Y and O samples)
# mean_mat <- mean_pseudo
# md <- metadata
subset_ge <- function(mean_mat, md){
  donors <- rownames(md)
  mean_mat.subset <- mean_mat[,colnames(mean_mat)%in%donors]
  return(mean_mat.subset)
}

# 4. Pick highly expressed genes - 50% expressed
## 4.1 Calculate gene expression frequency --> Mymic MAST::freq(), but faster (accessory)
# mat <- mean_mat
# na_rm=TRUE
freq_expr <- function(mat, na_rm){
  out <- apply(t(mat)>0, 2, mean, na.rm=na_rm)
  return(out)
}

## 4.2 Filter by expression frequency threshold (main)
# mean_mat <- mean_pseudo
# expr_th <- opt$expr
highly_expressed <- function(mean_mat, expr_th, na_rm = TRUE){
  n_genes <- nrow(mean_mat)
  
  # Calculate freq of expression by gene
  freq_expr_gene <- freq_expr(mean_mat, na_rm)
  
  # Set a minimum freq threshold
  ## select genes with at least some expression in some donor (filter out genes with no expression in any donor)
  genes.expr <- freq_expr_gene>0 
  pseudo_counts <- mean_mat[genes.expr,] 
  n_genes.expr <- nrow(pseudo_counts)
  
  ## select genes with a minimum of frequency across donors
  genes.expr_th <- freq_expr(pseudo_counts, na_rm) > expr_th
  pseudo_counts <- pseudo_counts[genes.expr_th,]
  n_genes.expr_th <- nrow(pseudo_counts)
  
  # Pick expressed genes
  genes_expr <- rownames(pseudo_counts)
  
  # Report
  print(paste0('# of initial genes: ', n_genes))
  print(paste0('# of genes with at least some expression in some donor: ', n_genes.expr))
  print(paste0('# of genes with a minimum of frequency across donors: ', n_genes.expr_th))
  
  return(pseudo_counts)
}

# 5. Pick highly variable genes - Highly variable genes (HVGs) were defined as the genes in the top two quartiles based on their squared coefficient of variation (CV^2 = variance / mean^2) calculated across all cells of each different cell-type.
## 5.1 Calculate CV (acessory)
# mat <- mean_mat
# na_rm=TRUE
cv_squared <- function(x){
  sd(x)/mean(x)^2
}
func_to_list <- function(l, f){
  vapply(l, f, FUN.VALUE=0.0)
}

## 5.2 Filter by CV threshold (main)
# mean_mat <- mean_pseudo.expr
# margin.idx = 1
# cv_th <- opt$cv
highly_variable <- function(mean_mat, cv_th, margin.idx = 1){
  n_genes <- nrow(mean_mat)
  
  # Split matrix into list by rows
  print('Splitting matrix to list by rows...')
  system.time(gene_vec.list <- asplit(mean_mat, margin.idx))
  
  # Apply CV function to the list of gene vectors
  print('Applying CV function to the list of gene vectors...')
  system.time(cv_squared_gene <- func_to_list(gene_vec.list, cv_squared))
  cv_squared_gene <- sort(cv_squared_gene, decreasing=T)
  
  # Set a minimum cv-squared threshold
  ## select genes with at least some variation in some donor (filter out genes with no variation in any donor)
  genes.var <- cv_squared_gene>0 
  pseudo_counts.var <- mean_mat[genes.var,] 
  n_genes.var <- nrow(pseudo_counts.var)
  
  ## select genes with a minimum of variation across donors
  ### calculate the value of the cv-squared quantile 
  cv_th.value <- unname(quantile(cv_squared_gene, probs = cv_th, na.rm = T)) 
  print(paste0('Quantile (probs): ', as.character(round(cv_th,3)), 
               ' --> value: ', as.character(round(cv_th.value,3))))
  
  ### select the highly variable genes
  genes.hvgs_th <- cv_squared_gene>cv_th.value
  pseudo_counts <- mean_mat[genes.hvgs_th,]
  n_genes.hvgs_th <- nrow(pseudo_counts)
  
  # Pick expressed genes
  hvgs <- rownames(pseudo_counts)
  
  # Report
  print(paste0('# of initial genes: ', n_genes))
  print(paste0('# of genes with at least some variation: ', n_genes.var))
  print(paste0('# of genes with a minimum of variation: ', n_genes.hvgs_th))
  
  return(pseudo_counts)
}

# 6. Make sure the donor ids from the metadata (aggregate_metadata) matches the ones in the expression data (geneExpr)
# mean_mat <- mean_pseudo.expr.hvgs
# md <- metadata
check_donor_ids <- function(mean_mat, md){
  if(!identical(colnames(mean_mat), rownames(md))){
    print('Matching donor ids from the metadata and gene expression...')
    idx <- match(colnames(mean_mat), rownames(md))
    md <- md[idx,]
  }
  return(md)
}

# 7. getContrast_by_comb for 'dreamer' function
# comb <- combinations[[1]]
# comb_param <- comb_parameter
# outdir <- out_dir
# form_model <- form
# vobjDream_obj <- vobjDream
# metadata_obj <- metadata
# contrast <- contrast
# eBayes_var = eBayes_var
getContrast_by_comb <- function(comb = NULL, comb_param = NULL, outdir, form_model, vobjDream_obj, metadata_obj, contrast, eBayes_var){
  # get contrast
  if(is.null(comb)){
    print('Continous variable...')
    L = getContrast(vobjDream_obj, form_model, metadata_obj, contrast)
    comb <- contrast
  }else{
    print('Categorical variable...')
    print(comb)
    metadata_obj[[contrast]] <- as.factor(metadata_obj[[contrast]])
    if(comb_param){
      print('The contrast variable has >2 levels...') #not tested
      coef_var <- rev(c(comb[[1]], comb[[2]]))
      L = tryCatch(
        {
          getContrast(vobjDream_obj, form_model, metadata_obj, coef_var)
        }, 
        error = function(e){
          getContrast(vobjDream_obj, form_model, metadata_obj, coef_var[[1]])
        })
    }else{
      print('The contrast variable has <=2 levels...')
      comb_levels <- gsub(contrast, '', comb)
      metadata_obj[[contrast]] <- factor(metadata_obj[[contrast]], 
                                         levels=comb_levels)
      coef_var <- comb[[2]]
      L = getContrast(vobjDream_obj, form_model, metadata_obj, coef_var)
    }
  }
  
  # dream()
  print('Dream() default: using the whole object...')
  
  ## fit contrast
  fit = dream(vobjDream_obj, form_model, metadata_obj, L)
  if(eBayes_var){
    fit = eBayes(fit)
  }
  
  ## grab the exact fit
  limma_res <- variancePartition::topTable(fit, coef=c('L1'), number=length(fit$F.p.value))
  limma_res <- limma_res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    arrange(adj.P.Val)
  limma_res <- limma_res[order(limma_res$P.Value),]
  
  # output: create list of variables
  vars_list <- list()
  vars_list[['form']] <- form_model
  vars_list[['vobjDream']] <- vobjDream_obj
  vars_list[['metadata']] <- metadata_obj
  vars_list[['L']] <- L
  vars_list[['fit']] <- fit 
  
  # save/write the result
  comb_tag <- paste(comb,collapse='_')
  limma_res.fn <- paste0(outdir, comb_tag)
  limma_res_tsv.fn <- paste0(limma_res.fn, '.tsv')
  print('Saving pseudobulk-DEA results with limma dream in: ')
  print(limma_res_tsv.fn)
  write.table(limma_res, limma_res_tsv.fn, sep = '\t', row.names = T)
  limma_res_rds.fn <- paste0(limma_res.fn, '.rds')
  print(limma_res_rds.fn)
  saveRDS(limma_res, limma_res_rds.fn)
  vars.fn <- paste0(limma_res.fn, '.vars.rds')
  print(vars.fn)
  saveRDS(vars_list, vars.fn)
  
  return(limma_res)
}

# 8. Limma dream function
# ge_dge = geneExpr
# contrast = opt$phenotype
# fixed_effects = fixed_effects.vars
# random_effects = random_effects.vars
# metadata = aggregate_metadata
# eBayes_var = opt$eBayes
# span_var = opt$span
# weights = opt$weights
# out_dir = out.dir
# counts = 'CellCount'
dreamer <- function(ge_dge, contrast, fixed_effects, random_effects, metadata, eBayes_var, span_var, weights, out_dir, counts = 'CellCount'){
  n_genes <- nrow(ge_dge)
  n_samples <- ncol(ge_dge)
  print(paste0('# of initial genes: ', n_genes))
  print(paste0('# of samples: ', n_samples))
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(4, "SOCK", progressbar=TRUE)
  register(param)
  
  # The variable to be tested must be a fixed effect
  ## defining the form
  contrast_fixed.fmla <- paste(c(contrast,fixed_effects),collapse='+')
  form_vars <- contrast_fixed.fmla
  contrast.var_out <- contrast
  if(!is.null(random_effects)){
    random_effects.fmla <- paste(paste0('(1|',random_effects,')'),collapse='+')
    form_vars <- paste(c(contrast_fixed.fmla,random_effects.fmla), collapse='+')
  }
  form_vars <- paste0('~',form_vars)
  form <- as.formula(form_vars)
  print(paste0('Fitting lmer: ',form_vars))
  print(paste0('Testing: ',contrast))
  
  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights(ge_dge, form, metadata, span = span_var)
  vobjDream.max <- apply(vobjDream$weights, 1, max)
  names(vobjDream.max) <- rownames(ge_dge)
  print('Check voomWithDreamWeights() -weights- maximum distribution:')
  print(summary(vobjDream.max))
  print('Genes with highest weights:')
  print(head(sort(vobjDream.max,decreasing=TRUE)))
  
  if(!is.null(weights)){
    print(paste0('Filtering genes with very high weights (top ', weights, ')...'))
    weights_probs <- 1-weights
    weights_th.value <- unname(quantile(vobjDream.max, probs = weights_probs, na.rm = T)) 
    print(paste0('Quantile (probs): ', as.character(weights_probs), 
                 ' --> value: ', as.character(round(weights_th.value,3))))
    genes.weights_th <- vobjDream.max<weights_th.value
    n_genes.removed <- unname(table(genes.weights_th)['FALSE'])
    # which(genes.weights_th==FALSE)
    
    print(paste0('Filtering gene expression matrix...'))
    ge_dge <- ge_dge[genes.weights_th,]
    
    print(paste0('Recomputing voomWithDreamWeights()...'))
    vobjDream = voomWithDreamWeights(ge_dge, form, metadata, span = span_var)
    vobjDream.max <- apply(vobjDream$weights, 1, max)
    names(vobjDream.max) <- rownames(ge_dge)
    print('Check voomWithDreamWeights() maximum distribution:')
    print(summary(vobjDream.max))
    print('Genes with highest weights:')
    print(head(sort(vobjDream.max,decreasing=TRUE)))
    
    # Report
    print(paste0('# of genes with high weights (filtered out): ', n_genes.removed))
    n_genes <- n_genes-n_genes.removed
  }
  print(paste0('# of tested genes: ', n_genes))
  
  # Fit the dream model on each gene
  # By default, uses the Satterthwaite approximation for the hypothesis test
  #fit = dream(vobjDream, form, metadata)
  
  ## use a contrast matrix
  combinations <- NULL
  if(is.factor(metadata[[contrast]])){
    contrast_levels <- levels(metadata[[contrast]])
    combinations <- split(combn(contrast_levels,2),  col(combn(contrast_levels,2)))
    combinations <- lapply(combinations, function(x) paste0(contrast,x))
    combinations <- unname(combinations)
    names(combinations) <- unlist(lapply(combinations, function(x) paste(x,collapse='_')))
    comb_parameter <- ifelse(nlevels(metadata[[contrast]])>2, TRUE, FALSE)
    # testing
    # combinations <- c(combinations, list(c("SexF", "SexM")))
  }
  
  ## apply getContrast_by_comb function
  if(!is.null(combinations)){
    getContrast.list <- lapply(combinations, function(i) getContrast_by_comb(comb = i,
                                                                             comb_param = comb_parameter,
                                                                             outdir = out_dir,
                                                                             form_model = form,
                                                                             vobjDream_obj = vobjDream,
                                                                             metadata_obj = metadata, 
                                                                             contrast = contrast, 
                                                                             eBayes_var = eBayes_var))
  }else{
    getContrast <- getContrast_by_comb(outdir = out_dir,
                                       form_model = form,
                                       vobjDream_obj = vobjDream,
                                       metadata_obj = metadata, 
                                       contrast = contrast,
                                       eBayes_var = eBayes_var)
    getContrast.list <- list()
    getContrast.list[[contrast]] <- getContrast
  }
  
  ## save
  suffix_fn <- paste0(contrast, '.combinations')
  getContrast.list_fn <- paste0(out_dir, suffix_fn, '.rds')
  print(paste0('Saving all combinations pseudobulk-DEA results in: ',getContrast.list_fn))
  saveRDS(getContrast.list, getContrast.list_fn)
  return(getContrast.list)
}

# 9. Get pseudobulk-DEGs from limma dream output
# comb <- names(dreamer.res)[1]
# dea_list <- dreamer.res
get_degs <- function(comb, dea_list){
  df <- dea_list[[comb]]
  degs.df <- df[df$adj.P.Val<=0.05,]
  return(degs.df)
}