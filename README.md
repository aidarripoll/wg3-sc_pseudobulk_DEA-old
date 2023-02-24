# sc-eQTLGen WG3 pipeline (II): sc- and pseudobulk-differential expression analysis 

We provide **two main scripts** to peform **differential expression analysis (DEA)** with human phenotypes (e.g., sex or age) using single-cell RNA-seq data (scRNA-seq) (i.e., 10x Genomics) at different levels:

* **[single-cell level (sc-DEA)](/scDEA_MAST_glmer.R)**: using the [MAST](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5) glmer implementation.
* **[pseudobulk level (pseudobulk-DEA)](/pseudobulkDEA_limmadream.R)**: using the [limma dream](https://academic.oup.com/bioinformatics/article/37/2/192/5878955) glmer implementation.

**Of note**: 
* This analysis is meant to be run on scRNA-seq data composed by **only one sample per donor**. In case you have more than one sample per donor (e.g., stimulated vs. non-stimulated samples from the same donor, etc...) you can run our **[subsetting script](/subset_by_metadata.R)** to select individuals with specifics characteristics (e.g., only non-stimulated samples, etc...). 

* To run these scripts you should have **successfully run** the following sc-eQTLGen consortium pipelines: **WG1**, **WG2** and **WG3 (I)** 

## Contact
If you have any questions or issues, feel free to open an issue or directly email Aida Ripoll-Cladellas (aida.ripoll@bsc.es)


-------

## Required Software
**R** >=4.1.2 version: You need to install the packages loaded in the:
* Two main DEA scripts: [sc-DEA](/scDEA_MAST_glmer.R) and [pseudobulk-DEA](/pseudobulkDEA_limmadream.R)
* Add-on subsetting script: [subsetting script](/subset_by_metadata.R)
* Additional scripts in the [scripts](/scripts/) directory

## Required Input
This section explains the input data and it’s structure to run the three scripts: [sc-DEA](/scDEA_MAST_glmer.R), [pseudobulk-DEA](/pseudobulkDEA_limmadream.R), and [subsetting script](/subset_by_metadata.R).

**Of note**: To follow better the explanations in the **Required Input** section, you can clone this repository and change your current working directory. 

```
git clone https://github.com/aidarripoll/wg3-sc_pseudobulk_DEA.git
cd wg3-sc_pseudobulk_DEA
```

### Test Data
We have provided some **testing inputs** in the **[inputs directory](inputs/)** that contains the B cells outputs (Azimuth's level 1) from WG3 (I) pipeline. 

**Of note**: These files have been anonymized and they are significantly down-sized and sub-sampled versions of the whole B cells outputs from WG3 (I). Specifically, the total number of cells is 663 from 40 donors, and the number of genes is 50 and 500 for the sc- and pseudobulk-DEA, respectively.

Here is the structure of the [testing input directory](/inputs/). This input directory (*inputs/*) should have the same structure as the WG3 (I) pipeline output directory. We will need only the following files since the other ones will be used for the eQTL calling pipeline in the WG3 (II):

**inputs/**    
```bash
|-- L1
|   |-- B.Exp.txt
|   |-- B.Qced.Normalized.SCs.Rds
|   |-- B.covariates.txt 
|-- donor_pool_stim.txt

````

### Required Data
**wg3-sc_pseudobulk_DEA/**  
```bash
|-- inputs
|   |-- L1
|   |   |-- B.Exp.txt
|   |   |-- B.Qced.Normalized.SCs.Rds
|   |   |-- B.covariates.txt 
|   |-- donor_pool_stim.txt
|-- pseudobulkDEA_covariates.tab
|-- scDEA_covariates.tab
```

#### Main input directory ([inputs/](/inputs/))
It contains a directory for each Azimuth's level, **[L1](/inputs/L1/)** or L2, with the main outputs per cell type from WG3 (I): 
* **${cell_type}.Exp.txt:** Pseudobulk gene expression matrix per donor-pool combination (PFlogPF normalization + mean on QC-filtered single-cell gene expression matrix)
* **${cell_type}.Qced.Normalized.SCs.Rds:** QC-filtered single-cell gene expression matrix
* **${cell_type}.covariates.txt:** Donor metadata (from the psam file in WG1 pipeline)

It also contains a file **[donor_pool_stim.txt](/inputs/donor_pool_stim.txt)** with the sample information associated to the donor, the pool and the stimulation condtion.

*Of note*:
* Tab separated.
* The header and the values of the columns should be the same as the ones in the [pseudobulk-](/inputs/L1/B.Exp.txt) and [sc-](/inputs/L1/Qced.Normalized.SCs.Rds)gene expression files, and in the [donor metadata file](/inputs/L1/B.covariates.txt).
* The [sample information file](/inputs/donor_pool_stim.txt) provided for testing has the following structure:

| Donor  | Pool | Stimulation  | 
| ------------- | ------------- | ------------- | 
| ID_1  | pilot3_lane1  | UT  | 
| ID_2  | pilot3_lane1  | UT  | 
| ID_3  | pilot3_lane1  | UT  | 
| ID_4  | pilot3_lane1  | UT  |
| ID_5  | pilot3_lane1  | UT  | 
| ID_6  | pilot3_lane1  | UT  |

#### DEA covariates files: sc-DEA ([scDEA_covariates.tab](/scDEA_covariates.tab)) and pseudobulk-DEA ([pseudobulkDEA_covariates.tab](/pseudobulkDEA_covariates.tab))
A priori, this file should not be modified. Each tsv file that has in the:
* 1st column (covariate): Covariates included in the model
* 2nd column (type): Fixed/random effect
* 3rd column (class): Categorical (factor) or quantitative (integer, double)

*Of note*:
* Tab separated.
* The values of the 1st column (covariate) should be the same as the ones in the [pseudobulk-](/inputs/L1/B.Exp.txt) and [sc-](/inputs/L1/Qced.Normalized.SCs.Rds)gene expression files, and in the [donor metadata file](/inputs/L1/B.covariates.txt).  
* This file must have this header. 
* The covariates files provided for testing have the following structure:

sc-DEA ([scDEA_covariates.tab](/scDEA_covariates.tab))
| covariate  | type | class  | 
| ------------- | ------------- | ------------- | 
| SEX  | fixed  | factor  | 
| age  | fixed  | integer  | 
| Donor_Pool  | random  | factor  | 

pseudobulk-DEA ([pseudobulkDEA_covariates.tab](/pseudobulkDEA_covariates.tab))
| covariate  | type | class  | 
| ------------- | ------------- | ------------- | 
| SEX  | fixed  | factor  | 
| age  | fixed  | integer  | 
| CellCount  | fixed  | integer  | 
| Donor_Pool  | random  | factor  |

### Optional Data

**wg3-sc_pseudobulk_DEA/**  
```bash
|-- donor_sample.tab
```
#### Sample information associated to the donor ([donor_sample.tab](/donor_sample.tab))
This file is only needed if you have more than one sample per donor (e.g., stimulated vs. non-stimulated samples from the same donor, etc...). In this case, you need to run the **[subsetting script](/subset_by_metadata.R)** to select individuals with specifics characteristics (e.g., only non-stimulated samples, etc...). For this analysis, we ask you to only select european individuals and non-stimulated samples.

This is a tsv file that has in the:
* 1st column (type): Donor or sample metadata information.
* 2nd column (covariate): Metadata variable. 
* 3rd column (level): Category or categories to select of the associated donor/sample-metadata variable.

*Of note*:
* Tab separated.
* The values of the 2nd column (covariate) should be the same as the ones in the [pseudobulk-](/inputs/L1/B.Exp.txt) and [sc-](/inputs/L1/Qced.Normalized.SCs.Rds)gene expression files, and in the [donor metadata file](/inputs/L1/B.covariates.txt).  
* This file must have this header. 
* The sample information file provided for testing has the following structure:

| type  | covariate | level  | 
| ------------- | ------------- | ------------- | 
| Donor  | Provided_Ancestry  | EUR  | 
| Sample  | Stimulation  | UT  | 

## Running the sc/pseudobulk-DEA
*Of note*: 

* If you have not done it yet, the first step would be to *clone** this repository and change your current working directory. 

```
git clone https://github.com/aidarripoll/wg3-sc_pseudobulk_DEA.git
cd wg3-sc_pseudobulk_DEA
```

* This sc/pseudobulk-DEA is meant to be run on scRNA-seq data composed by **only one sample per donor**. In case you have more than one sample per donor (e.g., stimulated vs. non-stimulated samples from the same donor, etc...) you should run the **[subsetting script](/subset_by_metadata.R)** to select individuals with specifics characteristics (e.g., only non-stimulated samples, etc...). For this analysis, we ask you to only select european individuals and non-stimulated samples (e.g, [donor_sample.tab](/donor_sample.tab)). 

* The **functions** called in the ** mandatory sc/pseudobulk-DEA scripts** ([sc-DEA](/scDEA_MAST_glmer.R) and in the [pseudobulk-DEA](/pseudobulkDEA_limmadream.R)), and in the **[optional subsetting script](/subset_by_metadata.R)** are defined in the [external scripts](/scripts/).

### Running the subsetting script (OPTIONAL)

**1.** Set common environmental variables:  
```
cell_level=L1
cell_type=B
donor_sample=donor_sample.tab #default
input_directory=inputs #default
output_directory=subset_by_metadata #default
```

**2.** Running the **[subsetting script](/subset_by_metadata.R)**:
```
Rscript subset_by_metadata.R -l $cell_level -c $cell_type -m $donor_sample -i $input_directory -o $output_directory
```

The output directory (**[subset_metadata/](/subset_metadata/))** contains the same files as the [original input directory from WG3 (I) pipeline](/inputs/) but only for the selected samples.

```bash
-- Provided_Ancestry.EUR
    -- Stimulation.UT
        |-- L1
        |   |-- B.Exp.txt
        |   |-- B.Qced.Normalized.SCs.Rds
        |   |-- B.covariates.txt
        |   |-- B.qtlInput.Pcs.txt
        |   |-- B.qtlInput.txt
        |-- donor_pool_stim.txt

```

In this case, you will need to redefine your `$input_directory` when running the **sc/pseudobulk-DEA scripts**.
```
input_directory=subset_metadata 
```

### Running the sc-DEA script

**1.** Set common environmental variables:  
```
cell_level=L1
cell_type=B
phenotype=SEX #or age
input_directory=inputs #default
output_directory=scDEA_MAST_glmer #default
```

*Of note:*
* Set the `$phenotype` variable to SEX and age (in separate runs).
* If you needed to run the previous section on subsetting the data, you should redefine your `$input_directory` variable and add the `$donor_sample` variable:

```
input_directory=subset_metadata 
donor_sample=donor_sample.tab 
```

**2.** Running the **[sc-DEA](/scDEA_MAST_glmer.R)**:

2.1. If you needed to run the previous section on subsetting the data:

```
Rscript scDEA_MAST_glmer.R -l $cell_level -c $cell_type -v $phenotype -m $donor_sample -i $input_directory -o $output_directory 
```

2.2. If you did not need to run the previous section on subsetting the data:

```
Rscript scDEA_MAST_glmer.R -l $cell_level -c $cell_type -v $phenotype -i $input_directory -o $output_directory
```

The output directory (**[scDEA_MAST_glmer/](/scDEA_MAST_glmer/))** has the following structure:
```bash
L1
|-- SEX
|   |-- B
|       |-- de_glmer_nagq0.degs.rds
|       |-- de_glmer_nagq0.rds
|       |-- pbmc_sca.rds
|       |-- pbmc_sca_raw.rds
|       |-- pbmc_so.rds
|-- age
    |-- B
        |-- de_glmer_nagq0.degs.rds
        |-- de_glmer_nagq0.rds
        |-- pbmc_sca.rds
        |-- pbmc_sca_raw.rds
        |-- pbmc_so.rds
```

*Of note*: If you have run 2.1, the output directory (**[scDEA_MAST_glmer/](/scDEA_MAST_glmer/)**) will contain some subdirectories with the same files. For example:

```bash
Provided_Ancestry.EUR
|-- Stimulation.UT
    |-- L1
        |-- SEX
        |   |-- B
        |       |-- de_glmer_nagq0.degs.rds
        |       |-- de_glmer_nagq0.rds
        |       |-- pbmc_sca.rds
        |       |-- pbmc_sca_raw.rds
        |       |-- pbmc_so.rds
        |-- age
            |-- B
                |-- de_glmer_nagq0.degs.rds
                |-- de_glmer_nagq0.rds
                |-- pbmc_sca.rds
                |-- pbmc_sca_raw.rds
                |-- pbmc_so.rds

```

### Running the pseudobulk-DEA script

**1.** Set common environmental variables:  
```
cell_level=L1
cell_type=B
phenotype=SEX #or age
input_directory=inputs #default
output_directory=pseudobulkDEA_limmadream #default
```

*Of note:*
* Set the `$phenotype` variable to SEX and age (in separate runs).
* If you needed to run the previous section on subsetting the data, you should redefine your `$input_directory` variable and add the `$donor_sample` variable:

```
input_directory=subset_metadata 
donor_sample=donor_sample.tab 
```

**2.** Running the **[pseudobulk-DEA](/pseudobulkDEA_limmadream.R)**:

2.1. If you needed to run the previous section on subsetting the data:

```
Rscript pseudobulkDEA_limmadream.R -l $cell_level -c $cell_type -v $phenotype -m $donor_sample -i $input_directory -o $output_directory 
```

2.2. If you did not need to run the previous section on subsetting the data:

```
Rscript pseudobulkDEA_limmadream.R -l $cell_level -c $cell_type -v $phenotype -i $input_directory -o $output_directory
```

The output directory (**[pseudobulkDEA_limmadream/](/pseudobulkDEA_limmadream/)**) has the following structure:
```bash
L1
|-- SEX
|   |-- B
|       |-- eBayes
|           |-- SEX.combinations.degs.rds
|           |-- SEX.combinations.rds
|           |-- SEXM_SEXF.rds
|           |-- SEXM_SEXF.tsv
|           |-- SEXM_SEXF.vars.rds
|-- age
    |-- B
        |-- eBayes
            |-- age.combinations.degs.rds
            |-- age.combinations.rds
            |-- age.rds
            |-- age.tsv
            |-- age.vars.rds
```

*Of note*: If you have run 2.1, the output directory (**[pseudobulkDEA_limmadream/](/pseudobulkDEA_limmadream/)**) will contain some subdirectories with the same files. For example:

```bash
Provided_Ancestry.EUR/
|-- Stimulation.UT
    |-- L1
        |-- SEX
        |   |-- B
        |       |-- eBayes
        |           |-- SEX.combinations.degs.rds
        |           |-- SEX.combinations.rds
        |           |-- SEXM_SEXF.rds
        |           |-- SEXM_SEXF.tsv
        |           |-- SEXM_SEXF.vars.rds
        |-- age
            |-- B
                |-- eBayes
                    |-- age.combinations.degs.rds
                    |-- age.combinations.rds
                    |-- age.rds
                    |-- age.tsv
                    |-- age.vars.rds

```

