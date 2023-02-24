# sc-eQTLGen WG3 pipeline (II): sc- and pseudobulk-differential expression analysis 

We provide **two main scripts** to peform **differential expression analysis (DEA)** with human phenotypes (e.g., sex or age) using single-cell RNA-seq data (scRNA-seq) (i.e., 10x Genomics) at different levels:

* **[single-cell level (sc-DEA)](scDEA_MAST_glmer.R)**: using the [MAST](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5) glmer implementation.
* **[pseudobulk level (pseudobulk-DEA)](pseudobulkDEA_limmadream.R)**: using the [limma dream](https://academic.oup.com/bioinformatics/article/37/2/192/5878955) glmer implementation.

**Of note**: 
* This analysis is meant to be run on scRNA-seq data composed by **only one sample per donor**. In case you have more than one sample per donor (e.g., stimulated vs. non-stimulated samples from the same donor, etc...) you can run our **[subsetting script](subset_by_metadata.R)** to select individuals with specifics characteristics (e.g., only non-stimulated samples, etc...). 

* To run these scripts you should have **successfully run** the following sc-eQTLGen consortium pipelines: **WG1**, **WG2** and **WG3 (I)** 

## Contact
If you have any questions or issues, feel free to open an issue or directly email Aida Ripoll-Cladellas (aida.ripoll@bsc.es)

-------

## Required Software
**R** >=4.1.2 version: You need to install the packages loaded in the:
* Two main DEA scripts: [sc-DEA](scDEA_MAST_glmer.R) and [pseudobulk-DEA](pseudobulkDEA_limmadream.R)
* Add-on subsetting script: [subsetting script](subset_by_metadata.R)
* Additional scripts in the [scripts](scripts/) directory

## Required Input
This section explains the input data and it’s structure to run the three scripts: [sc-DEA](scDEA_MAST_glmer.R), [pseudobulk-DEA](pseudobulkDEA_limmadream.R), and [subsetting script](subset_by_metadata.R).

**Of note**: To follow better the explanations in the **Required Input** section, you can clone this repository and change your current working directory. 

```
git clone https://github.com/aidarripoll/wg3-sc_pseudobulk_DEA.git
cd wg3-sc_pseudobulk_DEA
```

### Test Data
We have provided some **testing inputs** in the **[inputs directory](inputs/)** that contains the B cells outputs (Azimuth's level 1) from WG3 (I) pipeline. 

**Of note**: These files have been anonymized and they are significantly down-sized and sub-sampled versions of the whole B cells outputs from WG3 (I). Specifically, the total number of cells is 663 from 40 donors, and the number of genes is 50 and 500 for the sc- and pseudobulk-DEA, respectively.

Here is the structure of the [testing input directory](inputs/). This input directory (*inputs/*) should have the same structure as the WG3 (I) pipeline output directory. We will need only the following files since the other ones will be used for the eQTL calling pipeline in the WG3 (II):

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
|-- donor_sample.tab
|-- inputs
|   |-- L1
|   |   |-- B.Exp.txt
|   |   |-- B.Qced.Normalized.SCs.Rds
|   |   |-- B.covariates.txt 
|   |-- donor_pool_stim.txt
|-- pseudobulkDEA_covariates.tab
|-- scDEA_covariates.tab
```

#### Main input directory ([inputs/](inputs/))
It contains a directory for each Azimuth's level, **[L1](inputs/L1/)** or L2, with the main outputs per cell type from WG3 (I): 
* **${cell_type}.Exp.txt:** Pseudobulk gene expression matrix per donor-pool combination (PFlogPF normalization + mean on QC-filtered single-cell gene expression matrix)
* **${cell_type}.Qced.Normalized.SCs.Rds:** QC-filtered single-cell gene expression matrix
* **${cell_type}.covariates.txt:** Donor metadata (from the psam file in WG1 pipeline)

It also contains a file **[donor_pool_stim.txt](inputs/donor_pool_stim.txt)** with the sample information associated to the donor, the pool and the stimulation condtion.

*Of note*:
* Tab separated.
* The header and the values of the columns should be the same as the ones in the [pseudobulk-](inputs/L1/B.Exp.txt) and [sc-](inputs/L1/Qced.Normalized.SCs.Rds)gene expression files, and in the [donor metadata file](inputs/L1/B.covariates.txt).
* The [sample information file](inputs/donor_pool_stim.txt) provided for testing has the following structure:

| Donor  | Pool | Stimulation  | 
| ------------- | ------------- | ------------- | 
| ID_1  | pilot3_lane1  | UT  | 
| ID_2  | pilot3_lane1  | UT  | 
| ID_3  | pilot3_lane1  | UT  | 
| ID_4  | pilot3_lane1  | UT  |
| ID_5  | pilot3_lane1  | UT  | 
| ID_6  | pilot3_lane1  | UT  |

#### DEA covariates files: sc-DEA ([scDEA_covariates.tab](scDEA_covariates.tab)) and pseudobulk-DEA ([pseudobulkDEA_covariates.tab](pseudobulkDEA_covariates.tab))
A priori, this file should not be modified. Each tsv file that has in the:
* 1st column (covariate): Covariates included in the model
* 2nd column (type): Fixed/random effect
* 3rd column (class): Categorical (factor) or quantitative (integer, double)

*Of note*:
* Tab separated.
* The values of the 1st column (covariate) should be the same as the ones in the [pseudobulk-](inputs/L1/B.Exp.txt) and [sc-](inputs/L1/Qced.Normalized.SCs.Rds)gene expression files, and in the [donor metadata file](inputs/L1/B.covariates.txt).  
* This file must have this header. 
* The covariates files provided for testing have the following structure:

sc-DEA ([scDEA_covariates.tab](scDEA_covariates.tab))
| covariate  | type | class  | 
| ------------- | ------------- | ------------- | 
| SEX  | fixed  | factor  | 
| age  | fixed  | integer  | 
| Donor_Pool  | random  | factor  | 

pseudobulk-DEA ([pseudobulkDEA_covariates.tab](pseudobulkDEA_covariates.tab))
| covariate  | type | class  | 
| ------------- | ------------- | ------------- | 
| SEX  | fixed  | factor  | 
| age  | fixed  | integer  | 
| CellCount  | fixed  | integer  | 
| Donor_Pool  | random  | factor  |

### Optional Data
#### Metadata variables ([metadata_variables.tab](/metadata_variables.tab))
A tsv file that has in the:
* 1st column: Metadata variable name. 
* 2nd column: Metadata variable type. 
* 3rd and 4rd columns: minimum and maximum MADs. By default, *minimum*=1 and *maximum*=5.

*Of note*:
* Tab separated.
* This file must have this header.
* By default, the QC statistics will be summarized at the whole dataset. You can choose to summarize them by metadata variable.
* In case you have another type of metadata variable (e.g. stimulation condition), you could add them. For example, 'pathogen' in the 1st column (md_var) and 'condition' in the 2nd column (type).
* It is assumed that the metadata variable names are columns of the metadata file or metadata slot of the seurat object.
* The metadata variables file provided for the test dataset is the [metadata_variables.tab](/metadata_variables.tab) file: 

| md_var  | type |  
| ------------- | ------------- |  
| Pool  | donor  |  
| Assignment  | donor  |  
| predicted.celltype.l2  | cell  |  
| scpred_prediction  | cell  |  
| predicted.celltype.l1  | cell  |  


#### Downsampling file ([downsampling.tab](/downsampling.tab))
A tsv file that has in the:
* 1st column: Metadata variable name. 
* 2nd column: Number of cells to use for downsampling every level of the specified metadata variable.

*Of note*:
* Tab separated.
* It is assumed that the metadata variable name is a column of the metadata file or metadata slot of the seurat object.
* This file must have this header.
* By default, the QC statistics will be calculated using the whole dataset. You can choose to downsample the whole dataset to a specific number of cells *(n)* for each level of a specific metadata variable *(md_var)*.
* The downsampling file provided for the test dataset is the [downsampling.tab](/downsampling.tab) file: 

| md_var  | n |  
| ------------- | ------------- |  
| predicted.celltype.l1  | 100  |  

