# wg3-sc_pseudobulk_DEA

We provide **two main scripts** to peform **differential expression analysis (DEA)** with human phenotypes (e.g., sex or age) using single-cell RNA-seq data (scRNA-seq) (i.e., 10x Genomics) at different levels:

* **[single-cell level (sc-DEA)](scDEA_MAST_glmer.R)**: using the [MAST](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5) glmer implementation.
* **[pseudobulk level (pseudobulk-DEA)](pseudobulkDEA_limmadream.R)**: using the [limma dream](https://academic.oup.com/bioinformatics/article/37/2/192/5878955) glmer implementation.

**Of note**: 
* This analysis is meant to be run on scRNA-seq data composed by **only one sample per donor**. In case you have more than one sample per donor (e.g., stimulated vs. non-stimulated samples from the same donor, etc...) you can run our **[subsetting script](subset_by_metadata.R)** to select individuals with specifics characteristics (e.g., only non-stimulated samples, etc...). 

* To run these scripts you should have **successfully run** the following sc-eQTLGen consortium pipelines: **WG1**, **WG2** and **WG3 (I)** 

## Contact
If you have any questions or issues, feel free to open an issue or directly email Aida Ripoll-Cladellas (aida.ripoll@bsc.es)

## Required Software
**R** >=4.1.2 version: You need to install the packages loaded in the:
* Two main DEA scripts: [sc-DEA](scDEA_MAST_glmer.R) and [pseudobulk-DEA](pseudobulkDEA_limmadream.R)
* Add-on subsetting script: [subsetting script](subset_by_metadata.R)
* Additional scripts in the [scripts](scripts) directory

## Required Input
This section explains the input data and it’s structure to run the three scripts: [sc-DEA](scDEA_MAST_glmer.R), [pseudobulk-DEA](pseudobulkDEA_limmadream.R), and [subsetting script](subset_by_metadata.R).

*Of note*: To follow better the explanations in the **Required Input** section, you can clone this repository and change your current working directory.   
```
git clone https://github.com/aidarripoll/wg1-qc_filtering.git  
cd wg1-qc_filtering
```
