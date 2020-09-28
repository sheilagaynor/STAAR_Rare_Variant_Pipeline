# STAAR_Rare_Variant_Pipeline: Rare variant analysis methods for WGS data
Maintainer: Sheila Gaynor
Version: 0.1

## Description:
Workflow to perform aggregate rare variant tests for sequencing studies and genetic data. Implements the variant-Set Test for Association using Annotation infoRmation (STAAR) procedure, as well as SKAT, Burden, and ACAT tests for both continuous and dichotomous traits. The STAAR method incorporates qualitative functional categories and quantitative complementary functional annotations (Li and Li et al, 2020). The workflow accounts for population structure and relatedness, and scales for large whole genome sequencing studies.

## Functionality:
The workflow contains a few key steps. The workflow fits a null model for testing, incorporating the outcome, covariates, and kinship (optional). The workflow then uses the null model, genotypes, and aggregation units (optional) to run rare variant association analyses.

## Required inputs:
### Null model inputs:
- **pheno_file**: [file] File name of phenotype file
- **null_file **: [string] String for naming the null model file
- **sample_id**: [string] Column name of observation/id variable  
- **outcome**: [string] Column name of outcome variable  
- **outcome_type**: [string] Continuous or Dichotomous  

### Association test inputs:
- **geno_files**: [file] File name of GDS containing genotypes 
- **results_file **: [string] String for naming the results file

## Resulting output:
The workflow produces a copy of the null model (.Rds) and results of the aggregation test in a compressed file (.gz).


