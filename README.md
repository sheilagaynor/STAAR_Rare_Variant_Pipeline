# STAAR_Rare_Variant_Pipeline: Rare variant analysis methods for WGS data
Maintainer: Sheila Gaynor
Version: 0.1

##Description:
Workflow to perform aggregate rare variant tests for sequencing studies and genetic data. Implements the variant-Set Test for Association using Annotation infoRmation (STAAR) procedure, as well as SKAT, Burden, and ACAT tests for both continuous and dichotomous traits. The STAAR method incorporates qualitative functional categories and quantitative complementary functional annotations (Li and Li et al, 2020). The workflow accounts for population structure and relatedness, and scales for large whole genome sequencing studies.

##Functionality:
The workflow contains a few key steps. The workflow fits a null model for testing, incorporating the outcome, covariates, and kinship (optional). The workflow then uses the null model, genotypes, and aggregation units (optional) to run rare variant association analyses.

##Resulting output:
The workflow produces a copy of the null model (.Rds) and results of the aggregation test in a compressed file (.gz).
