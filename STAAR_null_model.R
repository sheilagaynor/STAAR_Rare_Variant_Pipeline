# Description: Generate a null model using the STAAR package which provides a wrapper of the GMMAT package.
# Inputs:
# pheno_file : file containing the outcome, covariates for the null model (.csv)
# null_file : file containing output from null model fitting via STAAR (.Rds)
# sample_name : column name in pheno_file for observation IDs (string)
# outcome_name : column name in pheno_file for outcome (string)
# outcome_type : type of variable outcome_name in phenofile is, 'continuous' or 'dichotomous' (string)
# covariate_names : comma-separated names of covariate variables in pheno_file to be treated as covariates (string)
# kinship_file : file containing the kinship matrix for related null model (.Rds or .Rdata)
# het_var_name : column name in pheno_file or group for heteroscedastic errors (string)

## Parse arguments
args <- commandArgs(T)

## Required arguments
# File inputs
pheno_file <- args[1]
null_file <- args[2]
# Variables
sample_name <- args[3]
outcome_name <- args[4]
outcome_type <- args[5]
# Optional inputs
covariate_names <- args[6]
kinship_file <- args[7]
het_var_name <- args[8]

#####################
# Functions for input processing
# Adapted from https://github.com/AnalysisCommons/genesis_wdl/blob/master/genesis_nullmodel.R 
get_family <- function(outcome_type_from_args) {
  if (toupper(outcome_type_from_args) == "CONTINUOUS"){
    family <- gaussian()
  } else if (toupper(outcome_type_from_args) == "DICHOTOMOUS"){
    family <- binomial()
  } else {
    stop_msg <- paste("Invalid outcome type provided:", outcome_type_from_args)
    stop(stop_msg)  }
  return(family)
}
get_kinship <- function(kinship_file_from_args){
  cat('Loading Kinship Matrix:',kinship_file_from_args,'\n')
  if (grepl('Rda',kinship_file_from_args,ignore.case=TRUE)){
    kins <- get(load(kinship_file_from_args))
  } else if (grepl('Rds',kinship_file_from_args,ignore.case=TRUE)){
    kins <- readRDS(kinship_file_from_args) 
  } else {
    kins <- as.matrix(read.csv(kinship_file_from_args,as.is=T,check.names=F,row.names=1))
  }
  cat('Loaded Kinship: no. rows:',nrow(kins),' no. cols:',ncol(kins),'\n')
  kins
}

# Load packages
suppressMessages(library(STAAR))

#####################
# Read phenotypes, covariates
pheno <- read.csv(pheno_file, header=TRUE, as.is=TRUE)
cat('Loaded Phenotypes: no. rows:',nrow(pheno),' no. cols:',ncol(pheno),'\n')
# Subset to complete cases for phenotype file, create null model formula
if ( covariate_names=='NA' & het_var_name=='NA' ){
  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name)]), c(sample_name,outcome_name)]
  null_model_char <- paste0(outcome_name, "~", 1)
} else if ( covariate_names=='NA' & het_var_name!='NA' ){
  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name,het_var_name)]), c(sample_name,outcome_name,het_var_name)]
  null_model_char <- paste0(outcome_name, "~", 1)
} else if ( covariate_names!='NA' & het_var_name=='NA' ){
  covar_split <- unlist(strsplit(covariate_names, split=","))
  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name,covar_split)]), c(sample_name,outcome_name,covar_split)]
  null_model_char <- paste0(outcome_name, "~", paste(covar_split, collapse="+"))
} else {
  covar_split <- unlist(strsplit(covariate_names, split=","))
  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name,het_var_name,covar_split)]), 
                           c(sample_name,outcome_name,het_var_name,covar_split)]
  null_model_char <- paste0(outcome_name, "~", paste(covar_split, collapse="+"))
}
cat('Complete Phenotypes: no. rows:',nrow(phenotype_input),' no. cols:',ncol(phenotype_input),'\n')

#####################
# Get kinship
if ( kinship_file!='NA' ){
  kinship_input <- get_kinship(kinship_file)
  shared_obs <- intersect(row.names(kinship_input), phenotype_input[,sample_name])
  if (length(shared_obs)==0){
    stop('No shared observations between phenotypes and kinship')
  } else {
    kinship_analysis <- kinship_input[rownames(kinship_input) %in% phenotype_input[,sample_name],rownames(kinship_input) %in% phenotype_input[,sample_name]]
    rm(kinship_input)
    phenotype_analysis <- phenotype_input[match(rownames(kinship_analysis),phenotype_input[,sample_name]),]
    cat('Matched Phenotypes with Kinship: no. rows:',nrow(phenotype_input),' no. cols:',ncol(phenotype_input),'\n')
  }
} else {
  phenotype_analysis <- phenotype_input
}

#####################
# Fit, save null model
if ( kinship_file=='NA' & het_var_name=='NA' ){
  cat('Fitting null model for unrelated samples, homogeneous variance')
  null_model <- STAAR::fit_null(as.formula(null_model_char), id = sample_name, 
                data = phenotype_analysis, family = get_family(outcome_type))
} else if ( kinship_file=='NA' & het_var_name!='NA' ){
  cat('Fitting null model for unrelated samples, heterogeneous variance')
  null_model <- STAAR::fit_null(as.formula(null_model_char), id = sample_name, 
                groups = het_var_name, data = phenotype_analysis, family = get_family(outcome_type))
  
} else if ( kinship_file!='NA' & het_var_name=='NA' ){
  cat('Fitting null model for related samples, homogeneous variance')
  null_model <- STAAR::fit_null_glmmkin(as.formula(null_model_char), id = sample_name, 
                data = phenotype_analysis, kins = kinship_analysis, family = get_family(outcome_type))
  
} else {
  cat('Fitting null model for related samples, heterogeneous variance')
  null_model <- STAAR::fit_null_glmmkin(as.formula(null_model_char), id = sample_name, 
                groups = het_var_name, data = phenotype_analysis, kins = kinship_analysis, family = get_family(outcome_type))
}
saveRDS(null_model, file=paste0(null_file,".Rds"))
