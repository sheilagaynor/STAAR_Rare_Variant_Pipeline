# Description: Run genomewide analysis using the STAAR package.
# Inputs:
# null_file : file containing output from null model fitting via STAAR (.Rds)
# geno_file : annotated GDS file containing the given annotation channels (.gds)
# annot_file : file containing annotations as input with columns 'chr', 'pos', 'ref', 'alt' (.Rds, .Rdata, .csv)
# results_file : string of name of results file output (string)
# agds_file : string indicating whether input geno is an agds file containing the annotations, 'None' if not provided (string)
# agds_annot_channels : comma-separated names of channels in agds to be treated as annotations (string)
# agg_file : file containing the aggregation units for set-based analysis with columns 'chr', 'pos', 'ref', 'alt', 'group_id' (.Rds, .Rdata, .csv)
# cond_file : file containing the variants to be conditioned upon with columns 'chr', 'pos', 'ref', 'alt' (.Rds, .Rdata, .csv)
# cond_geno_files : file containing the genotypes for conditional analysis; often same as geno_file (.gds)
# cand_file : file containing units/windows for candidate sets of interest with columns 'group_id' or 'chr', 'start', 'end' (.Rds, .Rdata, .csv)
# maf_thres : AF threshold below which variants will be considered in rare variant analysis, 0.05 default (numeric)
# mac_thres : AC threshold above which variants will be considered in rare variant analysis, 1 default (numeric)
# window_length : length of window for region-based analysis, 2000 default (numeric)
# step_length : length of overlap for region-based analysis, 1000 default (numeric)
# num_cores : number of cores to be used in parallelized analysis (numeric)
# num_chunk_divisions : for agg units, number of units to consider at a time within a parallel loop;
#                       for region-based, length of chunk for windows to consider at a time within a parallel loop (numeric)

## Parse arguments
args <- commandArgs(T)

## Required arguments
# File inputs
null_file <- args[1]
geno_file <- args[2]
annot_file <- args[3]
results_file <- args[4]
agds_file <- args[5]
agds_annot_channels <- args[6]
agg_file <- args[7]
cond_file <- args[8]
cond_geno_files <- args[9]
cand_file <- args[10]
# Analysis inputs
maf_thres <- as.numeric(args[11])
mac_thres <- as.numeric(args[12])
window_length <- as.numeric(args[13])
step_length <- as.numeric(args[14])
# Compute inputs
num_cores <- as.numeric(args[15])
num_chunk_divisions <- as.numeric(args[16])

#####################
# Function for input file processing
get_file <- function(file_from_args,file_type){
  # Read in file
  cat('Loading File: ',file_from_args,'\n')
  if (grepl('Rda$',kinship_file_from_args,ignore.case=TRUE) | grepl('Rdata$',kinship_file_from_args,ignore.case=TRUE)){
    file_in <- get(load(file_from_args))
  } else if (grepl('Rds$',kinship_file_from_args,ignore.case=TRUE)){
    file_in <- readRDS(file_from_args) 
  } else {
    file_in <- fread(file_from_args,stringsAsFactors=F,sep=',',header=T,data.table=F)
  }
  cat('Loaded ', file_type,' file: no. rows:',nrow(file_in),' no. cols:',ncol(file_in),'\n')
  if (sum(names(file_in) %in% c('CHR','CHROM','Chr','chrom'))>0){
    names(file_in)[names(file_in) %in%  c('CHR','CHROM','Chr','chrom')] = 'chr'	
    file_in$chr = sub('chr','',file_in$chr)
  }
  names(file_in) = tolower(names(file_in))
  # Check that required columns are available
  if (file_type=='Annotation' | file_type=='Conditional'){
    if ( !(sum(names(file_in) %in% c('chr','pos','ref','alt'))==4) ){
      stop(paste0(file_type, " file does not provide necessary input \n")) }
  }
  if (file_type=='Aggregation'){
    if ( !(sum(names(file_in) %in% c('chr','pos','ref','alt','group_id'))==5) ){
      stop(paste0(file_type, " file does not provide necessary input \n")) }
  }
  if (file_type=='Candidate'){
    if ( !(sum(names(file_in) %in% c('chr','start','end'))==3 | any(names(file_in)=='group_id')) ){
      stop(paste0(file_type, " file does not provide necessary input \n")) }
  }
  file_in
}
#https://github.com/UW-GAC/analysis_pipeline/blob/master/TopmedPipeline/R/aggregateList.R
#aggregateGRangesList function Credit to: smgogarten 
aggregateGRangesList <- function(variants) {
  stopifnot(all(c("group_id", "chr", "pos") %in% names(variants)))
  groups <- unique(variants$group_id)
  cols <- setdiff(names(variants), c("group_id", "chr", "pos"))
  GRangesList(lapply(setNames(groups, groups), function(g) {
    x <- variants[variants$group_id == g,]
    gr <- GRanges(seqnames=x$chr, ranges=IRanges(start=x$pos, width=1))
    mcols(gr) <- x[,cols]
    gr
  }))
}

# Load packages
suppressMessages(library(gdsfmt))
suppressMessages(library(SeqArray))
suppressMessages(library(STAAR))
suppressMessages(library(SeqVarTools))
suppressMessages(library(dplyr))
suppressMessages(library(doMC))
suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))

# Read in files: null model, genotypes
null_model <- readRDS(null_file)
geno_all <- seqOpen(geno_file)
cat('Read in provided null model and genotype files \n')
if (annot_file=='None' & agds_file=='None'){
  cat('No annotation file provided, proceeding without annotations \n')
} else if (annot_file!='None' & agds_file=='None') {
  annot_table <- get_file(annot_file,'Annotation')
  cat('Read in provided annotation file \n')
}
if (agg_file!='None'){
  cat('Aggregate testing based on grouping file will be used \n')
  agg_units <- get_file(agg_file,'Aggregation')
} else {
  window_length <- ifelse(window_length>0, window_length, 5000)
  step_length <- ifelse(step_length>0, step_length, window_length/2.0)
  cat(paste0('Sliding window testing will be done with window length: ',window_length,', step length: ',step_length, ' \n'))
}
if (cand_file!='None'){
  cand_in <- get_file(cand_file,'Candidate')
  cat('Read in provided candidate unit file \n')
  if (agg_file!='None' & any(names(cand_in)=='group_id')){
    agg_units <- agg_units[agg_units$group_id %in% cand_in$group_id,]
    cat("Candidate aggregate testing based on groups to be run \n")
  } else {
    cat("Candidate sliding window testing to be run \n")
  }
}

# Set up potential multi-core analysis
n_cores <- min(c(num_cores, parallel::detectCores(logical = TRUE)))

#####################
# Define sample and variants of interest
pheno_id <- as.character(null_model$id_include)
variant_id <- seqGetData(geno_all, "variant.id")
if (agds_file!='None'){
  filter <- seqGetData(geno_all, "annotation/filter")
  AVGDP <- seqGetData(geno_all, "annotation/info/AVGDP")
  SNVlist <- filter == "PASS" & AVGDP > 10 & isSNV(geno_all)
  rm(filter); rm(AVGDP)
  seqSetFilter(geno_all,sample.id=pheno_id,variant.id=variant_id[SNVlist])
} else {
  seqSetFilter(geno_all,sample.id=pheno_id,variant.id=variant_id)
}

#Below adapted from https://github.com/AnalysisCommons/genesis_wdl/blob/master/genesis_tests.R
variant_info <- variantInfo(geno_all, alleles = FALSE, expanded=FALSE)
chr <- variant_info$chr[1]
if(length(unique(chr)) > 1) stop("Multiple chromosomes detected; terminating \n")
chr = chr[1] 
#Get the aggregation units
if(agg_file!='None'){
  chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  aggVarList <- aggregateGRangesList(agg_units)
  n_chunk <- length(chunk(names(aggVarList),num_chunk_divisions))
} else {
  if (cand_file=='None'){
    #Get the genome chunks for windows
    grange_df <- data.frame(chr=chr, start=min(variant_info$pos), end=max(variant_info$pos))
    grange <- makeGRangesFromDataFrame(grange_df)
    grange$seg.length <- num_chunk_divisions
    #Get range data
    range_data <- do.call(c, lapply(seq_along(grange), function(i) {
      x <- grange[i]
      window.start <- seq(BiocGenerics::start(x), BiocGenerics::end(x), x$seg.length)
      GRanges(seqnames = seqnames(x), IRanges(window.start, width = x$seg.length))}))
    n_chunk  <- length(range_data)
  } else {
    #Get the genome chunks for windows from candidate input
    grange_df <- data.frame(chr=cand_in$chr, start=cand_in$start, end=cand_in$end)
    range_data <- makeGRangesFromDataFrame(grange_df)
    n_chunk <- length(range_data)
  }
}
seqClose(geno_all)

#####################
# Define, prepare conditional set
if (cond_file != 'None'){
  cond_in <- get_file(cond_file,'Conditional')
  cond_geno_vec <- unlist(strsplit(cond_geno_files,','))
  cond_matrix <- c()
  for (cond_geno_file in cond_geno_vec){
    cond_geno <- seqOpen(cond_geno_file)
    cond_geno_variant_info_all <- variantInfo(cond_geno, alleles = TRUE, expanded=FALSE)
    cond_geno_variant_info <- cond_geno_variant_info_all[cond_geno_variant_info_all$pos %in% cond_in$pos,]
    avail_cond_geno <- merge(cond_geno_variant_info, cond_in, by=c('chr','pos','ref','alt'))
    cond_var_list <- cond_geno_variant_info_all$variant.id[cond_geno_variant_info_all$variant.id %in% avail_cond_geno$variant.id]
    rm(cond_geno_variant_info_all); rm(cond_geno_variant_info); rm(avail_cond_geno)
    seqSetFilter(cond_geno,sample.id=pheno_id,variant.id=cond_var_list)
    cond_id.genotype.match <- match(pheno_id, seqGetData(cond_geno,"sample.id"))
    cond_genotypes <- seqGetData(cond_geno, "$dosage")
    cond_matrix <- cbind(cond_matrix, cond_genotypes[cond_id.genotype.match,])
    rm(cond_var_list); rm(cond_id.genotype.match); rm(cond_genotypes)
    seqClose(cond_geno)
  }
  cat('Prepared conditional variant lists: no. conditioning variants:',ncol(cond_matrix),'\n')
}


#####################
# Function for running tests on the large chunk
test_chunk <- function( indx ) {
  print(paste0('Iteration: ', indx))
  #Read in genotype data
  geno <- seqOpen(geno_file)
  seqSetFilter(geno,sample.id=pheno_id,verbose=F)	
  
  #First for gene based/agg unit based (candidate or full chromosome)
  #Agg unit option adapted from https://github.com/AnalysisCommons/genesis_wdl/blob/master/genesis_tests.R
  if(agg_file!='None'){
    if (agds_file!='None'){
      seqSetFilter(geno,sample.id=pheno_id,variant.id=variant_id[SNVlist],verbose=F)	
    }
    if( length(names(aggVarList)) < indx) return(data.frame())
    agg_var <- aggVarList[names(aggVarList) %in% chunk(names(aggVarList),num_chunk_divisions)[[indx]]]
    seqSetFilter(geno, variant.sel = agg_var, verbose = TRUE)
    #Subset to rare variants for efficiency
    geno_variant_id <- seqGetData(geno,'variant.id')
    if (length(geno_variant_id)>1){
      freq_vec <- seqAlleleFreq(geno)
      rare_freq_inc <- ((freq_vec <= maf_thres) | (1 - freq_vec <= maf_thres)) & (freq_vec != 0) & (freq_vec != 1)
      ct_vec <- seqAlleleCount(geno)
      rare_ct_inc <- (ct_vec > mac_thres)
      seqSetFilter(geno, variant.id=geno_variant_id[rare_freq_inc & rare_ct_inc], verbose = TRUE)
      rm(freq_vec); rm(rare_freq_inc); rm(rare_ct_inc)
      #Subset annotations for efficiency
      if(annot_file!='None'){
        geno_variant_info <- variantInfo(geno, alleles = TRUE, expanded=FALSE)
        geno_matching <- paste(geno_variant_info$chr, geno_variant_info$pos, geno_variant_info$ref, geno_variant_info$alt, sep='_')
        annot_matching <- paste(annot_table$chr, annot_table$pos, annot_table$ref, annot_table$alt, sep='_')
        annot_chunk <- annot_table[annot_matching %in% geno_matching,]
      }
      #Get iterator for looping tests
      if(length(seqGetData(geno, "variant.id"))>1){
        if(agg_file!='None'){
          iterator <- SeqVarListIterator(geno, agg_var, verbose = T)
        } 
        var_info_iter <- list(variantInfo(iterator))
        iter <- 2
        while(iterateFilter(iterator)) {
          var_info_iter[[iter]] <- variantInfo(iterator) ; iter <- iter + 1
        }
        var_info_iter <- Filter(function(x) nrow(x) > 0, var_info_iter)
        results <- c()
        for ( var_set in 1:length(var_info_iter)) {
          if (annot_file=='None' & agds_annot_channels=='None'){
            ###############################
            #Proceed without annotations
            #Subset to the genotypes of interest
            seqSetFilter(geno, sample.id=pheno_id, variant.id=var_info_iter[[var_set]]$variant.id, verbose = TRUE)
            # Match the genotype and phenotype ids
            id.genotype.match <- match(pheno_id, seqGetData(geno,"sample.id"))
            ## Get genotype	
            genotypes <- seqGetData(geno, "$dosage")
            genotypes <- genotypes[id.genotype.match,]
            pvalues <- 0
            if(cond_file=='None'){
              try(pvalues <- STAAR(genotypes,null_model))
            } else {
              try(pvalues <- STAAR_cond(genotypes,cond_matrix,null_model))
            }
            if(class(pvalues)=="list" & agg_file!='None'){
              results_temp <- c(chr, names(agg_var[var_set]), unlist(pvalues[-2]))
              results <- rbind(results,results_temp)
            } 
          } else if(annot_file=='None' & agds_annot_channels!='None') {
            ###############################
            #Proceed with aGDS annotations
            #Subset to the genotypes of interest
            seqSetFilter(geno, sample.id=pheno_id, variant.id=var_info_iter[[var_set]]$variant.id, verbose = TRUE)
            # Match the genotype and phenotype ids
            id.genotype.match <- match(pheno_id, seqGetData(geno,"sample.id"))
            ## Get genotype	
            genotypes <- seqGetData(geno, "$dosage")
            genotypes <- genotypes[id.genotype.match,]
            ## Get annotations
            annot_str_spl <- unlist(strsplit(agds_annot_channels, split=","))
            annot_tab <- c()
            for (kk in 1:length(annot_str_spl)){
              annot_tab <- cbind( annot_tab, seqGetData(geno_all, annot_str_spl[kk]))
            }
            annot_df <- data.frame(annot_tab)
            pvalues <- 0
            if(cond_file=='None'){
              try(pvalues <- STAAR(genotypes,null_model,annot_df))
            } else {
              try(pvalues <- STAAR_cond(genotypes,cond_matrix,null_model,annot_df))
            }
            if(class(pvalues)=="list" & agg_file!='None'){
              results_temp <- c(chr, names(agg_var[var_set]), unlist(pvalues[-2]))
              results <- rbind(results,results_temp)
            } 
          } else {
            ###############################
            #Proceed with annotations
            #Subset to the genotypes of interest
            seqSetFilter(geno, sample.id=pheno_id, variant.id=var_info_iter[[var_set]]$variant.id, verbose = TRUE)
            geno_iter_matching <- paste(var_info_iter[[var_set]]$chr, var_info_iter[[var_set]]$pos, var_info_iter[[var_set]]$ref, var_info_iter[[var_set]]$alt, sep='_')
            annot_iter_matching <- paste(annot_chunk$chr, annot_chunk$pos, annot_chunk$ref, annot_chunk$alt, sep='_')
            annot_iter <- annot_chunk[ annot_iter_matching %in% geno_iter_matching,]; setDT(annot_iter)
            variant_id_iter <- merge(var_info_iter[[var_set]], annot_iter, by=c('chr','pos','ref','alt'))
            if (!is.null(dim(annot_iter))){
              seqSetFilter(geno, sample.id=pheno_id, variant.id=variant_id_iter$variant.id, verbose = TRUE)
              #Match the genotype and annotation data
              geno_info_matching <- variantInfo(geno, alleles = TRUE, expanded=FALSE)
              geno_matching <- paste(geno_info_matching$chr, geno_info_matching$pos, geno_info_matching$ref, geno_info_matching$alt, sep='_')
              annot_matching <- paste(annot_iter$chr, annot_iter$pos, annot_iter$ref, annot_iter$alt, sep='_')
              id.variant.match <- match(annot_matching, geno_matching)
              annot_in <- as.data.frame(annot_iter[, c("chr","pos","ref","alt"):=NULL])
              # Match the genotype and phenotype ids
              id.genotype.match <- match(pheno_id, seqGetData(geno,"sample.id"))
              ## Get genotype	
              genotypes <- seqGetData(geno, "$dosage")
              genotypes <- genotypes[id.genotype.match,id.variant.match]
              pvalues <- 0
              if(cond_file=='None'){
                try(pvalues <- STAAR(genotypes,null_model,annot_in))
              } else {
                try(pvalues <- STAAR_cond(genotypes,cond_matrix,null_model,annot_in))
              }
              if(class(pvalues)=="list" & agg_file!='None'){
                results_temp <- c(chr, names(agg_var[var_set]), unlist(pvalues[-2]))
                results <- rbind(results,results_temp)
              } 
            }
          }
          seqResetFilter(geno)
        }
      }
    }
  }
  
  #Next for window based
  if(agg_file=='None' & cand_file=='None'){
    #Get variants in region
    range_data <- resize(range_data, width(range_data) + window_length, fix = "start")
    seg <- range_data[indx]
    #Extract variants in region
    variant_info_chunk <- variantInfo(geno, alleles = FALSE, expanded=FALSE)
    indx_vars <- (variant_info_chunk$pos>=start(seg@ranges)) & (variant_info_chunk$pos<= end(seg@ranges))
    if (agds_file!='None'){
      seqSetFilter(geno,sample.id=pheno_id,variant.id=variant_info_chunk$variant.id[SNVlist & indx_vars])
    } else {
      seqSetFilter(geno,sample.id=pheno_id,variant.id=variant_info_chunk$variant.id[indx_vars])
    }
    #Subset to rare variants for efficiency or break out
    geno_variant_id <- seqGetData(geno,'variant.id')
    if (length(geno_variant_id)>1) {
      freq_vec <- seqAlleleFreq(geno)
      rare_freq_inc <- ((freq_vec <= maf_thres) | (1 - freq_vec <= maf_thres)) & (freq_vec != 0) & (freq_vec != 1)
      ct_vec <- seqAlleleCount(geno)
      rare_ct_inc <- (ct_vec > mac_thres)
      seqSetFilter(geno, variant.id=geno_variant_id[rare_freq_inc & rare_ct_inc], verbose = TRUE)
      rm(freq_vec); rm(rare_freq_inc); rm(rare_ct_inc)
      # Match the genotype and phenotype ids
      id.genotype.match <- match(pheno_id, seqGetData(geno,"sample.id"))
      # Get genotype as matrix
      genotypes <- seqGetData(geno, "$dosage")
      genotypes <- genotypes[id.genotype.match,]
      geno_variant_rare_id <- variantInfo(geno, alleles = TRUE, expanded=FALSE)
      # Get annotations
      if(annot_file=='None' & agds_annot_channels!='None') {
        # Get annotations from agds
        annot_str_spl <- unlist(strsplit(agds_annot_channels, split=","))
        annot_tab <- c()
        for (kk in 1:length(annot_str_spl)){
          annot_tab <- cbind( annot_tab, seqGetData(geno, annot_str_spl[kk]))
        }
        annot_chunk <- data.frame(annot_tab)
      }  else if(annot_file!='None'){ ##CHECK THIS
        geno_matching <- paste(geno_variant_rare_id$chr, geno_variant_rare_id$pos, geno_variant_rare_id$ref, geno_variant_rare_id$alt, sep='_')
        annot_matching <- paste(annot_table$chr, annot_table$pos, annot_table$ref, annot_table$alt, sep='_')
        geno_annot_var <- intersect(geno_matching, annot_matching)
        annot_chunk <- annot_table[match(geno_annot_var,annot_matching),]
        genotypes <- genotypes[match(geno_annot_var,geno_matching),]
        geno_variant_rare_id <- variantInfo(geno, alleles = TRUE, expanded=FALSE)[match(geno_annot_var,geno_matching),]
      }
      #Get the genome chunks within this iteration for windows
      grange_df <- data.frame(chr=chr, start=start(seg@ranges), end=end(seg@ranges))
      grange <- makeGRangesFromDataFrame(grange_df)
      grange$seg.length <- step_length
      #Get range data
      range_data_chunk <- do.call(c, lapply(seq_along(grange), function(i) {
        x <- grange[i]
        window.start <- seq(BiocGenerics::start(x), BiocGenerics::end(x), x$seg.length)
        GRanges(seqnames = seqnames(x), IRanges(window.start, width = window_length))}))
      # Loop through the windows
      results <- c()
      for ( window_indx in 1:length(range_data_chunk)) {
        # Select the region from the geno matrix
        geno_region <- genotypes[,(geno_variant_rare_id$pos>= start(range_data_chunk[window_indx]@ranges)) & (geno_variant_rare_id$pos<= end(range_data_chunk[window_indx]@ranges))]
        # Select annotations from chunk matrix
        annot_region <- annot_chunk[(geno_variant_rare_id$pos>= start(range_data_chunk[window_indx]@ranges)) & (geno_variant_rare_id$pos<= end(range_data_chunk[window_indx]@ranges)),]
        pvalues <- 0
        if(cond_file=='None'){
          if (annot_file=='None' & agds_annot_channels=='None'){
            try(pvalues <- STAAR(geno_region,null_model))
          } else {
            try(pvalues <- STAAR(geno_region,null_model,annot_region))
          }
        } else {
          if (annot_file=='None' & agds_annot_channels=='None'){
            try(pvalues <- STAAR_cond(geno_region,cond_matrix,null_model))
          } else {
            try(pvalues <- STAAR_cond(geno_region,cond_matrix,null_model,annot_region))
          }
        }
        if(class(pvalues)=="list") {
          results_temp <- c(chr, start(range_data_chunk[window_indx]@ranges), end(range_data_chunk[window_indx]@ranges), unlist(pvalues[-2]))
          results <- rbind(results,results_temp)
        }
      }
    }
  }
  
  #Next for candidate window based
  if(agg_file=='None' & cand_file!='None'){
    #Get variants in region
    seg <- range_data[indx]
    #Extract variants in region
    variant_info_chunk <- variantInfo(geno, alleles = FALSE, expanded=FALSE)
    indx_vars <- (variant_info_chunk$pos>=start(seg@ranges)) & (variant_info_chunk$pos<= end(seg@ranges))
    if (agds_file!='None'){
      seqSetFilter(geno,sample.id=pheno_id,variant.id=variant_info_chunk$variant.id[SNVlist & indx_vars])
    } else {
      seqSetFilter(geno,sample.id=pheno_id,variant.id=variant_info_chunk$variant.id[indx_vars])
    }
    #Subset to rare variants for efficiency or break out
    geno_variant_id <- seqGetData(geno,'variant.id')
    if (length(geno_variant_id)>1) {
      freq_vec <- seqAlleleFreq(geno)
      rare_freq_inc <- ((freq_vec <= maf_thres) | (1 - freq_vec <= maf_thres)) & (freq_vec != 0) & (freq_vec != 1)
      ct_vec <- seqAlleleCount(geno)
      rare_ct_inc <- (ct_vec > mac_thres)
      seqSetFilter(geno, variant.id=geno_variant_id[rare_freq_inc & rare_ct_inc], verbose = TRUE)
      rm(freq_vec); rm(rare_freq_inc); rm(rare_ct_inc)
      # Match the genotype and phenotype ids
      id.genotype.match <- match(pheno_id, seqGetData(geno,"sample.id"))
      # Get genotype as matrix
      genotypes <- seqGetData(geno, "$dosage")
      genotypes <- genotypes[id.genotype.match,]
      geno_variant_rare_id <- variantInfo(geno, alleles = TRUE, expanded=FALSE)
      # Get annotations
      if(annot_file=='None' & agds_annot_channels!='None') {
        # Get annotations from agds
        annot_str_spl <- unlist(strsplit(agds_annot_channels, split=","))
        annot_tab <- c()
        for (kk in 1:length(annot_str_spl)){
          annot_tab <- cbind( annot_tab, seqGetData(geno, annot_str_spl[kk]))
        }
        annot_chunk <- data.frame(annot_tab)
      }  else if(annot_file!='None'){
        geno_matching <- paste(geno_variant_rare_id$chr, geno_variant_rare_id$pos, geno_variant_rare_id$ref, geno_variant_rare_id$alt, sep='_')
        annot_matching <- paste(annot_table$chr, annot_table$pos, annot_table$ref, annot_table$alt, sep='_')
        geno_annot_var <- intersect(geno_matching, annot_matching)
        annot_chunk <- annot_table[match(geno_annot_var,annot_matching),]
        genotypes <- genotypes[match(geno_annot_var,geno_matching),]
        geno_variant_rare_id <- variantInfo(geno, alleles = TRUE, expanded=FALSE)[match(geno_annot_var,geno_matching),]
      }
      #Genome chunks is the window of interest
      results <- c()
      for ( window_indx in 1:length(range_data)) {
        # Select the region from the geno matrix
        geno_region <- genotypes[,(geno_variant_rare_id$pos>= start(range_data_chunk[window_indx]@ranges)) & (geno_variant_rare_id$pos<= end(range_data_chunk[window_indx]@ranges))]
        # Select annotations from chunk matrix
        annot_region <- annot_chunk[(geno_variant_rare_id$pos>= start(range_data_chunk[window_indx]@ranges)) & (geno_variant_rare_id$pos<= end(range_data_chunk[window_indx]@ranges)),]
        pvalues <- 0
        if(cond_file=='None'){
          if (annot_file=='None' & agds_annot_channels=='None'){
            try(pvalues <- STAAR(geno_region,null_model))
          } else {
            try(pvalues <- STAAR(geno_region,null_model,annot_region))
          }
        } else {
          if (annot_file=='None' & agds_annot_channels=='None'){
            try(pvalues <- STAAR_cond(geno_region,cond_matrix,null_model))
          } else {
            try(pvalues <- STAAR_cond(geno_region,cond_matrix,null_model,annot_region))
          }
        }
        if(class(pvalues)=="list") {
          results_temp <- c(chr, start(range_data_chunk[window_indx]@ranges), end(range_data_chunk[window_indx]@ranges), unlist(pvalues[-2]))
          results <- rbind(results,results_temp)
        }
      }
    }
  }
  
  #Assemble results
  if(!exists("results")){
    results= data.frame();
  }
  seqClose(geno)
  results
}

#####################
# Function for running full analysis
#Below adapted from https://github.com/AnalysisCommons/genesis_wdl/blob/master/genesis_tests.R
run_analysis <- function( n_cores ){
  cat('Running Analysis with ', n_cores,' cores of ',num_cores,'\n')
  print(paste('Running in', n_chunk,' analysis units'))
  if (n_cores>1){
    doMC::registerDoMC(cores = n_cores)
    mc_options <- list(preschedule=FALSE, set.seed=FALSE)
    out <- foreach(i=1:n_chunk, .combine=rbind, .inorder=FALSE, .options.multicore = mc_options) %dopar% test_chunk(i)
  } else {
    out <- c()
    for (i in 1:n_chunk){
      out_temp <- test_chunk(i)
      out <- rbind(out, out_temp)
    }
  }
  cat("\nFinished analysis \n")
  out
}


#Clean up results
results <- run_analysis( n_cores )
if(!is.null(results))
{
  colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
  if(agg_file!='None'){
    colnames(results)[1:2] <- c('chr','group_id')
  } else {
    colnames(results)[1:3] <- c('chr','start_pos','end_pos')
  }
}


# Save output
write.csv(results, file=gzfile(paste0(results_file,'_chr',chr,".csv.gz")), row.names = F)