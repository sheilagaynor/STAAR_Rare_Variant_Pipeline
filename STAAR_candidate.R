# Description: Generate a null model using the STAAR package using a wrapper of the GMMAT package.
# Inputs:
# null_file : file containing output from null model fitting via STAAR (.Rds)
# geno_file : annotated GDS file containing the given annotation channels (.gds)
# agg_file : file variant masks, with a column for gene name and a column for unique variant identifier (.Rds)
# num_cores : number of cores for multi-threading when looping tests (numeric)

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
cand_file <- args[8]
# Analysis inputs
maf_thres <- as.numeric(args[9])
mac_thres <- as.numeric(args[10])
window_length <- as.numeric(args[11])
step_length <- as.numeric(args[12])

#####################
# Functions for input processing
get_agg <- function(agg_file_from_args){
  cat('Loading Agg File:',agg_file_from_args,'\n')
  if (grepl('Rda',agg_file_from_args,ignore.case=TRUE)){
    agg <- get(load(agg_file_from_args))
  } else if (grepl('Rds',agg_file_from_args,ignore.case=TRUE)){
    agg <- readRDS(agg_file_from_args) 
  } else {
    agg <- fread(agg_file_from_args,stringsAsFactors=F, sep=',', header=T,data.table=F)
  }
  cat('Loaded Agg Units: no. rows/units:',nrow(agg),' no. cols:',ncol(agg),'\n')
  names(agg)[names(agg) %in%  c('CHR','CHROM','Chr','chrom')] = 'chr'	
  agg$chr = sub('chr','',agg$chr)
  names(agg) = tolower(names(agg))
  agg
}
get_cand <- function(cand_file_from_args){
  cat('Loading Cand File:',cand_file_from_args,'\n')
  if (grepl('Rda',cand_file_from_args,ignore.case=TRUE)){
    cand <- get(load(cand_file_from_args))
  } else if (grepl('Rds',cand_file_from_args,ignore.case=TRUE)){
    cand <- readRDS(cand_file_from_args) 
  } else {
    cand <- fread(cand_file_from_args,stringsAsFactors=F, sep=',', header=T,data.table=F)
  }
  cat('Loaded Candidate Units: no. candidate sets:',nrow(cand),'\n')
  if (sum(names(cand) %in% c('CHR','CHROM','Chr','chrom'))>0){
    names(cand)[names(cand) %in%  c('CHR','CHROM','Chr','chrom')] = 'chr'	
    cand$chr = sub('chr','',cand$chr)
  }
  names(cand) = tolower(names(cand))
  if ( !(sum(names(cand) %in% c('chr','start','end'))==3 | any(names(cand)=='group_id')) ){
    stop("Candidate file does not provide necessary input \n")
  }
  cand
}
#https://github.com/UW-GAC/analysis_pipeline/blob/master/TopmedPipeline/R/aggregateList.R
#Credit to: smgogarten 
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
get_annot <- function(annot_file_from_args){
  cat('Loading Annot File:',annot_file_from_args,'\n')
  if (grepl('Rda',annot_file_from_args,ignore.case=TRUE)){
    annot <- get(load(annot_file_from_args))
  } else if (grepl('Rds',annot_file_from_args,ignore.case=TRUE)){
    annot <- readRDS(annot_file_from_args) 
  } else {
    annot <- fread(annot_file_from_args,stringsAsFactors=F, sep=',', header=T,data.table=F)
  }
  cat('Loaded Annot Units: no. rows/units:',nrow(annot),' no. cols:',ncol(annot),'\n')
  names(annot)[names(annot) %in%  c('CHR','CHROM','Chr','chrom')] = 'chr'	
  annot$chr = sub('chr','',annot$chr)
  names(annot) = tolower(names(annot))
  annot
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

# Read in files: null model, genotypes, candidate set
null_model <- readRDS(null_file)
geno_all <- seqOpen(geno_file)
cand_in <- get_cand(cand_file)
if (annot_file=='None' & agds_file=='None'){
  cat('No annotation file provided, proceeding without annotations \n')
} else if (annot_file!='None' & agds_file=='None') {
  annot_table <- get_annot(annot_file)
  cat('Read in provided annotation file \n')
}
if (agg_file!='None' & any(names(cand_in)=='group_id')){
  cat('Aggregate testing based on grouping file will be used \n')
  agg_units_all <- get_agg(agg_file)
  agg_units <- agg_units_all[agg_units_all$group_id %in% cand_in$group_id,]
} 

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
}
seqClose(geno_all)


#####################
# Run tests
# Below adapted from https://github.com/AnalysisCommons/genesis_wdl/blob/master/genesis_tests.R

  #Read in genotype data
  geno <- seqOpen(geno_file)
  seqSetFilter(geno,sample.id=pheno_id,verbose=F)	
  #Subset genotype data to chunk
  if (agg_file!='None' & any(names(cand_in)=='group_id')){
    agg_var <- aggVarList
    seqSetFilter(geno, variant.sel = agg_var, verbose = TRUE)
  } else {
    grange_df <- data.frame(chr=chr, start=cand_in$start, end=cand_in$end)
    seg <- makeGRangesFromDataFrame(grange_df)
    seqSetFilter(geno, variant.sel = seg, verbose = TRUE)
  }
  #Subset to rare variants for efficiency
  geno_variant_id <- variant_info$variant_id
  freq_vec <- seqAlleleFreq(geno)
  rare_inc <- ((freq_vec <= maf_thres) | (1 - freq_vec <= maf_thres)) & (freq_vec != 0) & (freq_vec != 1)
  seqSetFilter(geno, variant.id=geno_variant_id[rare_inc], verbose = TRUE)
  rm(freq_vec); rm(rare_inc)
  #Subset annotations for efficiency
  if(annot_file!='None'){
    geno_variant_info <- variantInfo(geno, alleles = TRUE, expanded=FALSE)
    geno_matching <- paste(geno_variant_info$chr, geno_variant_info$pos, geno_variant_info$ref, geno_variant_info$alt, sep='_')
    annot_matching <- paste(annot_table$chr, annot_table$pos, annot_table$ref, annot_table$alt, sep='_')
    annot_chunk <- annot_table[annot_matching %in% geno_matching,]
  }
  
  #Get iterator for looping tests
  if(length(seqGetData(geno, "variant.id"))>0){
    if(agg_file!='None'){
      iterator <- SeqVarListIterator(geno, agg_var, verbose = T)
    } else {
      iterator <- SeqVarRangeIterator(seqData, variantRanges=seg, verbose = T)
    }
    var_info_iter <- list(variantInfo(iterator))
    iter <- 2
    while(iterateFilter(iterator)) {
      var_info_iter[[iter]] <- variantInfo(iterator) ; iter <- iter + 1
    }
    var_info_iter <- Filter(function(x) nrow(x) > 0, var_info_iter)
    results <- c()
    for ( var_set in 1:length(var_info_iter)) {
      if (annot_file=='None'){
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
          try(pvalues <- STAAR(genotypes,null_model))
          if(class(pvalues)=="list" & agg_file!='None'){
            results_temp <- c(chr, names(agg_var[var_set]), unlist(pvalues[-2]))
            results <- rbind(results,results_temp)
          } else if(class(pvalues)=="list") {
            results_temp <- c(chr, start(iterator@variantRanges@ranges[var_set]), end(iterator@variantRanges@ranges[var_set]), unlist(pvalues[-2]))
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
        try(pvalues <- STAAR(genotypes,null_model,annot_df))
        if(class(pvalues)=="list" & agg_file!='None'){
          results_temp <- c(chr, names(agg_var[var_set]), unlist(pvalues[-2]))
          results <- rbind(results,results_temp)
        } else if(class(pvalues)=="list") {
          results_temp <- c(chr, start(iterator@variantRanges@ranges[var_set]), end(iterator@variantRanges@ranges[var_set]), unlist(pvalues[-2]))
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
            try(pvalues <- STAAR(genotypes,null_model,annot_in))
            if(class(pvalues)=="list" & agg_file!='None'){
              results_temp <- c(chr, names(agg_var[var_set]), unlist(pvalues[-2]))
              results <- rbind(results,results_temp)
            } else if(class(pvalues)=="list") {
              results_temp <- c(chr, start(iterator@variantRanges@ranges[var_set]), end(iterator@variantRanges@ranges[var_set]), unlist(pvalues[-2]))
              results <- rbind(results,results_temp)
            }
          }
        }
        seqResetFilter(geno)
      }
    }
  #Assemble results
  if(!exists("results")){
    results= data.frame();
  }


#Clean up results
if(!is.null(results))
{
  colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
  if(agg_file!='None'){
    colnames(results)[1:2] <- c('chr','group_id')
  } else {
    colnames(results)[1:3] <- c('chr','star_pos','end_pos')
  }
}


# Save output
write.csv(results, file=gzfile(paste0(results_file,'_chr',chr,".csv.gz")), row.names = F)