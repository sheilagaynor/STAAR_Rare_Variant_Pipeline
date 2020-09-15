# Set base image
FROM bioconductor/bioconductor_docker:RELEASE_3_11

# Install dependencies
RUN Rscript -e 'install.packages(c("BiocManager","Rcpp","Matrix","RcppArmadillo","readr","data.table","dplyr","doMC"))'
RUN Rscript -e 'BiocManager::install(c("SeqArray","gdsfmt","SeqVarTools","foreach","GMMAT","CompQuadForm","GENESIS","TxDb.Hsapiens.UCSC.hg38.knownGene"))'

# Install STAAR R package from source
COPY STAAR_0.9.5.tar.gz /STAAR_0.9.5.tar.gz
RUN Rscript -e 'install.packages("STAAR_0.9.5.tar.gz", repos=NULL, type="source")'

# Copy in R scripts
COPY STAAR_null_model.R STAAR_genetic_region.R /
