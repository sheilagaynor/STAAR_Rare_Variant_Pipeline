# Set base image
FROM uwgac/r-3.6.3-mkl:latest

# Install dependencies
RUN Rscript -e 'install.packages(c("BiocManager","Rcpp","Matrix","RcppArmadillo","readr","data.table","dplyr","doMC"))'
RUN Rscript -e 'BiocManager::install(c("SeqArray","gdsfmt","SeqVarTools","foreach","GMMAT","CompQuadForm","GENESIS","TxDb.Hsapiens.UCSC.hg38.knownGene"))'

# Install STAAR v0.9.5 R package from source
COPY STAAR_0.9.5.tar.gz /STAAR_0.9.5.tar.gz
RUN Rscript -e 'install.packages("STAAR-0.9.5.tar.gz", repos=NULL, type="source")'

# Copy in pipeline from GitHub
RUN git clone https://github.com/sheilagaynor/STAAR_Rare_Variant_Pipeline.git
