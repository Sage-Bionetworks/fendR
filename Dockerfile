FROM rocker/tidyverse

RUN Rscript -e "install.packages('argparse')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))"
RUN Rscript -e "devtools::install_github('sgosline/PCSF')"
RUN Rscript -e "install.packages("viper")

RUN Rscript -e "BiocManager::install('preprocessCore', version = '3.8')"

RUN RCMD build fendR
RUN RCMD INSTALL fendR_0.1.tar.gz
