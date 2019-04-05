FROM rocker/tidyverse

MAINTAINER Sara Gosline

RUN apt-get install -y net-tools
RUN apt-get update -qq && apt-get -y install libffi-dev

RUN Rscript -e "install.packages('argparse')"
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages('synapser', repos=c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" -e "biocLite('viper')" -e "biocLite('topGO')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" -e "biocLite('org.Hs.eg.db')"
RUN Rscript -e "devtools::install_github('sgosline/PCSF')" -e "devtools::install_github('Sage-Bionetworks/fendR')"