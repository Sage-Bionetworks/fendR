# fendr: From Expression (to) Novel Drug Response
This package represents our attempts to identify novel drugs and/or targets from transcriptomic and/or genomic data. 

## Installation
To get this to run we have not built a package yet but you can clone the repository as follows:
```
git clone https://github.com/Sage-Bionetworks/fendR/
```

Once you open the Rstudio project you can install the required packages:
```
install.packages('devtools')
devtools::install_github('sgosline/PCSF')
source("https://bioconductor.org/biocLite.R")
biocLite("viper")
```

## Ongoing work
We are currently experimenting with:

1- Incorporating @allaway's drug-target network to quantify putative drug targets of novel drugs
2- Mapping gene expression data to proteins in the protein and or drug-target interaction networks

See [development](./dev/) Directory for latest plans and updates
