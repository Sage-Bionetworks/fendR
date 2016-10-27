# Feature-Engineering-by-Network-Domain-fendR
R package that incorporates Troyanskaya/Guan's prediction modeling approach into a flexible package that can examine multiple networks and predictive data. 

##Comparing Networks
Through the incorporation of multiple networks we can evaluate how various types of networks aid in prediction of basic features.
* Physical networks: networks based on protein-protein interaction networks
* Associative networks: networks based on co-expression of transcripts across samples
* Integrated networks: these networks are some combination of various features, e.g. Troyanskaya

##Algorithm Comparison
In addition to the gene-based approach there are additional algorithms that can infer subnetworks from a set of training data. 

###Step 1: mapping features to network
####Gene-based approaches
YFG used a network approach to alter gene-specific weightings based on a training set of genes known to be associated with a phenotype

####Graph-based approaches
Numerous graph-reduction algorithms have used optimization strategies to identify subnetworks of interest from a set set of genes related to a particular phenotype:
* [HotNet2](https://github.com/raphael-group/hotnet2)
* [Forest/OmicsIntegrator](https://github.com/fraenkel-lab/OmicsIntegrator)

If we use these algorithms to infer subnetworks, we can update our clasification strategy to map a gene set of interest to the subnetworks that have been shown to be predictive.

###Step 2: using features to infer phenotype 
This can be done using basic classification or regression approaches, all of which can be compared. 
* Random Forests
* SVMs
* Elastic Net

##Datasets to evaluate
To start we should use genetic mutation data to infer response to drug response. By using datasets that also have mRNA expression data we can infer the Correlative networks described above.
* CCLE/CTRP Data: An effort by the Broad to molecularly profile and chemically probe various cancer cell lines.  That data has been downloaded into Synapse already and can be found [here](https://www.synapse.org/#!Synapse:syn5889324)
* Sanger Data: This data has not been collected yet but is a parallel response
* ROSMap Data: predicting cognitive decline based on genotype.

To go beyond cell lines we can also analyze the Novartis drug response of PDX's from TCGA. 

