# Feature-Engineering-by-Network-Domain-fendR
R package that incorporates Troyanskaya/Guan's prediction modeling approach into a flexible package that can examine multiple networks and predictive data. 

##Comparing Networks
Through the incorporation of multiple networks we can evaluate how various types of networks aid in prediction of basic features.
* Physical networks: networks based on protein-protein interaction networks
* Correlative networks: networks based on co-expression of transcripts across samples
* Integrated networks: these networks are some combination of various features, e.g. Troyanskaya

##Algorithm Comparison
In addition to the gene-based approach there are additional algorithms that can infer subnetworks from a set of training data. 
`Sara had more details here`

##Datasets to evaluate
To start we should use genetic mutation data to infer response to drug response. By using datasets that also have mRNA expression data we can infer the Correlative networks described above.
* CCLE/CTRP Data: An effort by the Broad to molecularly profile and chemically probe various cancer cell lines.  That data has been downloaded into Synapse already and can be found [here](https://www.synapse.org/#!Synapse:syn5889324)
* Sanger Data: This data has not been collected yet but is a parallel response

To go beyond cell lines we can also analyze the Novartis drug response of PDX's from TCGA. 

