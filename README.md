# Feature-Engineering-by-Network-Domain-fendR
The goal of this project is to explore the ability to incorporate prior-knowledge networks to boost the predictive power of algorithms. We aim to explore various classes of networks, various ways of mapping data to these networks and then using existing algorithms to use network-mapped data to predict a phenotype.

##Comparing Networks
Through the incorporation of multiple networks we can evaluate how various types of networks aid in prediction of basic features.
* Physical networks: networks based on protein-protein interaction networks
* Associative networks: networks based on co-expression of transcripts across samples from relevant data sources.  [Metanetwork package](https://github.com/blogsdon/metanetwork), [Metanetwork Synapse Integration](https://github.com/blogsdon/metanetworkSynapse)
* Integrated networks: these networks are some combination of various features, e.g. Troyanskaya

##Algorithm Comparison
In addition to the gene-based approach there are additional algorithms that can infer subnetworks from a set of training data.

###Step 1: mapping features to network
####Gene-based approaches
YFG used a network approach to alter gene-specific weightings based on a training set of genes known to be associated with a phenotype. These approaches are described in paper [here](http://dx.plos.org/10.1371/journal.pcbi.1000991) and [here](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002694).

####Graph-based approaches
Numerous graph-reduction algorithms have used optimization strategies to identify subnetworks of interest from a set set of genes related to a particular phenotype:
* [HotNet2](https://github.com/raphael-group/hotnet2): a network diffusion approach that can identify enriched subnetworks of genes given an underlying network and a set of weighted features such as mutations. The manuscript describing the algorithm is [here](http://www.nature.com/ng/journal/v47/n2/abs/ng.3168.html).
* [Forest/OmicsIntegrator](https://github.com/fraenkel-lab/OmicsIntegrator): this is another network reduction approach that is originally based on the Prize-Collecting Steiner Forest problem originally published [here](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002887).  The algorithm was moved from a deterministic to randomized approach and packaged in a larger python package that is described [here](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004879).

If we use these algorithms to infer subnetworks, we can update our clasification strategy to map a gene set of interest to the subnetworks that have been shown to be predictive.

###Step 2: using features to infer phenotype 
This can be done using basic classification or regression approaches, all of which can be compared. Our general consenus is that the choice of algorithm won't make a difference, but it could be wortwhile to confirm this.
* Random Forests
* SVMs
* Elastic Net
* [Goshawk](https://github.com/blogsdon/goshawk)

##Synapse Project
http://www.synapse.org/fendR

##Datasets to evaluate
To start we should use genetic mutation data to infer response to drug response. By using datasets that also have mRNA expression data we can infer the Correlative networks described above.
* CCLE/CTRP Data: An effort by the Broad to molecularly profile and chemically probe various cancer cell lines.  That data has been downloaded into Synapse already and can be found [here](https://www.synapse.org/#!Synapse:syn5889324). An effort to predict drug sensitivity from this data is described [here](http://cancerdiscovery.aacrjournals.org/content/5/11/1210.long). Other CTRP modeling approaches can be found on the [CTD2 data portal](https://ocg.cancer.gov/programs/ctd2/data-portal)
* Sanger Data: This data has not been collected yet but is a parallel approach described [here](http://www.cell.com/cell/fulltext/S0092-8674(16)30746-2).
* ROSMAP Data: predicting cognitive decline based on genotype. [ROSMAP Study](https://www.synapse.org/#!Synapse:syn3219045)

To go beyond cell lines we can also analyze the Novartis drug response of PDX's from TCGA.
