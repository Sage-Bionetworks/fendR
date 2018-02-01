# R utilities for drug response prediction

We have experimented with a number of approaches to predict drug response from related data.

# Omic data weighted with gene target data
This approach adapted from the Guan lab uses a network to weight nodes by a combination of basal expression values of the corresponding gene together with the proximity of that gene to known targets of a drug. The approach then builds a statistical model with these augmented gene features. 

# Network reduction using physical drug-target-protein networks
As an alternate approach, we are experimenting with an approach that uses gene lists of interest derived from differential expression and carries out the following analysis:
1. Identify transcription factors that give rise to the gene expression changes
2. Maps transcription factor activity to protein-protein interaction network
3. Merge protein interaction network with drug target network
4. Identify proteins-drug interactions that best interrupt the active proteins
