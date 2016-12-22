#Ongoing development

This directory contains the files and information used to develop the fendR package. This directory has been added to the [.Rbuildignore][../.Rbuildignore] file.

##Development code
The files in here are used to develop and test various aspects. We can create sub-directories as need be. So far this is limited to:
* [formatCcleData.R][formatCcleData.R]: a script that pulls data from cBioPortal and pushes matrices to Synapse and to inst directory
* [exampleForestRun.R][exampleForestRun.R]: a script that runs through the forest implementation of fendR

##Development Strategy
The goal of the architecture is to create one S4 class [fendR][../R/fendR.R] that is populated with the data objects as needed and can perform the three (so far) basic operations on those algorithms. 

Actual implementation of these methods will be done by sub-classes, of which [forestFendR][../R/forestFendR.R] is the first example. This class, for example, has an additional 'slot' for the executable path to the Forest code and will also have Forest-specific implementations of these methods. 

##Open Items
There are still a few open items including:
1. How to include additional datasets, and where
2. How to include various modeling (e.g. elasticNet vs SVM) approaches - would these be sub-classes of the network approach classes?

I think these issues can be handled after we have a working prototype.



