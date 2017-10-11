# Tutorial about the use of  BASiCS

This repository contains a tutorial for the analysis of single-cell RNA-seq datasets using BASiCS (Vallejos et al, 2015). 


**NOTE: THIS TUTORIAL WAS CREATED USING AN OLD VERSION OF BASiCS AND ITS CODE IS OUTDATED AND IS NO LONGER MAINTAINED. FOR UP-TO-DATE DOCUMENTATION, PLEASE REFER TO THE VIGNETTE AVAILABLE AT THE [BASiCS GITHUB REPOSITORY](https://github.com/catavallejos/BASiCS) OR TO OUR [WIKI PAGE](https://github.com/catavallejos/BASiCS/wiki).**



## Example dataset

We will illustrate the use of BASiCS using the mouse Embryonic Stem cells dataset described in Islam et al (2014). The data is provided in the folder 'Data' of this repository.

- **Expression counts**: GSE46980_CombinedMoleculeCounts.tab (source: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46980)

- **Quality control information**: 187_3lanes_CA.txt (source: Sten Linnarsson)

- **Input molecules of spike-in genes**: SilverBulletCTRLConc.txt (source: Sten Linnarsson).

## Downloads

Before the tutorial, please make sure to download

* The 3 data files contained in the 'Data' folder of this repository
* The 5 .txt files contained in the 'Chains2016' folder of https://www.dropbox.com/sh/0efdyln7moqtypo/AAApP0g2PamfYV4pEqDubfXva?dl=0 

## Software requirements

This tutorial was prepared using the statistical software R (https://www.R-project.org). During the tutorial we will use RStudio as an interface. Before the practical session, please make sure you install `R` version 3.3.0 or superior. Once R has been installed, the following libraries have to be installed:

* From BioConductor
	+ BiocGenerics
		- `source("http://bioconductor.org/biocLite.R"); biocLite("BiocGenerics")`
	+ scran	(this requires R 3.3.0 or superior and the latest release the BioConductor)
	 	- `source("http://bioconductor.org/biocLite.R"); biocLite("scran")`
    
* From CRAN
	+ data.table  
	+ Rcpp
	+ methods
	+ coda
	+ devtools
		- `install.packages("PACKAGE-NAME")`

* From Github
	+ BASiCS
		- `library(devtools); install_github('catavallejos/BASiCS')` 

## The tutorial

The script and results of this tutorial are available in the folder 'Materials'.


## Contact

Any doubts please contact cvallejos 'at' turing.ac.uk


## References

* Catalina A. Vallejos, John C. Marioni and Sylvia Richardson (2015)

BASiCS: Bayesian Analysis of Single-Cell Sequencing Data. *PLOS Computational Biology*

* Saiful Islam, Amit Zeisel, Simon Joost, Gioele La Manno, Pawel Zajac, Maria Kasper, Peter Lönnerberg and Sten Linnarsson (2014)

Quantitative single-cell RNA-seq with unique molecular identifiers. *Nature Methods*
