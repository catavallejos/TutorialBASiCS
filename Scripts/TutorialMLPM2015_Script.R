#################################################################
# BASiCS tutorial for the analysis of scRNA-seq datasets ########
# Catalina A. Vallejos (catalina.vallejos@mrc-bsu.cam.ac.uk) ####
#################################################################

###############################
# Setting up the R session ####
###############################

# To clean the existig R environment use
rm(list = ls())

# If not yet installed, install "BiocGenerics" from BioConductor using
#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocGenerics")

# If not yet isntalled, install BASiCS from Github using 
#library(devtools)
#install_github('catavallejos/BASiCS')

# Loading required libraries
library(BASiCS)
library(data.table) 

##########################
# Data pre-processing ####
##########################

# Directory where the data is stored (MODIFY as required)
data.path = "/Users/catalinavallejos/Documents/MRC/Projects/SCE/LaTeX/BASiCS/AnalysisMouseESC/"

# To read the matrix of expression counts (excluding metadata) use:
Counts <- read.table(file.path(data.path,
                               "GSE46980_CombinedMoleculeCounts.tab"),
                     skip = 7, sep = "\t", 
                     colClasses = c(rep("NULL", 7),
                                    rep("numeric", 96)))
dim(Counts)

# Gene names are given by 
Genes <- read.table(file.path(data.path,
                              "GSE46980_CombinedMoleculeCounts.tab"),
                    skip = 7, sep = "\t",
                    colClasses = c("character", 
                                   rep("NULL", 102)))[,1]

# Attaching gene names as 'rownames' for the matrix of expression counts
rownames(Counts) <- Genes

# We can also fix names that are not in the right format e.g.
Genes[Genes == "1-Sep"] = "Sept1" 

# Cell ids (to be used for quality control)
Cells <- read.table(file.path(data.path,
                              "GSE46980_CombinedMoleculeCounts.tab"),
                    skip = 5, nrows = 1, header = F)[-1]
Cells <- as.vector(t(Cells))
colnames(Counts) <- Cells 

# What does the data look like?
head(Counts[, 1:10], n = 10)

# Quality control: filtering of cells
QC_Info <- read.table(file.path(data.path,
                                "187_3lanes_CA.txt"),
                      header = TRUE)
GoodCells <- QC_Info$Well[QC_Info$GoodCell==1]

# Quality control: cells that seem to be MEFs
MEF <- c("D02", "E02", "A06", "H07", "D08", 
         "A09", "G10", "F12", "G12") 
GoodCells <- GoodCells[!(GoodCells %in% MEF)]

# Quality control: removing poor quality cells
CountsQC <- subset(Counts, select = Cells %in% GoodCells)
dim(CountsQC)

# Filtering of genes: total number of counts per gene
TotCountsPerGene = rowSums(CountsQC)
sum(TotCountsPerGene == 0)

# Filtering of genes: removing genes with less than 1 counts (on average) per cell
GenesInclude = I(TotCountsPerGene >= 41)
CountsQC = as.matrix(CountsQC[GenesInclude, ]) 
dim(CountsQC)

#####################
# Running BASiCS ####
#####################

# Indicator of technical spike-in genes
TechQC = grepl("SPIKE", rownames(CountsQC))

# Reading spike-in genes information
SpikeInfo <- read.table(file.path(data.path,
                                  "SilverBulletCTRLConc.txt"),
                        sep = "\t", header = TRUE,
                        colClasses = c(rep("NULL", 5),
                                       "character",
                                       rep("NULL", 3),
                                       "numeric"))
head(SpikeInfo, n = 4)

# Removing spike-in genes that did not pass the inclusion criteria
SpikeInfoQC <- SpikeInfo[SpikeInfo$Name %in% rownames(CountsQC)[TechQC],]

# Creating the BASiCS_Data object
Data = newBASiCS_Data(Counts = CountsQC, 
                      Tech = TechQC, 
                      SpikeInfo = SpikeInfoQC)


# For illustration, a very short run of BASiCS
MCMC_Output <- BASiCS_MCMC(Data, 
                           N = 100, 
                           Thin = 2, 
                           Burn = 50)

# Loading pre-computed longer chains 
# To make this faster, we use 'fread' from the library data.table
# If you don't have this library, you can replace 'fread' by 'read.table'
# However, reading the files will take much longer!
chains.path = "/Users/catalinavallejos/Documents/MRC/Teaching/MLPM2015/Tutorial/"
ChainMu=fread(file.path(chains.path,"chain_mu_BASiCS.txt"), header = TRUE)
ChainDelta=fread(file.path(chains.path,"chain_delta_BASiCS.txt"), header = TRUE)
ChainPhi=fread(file.path(chains.path,"chain_phi_BASiCS.txt"), header = TRUE)
ChainS=fread(file.path(chains.path,"chain_s_BASiCS.txt"), header = TRUE)
ChainNu=fread(file.path(chains.path,"chain_nu_BASiCS.txt"), header = TRUE)
ChainTheta=fread(file.path(chains.path,"chain_theta_BASiCS.txt"), header = TRUE)

# Creating a BASiCS_Chain object
MCMC_Output = newBASiCS_Chain(mu = as.matrix(ChainMu),
                              delta = as.matrix(ChainDelta),
                              phi = as.matrix(ChainPhi),
                              s = as.matrix(ChainS),
                              nu = as.matrix(ChainNu),
                              theta = as.matrix(ChainTheta))


# Traceplots
plot(MCMC_Output, Param = "mu", Gene = 1)
plot(MCMC_Output, Param = "delta", Gene = 1)
plot(MCMC_Output, Param = "phi", Cell = 1)
plot(MCMC_Output, Param = "s", Cell = 1)
plot(MCMC_Output, Param = "nu", Cell = 1)
plot(MCMC_Output, Param = "theta")

# Histograms
hist(displayChainBASiCS(MCMC_Output, Param = "mu")[,1])

# Cumulative means plot
plot(cumsum(displayChainBASiCS(MCMC_Output, Param = "mu")[,1])/1:1000, type = "l")

# Autocorrelation plots
acf(displayChainBASiCS(MCMC_Output, Param = "mu")[,1])

# Convergence diagnostics: Geweke
require(coda)
ChainMu_mcmc = mcmc(displayChainBASiCS(MCMC_Output, Param = "mu"))
GewekeDiag = geweke.diag(ChainMu_mcmc) # This will take a lot time because there are ~8000 genes
hist(GewekeDiag$z)

# Now, using the average across all mu[i]'s
AvMu = apply(displayChainBASiCS(MCMC_Output, Param = "mu"), 1, mean)
plot(AvMu, type = "l")
plot(cumsum(AvMu)/1:1000, type = "l")
acf(AvMu)
geweke.diag(mcmc(AvMu))

########################################
# Post-processing of BASiCS results ####
########################################

# Histogram and/or density plot
hist(displayChainBASiCS(MCMC_Output, Param = "mu")[,1], 
     main = displayGeneNames(Data)[1],
     xlab = expression(paste("Expression rate ",mu[i])))
plot(density(displayChainBASiCS(MCMC_Output, Param = "mu")[,1]),
     main = displayGeneNames(Data)[1],
     xlab = expression(paste("Expression rate ",mu[i])))


# Medians and 95% HPD intervals
MCMC_Summary <- Summary(MCMC_Output)
head(displaySummaryBASiCS(MCMC_Summary, Param = "mu"))

plot(density(displayChainBASiCS(MCMC_Output, Param = "mu")[,1]),
     main = displayGeneNames(Data)[1],
     xlab = expression(paste("Expression rate ",mu[i])))
abline(v = displaySummaryBASiCS(MCMC_Summary, Param = "mu")[1,], lty = 2)

plot(MCMC_Summary, Param = "mu", Genes = 1:100, main = "First 100 genes")
plot(MCMC_Summary, Param = "delta", Genes = 1:100, main = "First 100 genes")
plot(MCMC_Summary, Param = "phi", main = "All cells")

plot(MCMC_Summary, Param = "s", main = "All cells")
abline(h = 0.41, lty = 2)

plot(MCMC_Summary, Param = "nu", main = "All cells")
plot(MCMC_Summary, Param = "theta")

# Other summary: posterior means
MeanMu = apply(displayChainBASiCS(MCMC_Output, Param = "mu"), 2, mean)
head(MeanMu)

# Posterior medians of gene-specific parameters
plot(MCMC_Summary, Param = "mu", Param2 = "delta", log = "x", col = 8)

# Variance decomposition
VD = BASiCS_VarianceDecomp(Data, MCMC_Output)
# Most variability
head(VD)
# Least variability
tail(VD)

# Highly variable genes detection
DetectHVG <- BASiCS_DetectHVG(Data, MCMC_Output, 
                              VarThreshold = 0.79, Plot = TRUE)
head(DetectHVG$Table)

# Lowly variable genes detection
DetectLVG <- BASiCS_DetectLVG(Data, MCMC_Output, 
                              VarThreshold = 0.41, Plot = TRUE)
head(DetectLVG$Table)

#####################
# And beyond ... ####
#####################

# To include batch-effect correction
Data = newBASiCS_Data(Counts = CountsQC, 
                      Tech = TechQC, 
                      SpikeInfo = SpikeInfoQC, 
                      BatchInfo = c(rep(1, times = 20), rep(2, times = 21)))      

# Extracting denoised and normalised expression rates (to be used for downstream analyses)
DR = BASiCS_DenoisedRates(Data, MCMC_Output)
head(round(DR[, 1:10],1), n = 10)