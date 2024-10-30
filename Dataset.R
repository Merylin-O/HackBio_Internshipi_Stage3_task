#Install the latest version of BiocManager, which helps install bioconductor tools
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Install the latest version
BiocManager::install(version = "3.18") # Install the latest version of Bioconductor
library("BiocManager")

#Install TCGAbiolinks, a package designed to access, analyze, and visualize data from The Cancer Genome Atlas (TCGA) project.
BiocManager::install("TCGAbiolinks")
library("TCGAbiolinks")

#Install limma & edgeR commonly used for analyzing gene expression data, especially for RNA-seq and microarray data, to find differentially expressed genes (DEGs).
BiocManager::install("limma")
BiocManager::install("edgeR")
library("limma")
library("edgeR")

#Install EDASeq for the preprocessing and normalization of RNA-seq data.
BiocManager::install("EDASeq")
library("EDASeq")

#Install gplots for plotting graphs
install.packages('gplots')
library("gplots")

#Install sesameData, which provides preprocessed and annotation data to be used with the sesame package, which is designed for analyzing DNA methylation data, especially from Illumina methylation arrays (e.g., EPIC and 450K).
BiocManager::install("sesameData")
library("sesameData")

# Install SummarizedExperiment, which provides a standardized data structure for storing and managing high-throughput genomic data
BiocManager::install("SummarizedExperiment")
library("SummarizedExperiment")

#Get an Overview of the data on Melanoma
getProjectSummary("TCGA-SKCM")

#Check the other datasets available to work on
?GDCquery

#Query the Data
SKCM <- GDCquery(project = "TCGA-SKCM",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification")

#Download the dataset
GDCdownload(SKCM)

#Prepare the dataset
SKCM.Data <- GDCprepare(SKCM)

#Check the output of the data
head(SKCM.Data)
View(SKCM.Data)

#Explore the data
SKCM.Data$gender
SKCM.Data$tumor_descriptor
SKCM.Data$ #Insert anything after the dollar sign to check a subset you want, like race, gender, etc
  
  #Check the Subsets of the metadata and add them in one data frame
  SKCMMETA <- data.frame("gender" = SKCM.Data$gender,
                         "tumor_type" = SKCM.Data$tumor_descriptor,
                         "Barcode" = SKCM.Data$barcode)

#Extract the raw data from the prepared dataset
SKCMRaw <- assays(SKCM.Data)

#Select unstranded data
dim(SKCMRaw$unstranded)
View(SKCMRaw$unstranded)

#Downsize data to 20 primary and metastatic data + male and female
SelectedBarcodes <- c(subset(SKCMMETA, tumor_type == "Metastatic" & gender == "male")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Primary" & gender == "male")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Metastatic" & gender == "female")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Primary" & gender == "female")$Barcode[c(1:10)])

SelectedData <- SKCMRaw$unstranded[, c(SelectedBarcodes)]
dim(SelectedData)
View(SelectedData) #To check for consistency in expression levels

#Data normalization and filtering
normalized <- TCGAanalyze_Normalization(tabDF = SelectedData, geneInfo = geneInfoHT, method = "geneLength")

#Then filter
filtered <- TCGAanalyze_Filtering(tabDF = normalized,
                                  method = "quantile",
                                  qnt.cut = 0.25)

View(filtered)
dim(filtered)