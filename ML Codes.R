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

#Prepare the gene expression dataset
SKCM.Data <- GDCprepare(SKCM)
#Check the output of the data
head(SKCM.Data)
View(SKCM.Data)
#Check the Subsets of the metadata and add them in one data frame
SKCMMETA <- data.frame("gender" = SKCM.Data$gender,
                       "tumor_type" = SKCM.Data$tumor_descriptor,
                       "Barcode" = SKCM.Data$barcode)

#Extract the raw data from the prepared dataset
SKCMRaw <- assays(SKCM.Data)

#Select unstranded data
dim(SKCMRaw$unstranded)
View(SKCMRaw$unstranded)

#Downsize data to 10 primary and metastatic data + male and female
SelectedBarcodes <- c(subset(SKCMMETA, tumor_type == "Metastatic" & gender == "male")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Primary" & gender == "male")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Metastatic" & gender == "female")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Primary" & gender == "female")$Barcode[c(1:10)])

SelectedData <- SKCMRaw$unstranded[, c(SelectedBarcodes)]
dim(SelectedData)
View(SelectedData) #To check for consistency in expression levels

#Normalize
normalized <- TCGAanalyze_Normalization(tabDF = SelectedData, geneInfo = geneInfoHT, method = "geneLength")

#Then filter
filtered <- TCGAanalyze_Filtering(tabDF = normalized,
                                  method = "quantile",
                                  qnt.cut = 0.25)

View(filtered)
dim(filtered)

# Create labels for classification (tumor_type)
labels <- SKCMMETA[SKCMMETA$Barcode %in% colnames(filtered), "tumor_type"]

# Ensure labels are factors
labels <- factor(labels)

# View the first few labels to confirm
head(labels)

# Split the data (80% for training, 20% for testing)
set.seed(123) # For reproducibility
trainIndex <- createDataPartition(labels, p = 0.8, list = FALSE)

# Training and testing datasets
train_data <- filtered[, trainIndex]
test_data <- filtered[, -trainIndex]

# Training and testing labels
train_labels <- labels[trainIndex]
test_labels <- labels[-trainIndex]

# Install the class package if necessary
if (!requireNamespace("class", quietly = TRUE)) install.packages("class")
library(class)

# Set number of neighbors
k <- 5

# Run k-NN
predicted_labels <- knn(train = t(train_data), test = t(test_data), cl = train_labels, k = k)

# View the predicted labels
table(predicted_labels, test_labels)
# Confusion matrix to evaluate performance
confusion_matrix <- confusionMatrix(predicted_labels, test_labels)

# Print confusion matrix
print(confusion_matrix)
accuracies <- c()

for (k in 1:20) {
  predicted_labels <- knn(train = t(train_data), test = t(test_data), cl = train_labels, k = k)
  confusion_matrix <- confusionMatrix(predicted_labels, test_labels)
  accuracies[k] <- confusion_matrix$overall['Accuracy']
}

# Plot accuracy vs. K
plot(1:20, accuracies, type = "b", xlab = "K", ylab = "Accuracy")
print(confusion_matrix)
accuracies <- c()
