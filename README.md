# **Gender-Based Comparative Gene Expression and Functional Enrichment in Melanoma Samples**
 
## **Description**
This project explores the role of sex in Melanoma. It is split into multiple parts
1. Differential Expression Analysis
2. Functional enrichment analysis
3. k-NN Analysis

Dataset:TCGA-SKCM

*Biomarker Discovery*

You will need the following packages in RStudio:
![68747470733a2f2f6c68372d72742e676f6f676c6575736572636f6e74656e742e636f6d2f646f63737a2f41445f346e5866376371615537676b4363524e70654232643943662d583859664b6d646164324a5f78346733566d50325557666c75637830554b4434506774497](https://github.com/user-attachments/assets/ce2b1181-4157-4b3f-bec5-3c1f3dff3600)

These packages can be installed using install.packages()

### **Preprocessing**
1. Dataset is uploaded directly onto RStudio using TCGABiolink, then preproccessed by:
  1. Subsetting into the designated groups 
     ![68747470733a2f2f6c68372d72742e676f6f676c6575736572636f6e74656e742e636f6d2f646f63737a2f41445f346e58633544567153466e654d674f72584f587138415f6e4a52466c5152487765584b37546c484b4e514c3972727272574535796c4c5f3138304c67566](https://github.com/user-attachments/assets/f8abf333-3a15-4d50-b310-f8c9c9ae9e9b)
  2. Normalization
  3. Filteration with a quantile cut off 0.25 

### **Differential Expressionn Analysis**
1. Metastatic and primary data is grouped together
2. Treatment level DEA is performed
3. Ouput is used as data to generate heatmaps comparing the clustering of males vs females in both tumor types

### **Functional Enrichment Analysis**
1. The up and downregulated genes for all 4 subsets are extracted individually
2. The data undergoes functional enrichment analysis (particularly gene ontology), selecting the top 5 pathways.
3. Results are visualized in a barplot form and the PDF is saved on the computed

*Biomarker Discovery*
1. Perform k-NN analysis  to classify melanoma samples
     - The code was implemented by using the class package in R.
     - 80/20 train-test split is employed.
     - The k-NN algorithm execution with varying values of (k) (from 1 to the number of training samples) to optimize accuracy.
2. A confusion matrix is generated 
3. An analysis vs K plot are generated

