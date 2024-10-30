# **Gender-Based Comparative Gene Expression and Functional Enrichment in Melanoma Samples** 

Authors: Nour Nahtay (@NourNahtay), Igwebuike Oluchukwu Vivian (@igwebuikevee0000), Rutuja Pangare (@rutuja0502), Merylin Ogunlola (@MerylinO), Onare Opeyemi Mary (@Onare),  Titilola Shittu (@lola), Princess Beatrice Sunday-Jimmy (@LaidyCharm)

Melanoma is skin cancer that originates in the melanocytes or the melanin-producing skin cells. It occurs due to UV-induced DNA damage, leading to tumorigenesis and uncontrollable proliferation. (_Davis, Lauren E et al, 2019)_ The global melanoma burden is projected to escalate to 510,000 cases and 96,000 deaths annually by 2040, emphasizing the urgent need for continued research to combat this lethal disease.

***Sex Differences and Outcomes in Melanoma Patients***

Studies show that males have less survival rates than females despite having more female cases, suggesting that the sex-based differences between both sexes play a major role in dictating survival. (_Arnold et al, 2022)_

***Analysis Pipeline***

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXc0YzHDArNngwGWgrK_Yh-fWEVTeWR966UZIhEJykfk8tXauX5FItCe7nEUxTJ5vjrU6osSc2qkhEaqtgmERiNYd25-ebvBpJX8GkMi0CMOkgKSMAplKlQhiwMeUUX8xLzjFLgT7i1AcXWa7vh2PgsLCOsI?key=GNpt1ZRV2jJ52mDok30oRw)****

Figure 1: An Overview of the Pipeline Used for the Biomarker Discovery Portion of the Analysis 

***Preprocessing Data***

The transcriptome profiles from the TCGA-SKCM project were loaded, subsetted, normalized, then filtered with a quantile cut of 0.25 to eliminate background noise.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXc5DVqSFneMgOrXOXq8A_nJRFlQRHweXK7TlHKNQL9rrrrWE5ylL_180LgVcI3PQb9tRbICWypPe9ouxBBY0aCqsn4Z4fBNQp_Ri2ZH_fuNP236PEz5_VZ3eIcx3S-BMf9iXashIeJY6xARFka7pO7wKlnx?key=GNpt1ZRV2jJ52mDok30oRw)

Figure 2: The 4 Subsets Used for the Analysis (n=40)

***1. Biomarker Discovery Pipeline***

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXf7cqaU7gkCcRNpeB2d9Cf-X8YfKmdad2J_x4g3VmP2UWflucx0UKD4PgtIy6J3HViQTaWS2hK6184d2OudvF231aoYXuttr0YorN91y4RSH0JZjVYuoEHvI36vqumMRTcgLQoIQE9Zxjbhe0XnVAWvh80?key=GNpt1ZRV2jJ52mDok30oRw)****

Figure 3: An Overview of the Packages Used for Biomarker Discovery Pipeline

1. ***_Differential Expression Analysis (DEA)_***

The data is split into 2 tumor-based groups to compare the DEGs in both sexes. DEA was performed with a FDR cutoff of <0.01 and a logFC of 2 for accuracy. The output was then used for comparison using heatmaps.

2. ***_Functional Enrichment Analysis (Gene Ontology)_***

The up and downregulated genes for the 4 groups were isolated with a logFC cut off of 2 (up) and -2 (down), then prepped for gene ontology. The top 5 pathways were then visualized for comparison.

***Result Interpretation for Biomarker Discovery Findings***

DEA revealed 1,016 DEGs in primary samples and 1,014 DEGs in the metastatic samples with distinct clustering between both sexes in both tumor types, highlighting gender-specific differences. Despite having overall differing pathways, the top 5 functional enrichment results were the same for both sexes in both tumor types.

***2. Machine Learning*** 

Classification took place using k-NN analysis, using the class package and the 80/20 train-test split. The k-NN algorithm was executed with k values from 1 to sample size and a confusion matrix was used for a thorough evaluation of a classification model.

***Results and Interpretation***

The model classified only 25% of the sample, indicating low prediction rates with the P-Value of 0.9648 showing that is not significantly better than random guessing, making it unreliable. The model predicted 75% of samples to be metastatic, indicating bias and inability to predict non-metastatic samples. _(Gupta et. al, 2022)_

***Conclusion and Future Direction***

The bias towards metastatic samples shows the need for feature tuning. Accuracy could be increased by investigating different techniques, such as random forests.

Analyses provided insights into gene expression patterns, confirming sex-based differences. Pathways such as atherosclerosis being interlinked with melanoma suggest further exploration of the interplay between melanoma and other disorders in the etiology of cancer. 

Overall, increasing sample size and including various demographics can provide valuable insights and improved accuracy of analyses

***Citations:***

1. _Davis, Lauren E et al. “Current state of melanoma diagnosis and treatment.” Cancer biology & therapy vol. 20,11 (2019): 1366-1379. doi:10.1080/15384047.2019.1640032_

2. Arnold, Melina et al. “Global Burden of Cutaneous Melanoma in 2020 and Projections to 2040.” _JAMA dermatology_ vol. 158,5 (2022): 495-503. doi:10.1001/jamadermatol.2022.0160

3. Gupta, Sunil, et al. "Comparing the performance of machine learning algorithms using estimated accuracy." _Measurement: Sensors_ 24 (2022): 100432.

**Word Count (Excluding in-text citations):** 403


