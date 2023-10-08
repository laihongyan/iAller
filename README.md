# iAller
 

# Introduction
Allergy is an autoimmune disorder described as an undesirable response of the immune system to typically innocuous substance in the environment. Studies have shown that the ability of proteins to trigger allergic reactions in susceptible individuals can be evaluated by bioinformatics tools. However, developing computational methods to accurately identify new allergenic proteins remains a vital challenge. This work aims to propose a machine learning model based on multi-feature fusion for predicting allergenic proteins efficiently. Firstly, we prepared a benchmark dataset of allergenic and non-allergenic protein sequences and pretested on it with a machine-learning platform. Then, three preferable feature extraction methods, including amino acid composition (AAC), dipeptide composition (DPC) and composition of k-spaced amino acid pairs (CKSAAP) were chosen to extract protein sequence features. Subsequently, these features were fused and optimized by Pearson correlation coefficient (PCC) and principal component analysis (PCA). Finally, the most representative features were picked out to build the optimal predictor based on random forest (RF) algorithm. Performance evaluation results via 5-fold cross-validation showed that the final model, called iAller, could precisely distinguish allergenic proteins from non-allergenic proteins. The prediction accuracy and AUC value for validation dataset achieved 91.4% and 0.97, respectively. This model will provide guide for users to identify more allergenic proteins.
data: FASTA file of all protein sequences in this work

# File Description
1. iAller.py: the main program

2. config1: parameter setting file

3. pubscripts: codes for reading FASTA file of protein sequences

4. descproteins: codes for transforming sequences into feature vectors with AAC, DPC, CKSAAP methods

5. featureselection: codes for feature selection with Pearson correlation coefficient (PCC) method

6. dimreduction: codes for dimensionality reduction with principal component analysis (PCA) method

7. machinelearning: codes for constructing, testing and saving RF model 
