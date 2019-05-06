# Master Thesis
External Validation and characterization of Cancer subtypes using SBC

## Objectives
1. Perform an external validation of the SBC based clusters on the Verhaak study ensuring distinct survival curves for the predicted (as well as the training) clusters, using external datasets from patients diagnosed with Glioblastoma.

2. Characterise the obtained clusters and the corresponding signature looking at other data modalities for consistent patterns (e.g. Somatic Mutation Data, Copy Number, Pathway Enrichment, etc.)

## Preliminary results
### Data preprocessing
- DownloadTCGA.R is the code to download the GBM Data from the TCGA Firehose Database and also to obtain validated miRNA targets from a list of miRNA
- The python notebook Data_preparation.ipynb is contains the preprocessing of the data(where I unify patient identifiers, verify the data is complete and in the right scale and check for the presence of the Verhaak signature and samples in our dataset), initial visualization and selection of training and testing sets*
### Training and testing model (Preliminary results)
The folder R_Scripts_First_Results_Camila contains all the necessary scripts to train and run the model, and produce an R Notebook with the results.
- The file Preliminary_Results.Rmd produces the Preliminary_Results.html which is the result run only on the training set of the first split*
- The file Results.Rmd produces the Results.html which is the result from the SBC model train and tested on the whole dataset using the first split of the data

*First split was 160 samples as training(the ones which were already classified according to the 4 Verhaak labels), 338 samples as testing.

