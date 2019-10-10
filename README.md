# Master Thesis
## External validation and characterisation of cancer subtypes using SBC
### [Website with consolidated results](https://camiladuitama.github.io/MasterThesis/)
The background of this master thesis is the SBC, a model that infers clinically relevant cancer subtypes, by jointly clustering molecular data along with survival datain a semi-supervised manner.The original paper and the supplementary material are on pdf version in this repository inside the folder called [papers](/papers/). A graphical representation of the model is this:

![SBC](/images/Graphical_model_SBC.png)

SBC's main features are:
- Fully bayesian approach as omics data contains a lot of noise with p >> n.
- Dirichlet Process prior to automatically infer the number of clusters.
- Molecular Data modelled as a Hierarchical Multivariate Gaussian Distribution (Mixture model).
- Survival time is modelled as Log-linear (Accelerated Failure Time) distribution with molecular covariates (Mixture model).
- L-1 regularization for the covariates of the Survival Model (Bayesian Lasso).

### Abstract
The Survival Based Bayesian Clustering (SBC) model developed by Ahmad and Fröhlich (2017), infers clinically relevant cancer subtypes, by jointly clustering molecular data along with survival data. Originally, the model was tested on a a Breast Cancer (Van De Vijver et al., 2002) and a Glioblastoma Multiforme (GBM) (Verhaak et al., 2010) data set, without any further external validation. The objective of this master thesis was to perform an external validation of the SBC, a goal that entailed two major tasks: a rigorous feature engineering and selection process that improved the known predictive ability of the model, and the characterisation of the obtained clusters and corresponding signature by delving into other types of clinical and omics data such as Copy Number Variation and miRNA.
The TCGA-GBM data set was retrieved using the Bioconductor package RTCGAToolbox and after data preprocessing, appropriate normalisation and correction for sample selection bias, a combined patient cohort of 421 samples was obtained (160 patients for the training and 261 patients for the validation set). Various feature engineering and selection techniques were explored. Every SBC model fit was done using Gibbs sampling. The best feature engineering and selection approaches were the Block HSIC-Lasso model for mRNA-based selection and a Penalized Accelerated Failure Time model on a collection of oncogenic gene sets for pathway-based selection. In both cases there was an improvement of the initial Predictive C-Index (Block HSIC-Lasso feature selection = +1.5%, PAFT feature selection = +27.6%) and Recovery C-Index (Block HSIC-Lasso feature selection = +8.7%, PAFT feature selection = +5.0%). 
The work done in this master thesis is a step forward in the validation of the SBC model on an external data set such as the TCGA-GBM patient cohort.

### Data preprocessing
- [DownloadTCGA.R](DownloadTCGA.R) is the code to download the GBM Data from the TCGA Firehose Database (Using TCGAToolbox) and the GDC Data Portal (Using TCGABiolinks). It also contains the code to get the validated miRNA targets from a list of miRNA.
- The python notebook [Data_preparation-RTCGA.ipynb](Data%20preparation-RTCGA.ipynb) is contains the preprocessing of the data(where I unify patient identifiers, verify the data is complete and in the right scale and check for the presence of the Verhaak signature and samples in our dataset), initial visualization and selection of training and testing sets*.

### Results
The folder [R_Scripts_First_Results](/R_Scripts_First_Results/) contains all the necessary scripts to train and run the model, and produce an R Notebook with the results.

All the  models trained as a result of different feature selection and pre-processing methods are presented [HERE](https://camiladuitama.github.io/MasterThesis/)

### Written Thesis
If you would like to read my master thesis you can download it from [HERE](https://github.com/CamilaDuitama/MasterThesis/blob/master/Master%20thesis%20final%20version.pdf)

