library(RTCGAToolbox)
library(RTCGA)
library(TCGAbiolinks)
library(SummarizedExperiment)

#With TCGABiolinks with reference genome hg38
GBM_clinical <- GDCquery_clinic(project = "TCGA-GBM", type = "clinical")
patient_barcodes<-GBM_clinical$submitter_id

query<- GDCquery(project = "TCGA-GBM", 
                             data.category = "Transcriptome Profiling", 
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - FPKM")
GDCdownload(query)
expdat <- GDCprepare(query = query,
                       save = TRUE, 
                       save.filename = "GBM_Biolinks.rda")
GBM_FPKM_TCGA<-assays(expdat)

####Using RTCGA Toolbox
stddata <- getFirehoseRunningDates()
stddata
GBMData <- getFirehoseData (dataset="GBM", runDate="20160128",
                            clinic=TRUE,, mRNAArray=TRUE,miRNAArray=TRUE,Methylation=TRUE,
                            CNASNP=TRUE,CNASeq=TRUE, RNAseqNorm=TRUE)
GBM_RTCGA_RNASeq<-getData(GBMData,"RNASeq2GeneNorm")
GBM_RTCGA_mRNAArray<-GBMData@mRNAArray[[3]]@DataMatrix
