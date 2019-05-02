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
# clinical,miRNASeqGene(0*0), RNASeq2GeneNorm(0*0),Methylation,
#mRNAArray,miRNAArray,GISTIC,RNAseqNorm,RNAseq2Norm
GBMData <- getFirehoseData (dataset="GBM", runDate="20160128",
                            clinic=TRUE,mRNASeq= TRUE, mRNAArray=TRUE,miRNAArray=TRUE,Methylation=TRUE,
                            GISTIC =TRUE, RNAseqNorm= TRUE, RNAseq2Norm=TRUE,miRNASeqGene=TRUE)
clinicData <- getData(GBMData,"clinical")
mRNAArray <- getData(GBMData,"mRNAArray",platform = )
write.table(GBMData@clinical, file = "Clinical.csv")
write.table(GBMData@miRNASeqGene,file="miRNASeqGene.csv")
write.table(GBMData@RNASeqGene,file="RNASeqGene.csv")
write.table(GBMData@RNASeq2GeneNorm,file="RNASeq2GeneNorm.csv")
write.csv(GBMData@mRNAArray[[3]]@DataMatrix,file="gdac.broadinstitute.org_GBM.Merge_transcriptome__ht_hg_u133a__broad_mit_edu__Level_3__gene_rma__data.Level_3.2016012800.0.0.csv")
write.csv(GBMData@miRNAArray[[1]]@DataMatrix,file="gdac.broadinstitute.org_GBM.Merge_mirna__h_mirna_8x15k__unc_edu__Level_3__unc_DWD_Batch_adjusted__data.Level_3.2016012800.0.0.csv")
write.csv(GBMData@Methylation[[1]]@DataMatrix,file="gdac.broadinstitute.org_GBM.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.csv")
write.csv(GBMData@GISTIC@AllByGene,file="GISTIC_AllByGene.csv")

diffGeneExprs <- getDiffExpressedGenes(dataObject=GBMData,DrawPlots=TRUE,
                                       adj.method="BH", adj.pval=0.05, raw.pval=0.05, logFC=2, hmTopUpN=100,
                                       hmTopDownN=100)

#Code to use multiMiR to get the validated targets of a list of miRNA
library(multiMiR)
library(readr)
miRNA_GBM_Firebrowse <- read_csv("miRNA_GBM_Firebrowse.csv", 
                                   +     col_names = FALSE)
result<-list()
for (i in names(miRNA_GBM_Firebrowse))
{
  example1 <- get_multimir(mirna = 'hsa-miR-18a-3p', summary = TRUE)
  result <- c(result, example1@data$target_symbol) 
}
result<-unlist(result, recursive=FALSE)
write.table(result,file="miRNA_targets_GBM_Firebrowse.csv", row.names = FALSE, col.names = FALSE, sep = ",")



