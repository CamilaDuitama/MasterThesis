### This file calculates  SBC gene signature ############################
#### It's based on the idea of Univariate testing of Survival Data features

######## Input ###############
###For now the rows where there is no days_to_death data will be removed from ALL datasets, until we find an imputation method####
Clinical_TrainingSet <- read_csv("/Volumes/GoogleDrive/My Drive/Documents/Life Science Informatics Master/Master thesis/Data/First_data_split/TrainingSet/Clinical_TrainingSet.csv")
mRNAArray_TrainingSet <- read_csv("/Volumes/GoogleDrive/My Drive/Documents/Life Science Informatics Master/Master thesis/Data/First_data_split/TrainingSet/mRNAArray_TrainingSet.csv")
indexes<-which(is.na(Clinical_TrainingSet$days_to_death), arr.ind=TRUE)
Clinical_TrainingSet<-Clinical_TrainingSet[-indexes, ]
mRNAArray_TrainingSet<-mRNAArray_TrainingSet[-indexes,]
#### Overall Survival Time Vector (N*1)
time<-Clinical_TrainingSet[1:77,]$days_to_death
time<-log(time)
########### Event or Not Vector (N * 1) ####
censoring<-Clinical_TrainingSet[1:77,]$vital_status
######### mRNA (or miRNA expression values) (N*D) ###
mRNAArray_TrainingSet$X1<-NULL
Y.pre.train<-mRNAArray_TrainingSet[1:77,]

######## Prefiltering of the Genes ############################### ###########################
######## Using Univariate Cox's Model for Ranking of the Survival Power of the Genes ########
  

  surv.obj <- Surv(time,censoring)
  coeff.sig <- c(0)
  
  pvalue.sig <- c(0)
  
  
  calcCox = function(x){
    q1 <- unlist(summary(coxph(surv.obj ~ ., data = as.data.frame(x))))
    return(q1$logtest.pvalue)
  }
  
  
  pvalue.sig <- apply(Y.pre.train,2,calcCox)
  
  
  ###### Adjusting p-values for Multiple Test Correction
  pvalue.sig.adj <- p.adjust(pvalue.sig, method = "fdr")
  
  #### As the number of features are quite variable choose first a very loose cut-off 
  
  signature.loose <- colnames(Y.pre.train)[(pvalue.sig.adj < 0.55)] 
  
  ### Combined the P-values
  pvalue.combined <- (pvalue.sig.adj) 
  names(pvalue.combined) <- colnames(Y.pre.train)
  ## Sort it
  pvalue.combined.sort <- sort(pvalue.combined)
  ## Only select those genes which are loosely in the signature
  pvalue.combined.adj <- pvalue.combined.sort[names(pvalue.combined.sort) %in% signature.loose]
  
  
  ### Take the top 50 genes ####
  signature.sbc <- names(pvalue.combined.adj[1:50])
signature.sbc  
  
 
  