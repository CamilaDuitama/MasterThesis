
###############################################################################
############ This file takes in training data for Glioblastoma #################
################################################################################
#rm(list =ls())
################################################
######## Load the Data #########################

######## Input ###############
#### Overall Survival Time Vector (N*1)
library(readr)
time<-Clinical_TrainingSet[1:77,]$days_to_death
time.new<-Clinical_TrainingSet[78:154,]$days_to_death
time<-log(time)
time.new<-log(time.new)
########### Event or Not Vector (N * 1) ####
censoring<-Clinical_TrainingSet[1:77,]$vital_status
censoring.new<-Clinical_TrainingSet[78:154,]$vital_status

######### Signature of SBC ###############
signature.sbc

######## Verhaak Signature ###############
Verhaak_gene_signature <-read_csv("/Volumes/GoogleDrive/My Drive/Documents/Life Science Informatics Master/Master thesis/Data/Verhaak_Data/Verhaak_GS_according_to_Training_Set_Split1.csv", 
                                  col_names = FALSE, col_types = cols(X1 = col_skip()))
Verhaak_gene_signature<-Verhaak_gene_signature[-indexes,]
signature.vk<-Verhaak_gene_signature$X2
######## Verhaak Labels #################
labels.vk<-Clinical_TrainingSet$Subtype

########################################################################
###### We prepare the Data Structures Needed for the Running of the SBC ####
#############################################################################

Y <- data.matrix(Y.pre.train[,as.character(signature.sbc)][1:77,])
Y.new<-  data.matrix(mRNAArray_TrainingSet[,as.character(signature.sbc)][78:154,])
c.true <- as.factor(labels.vk[1:77])
c.true.new<-as.factor(labels.vk[78:154])
c.true <- unclass(c.true)
c.true.new<-unclass(c.true.new)

###### Scaling the data #########
Y <- scale(Y, center = TRUE, scale = TRUE)

######
N <- nrow(Y)
D <- ncol(Y)
smod <-  Surv(exp(time), censoring)

##### Initial number of clusters
k =4
F =k


##################### STATE OF THE ART TECHNIQUES #################################
##################### BASIC METHODS + SOME ADVANCED METHODS ########################

source('Comparisonx.R')
Comparisonx()

source('ComparisionFLX.R')
ComparisionFLX()

#source('ComparisionPReMiuM.R')
#ComparisionPReMiuM()


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 150
iter.thin  = 5
Nps = as.integer(iter/ iter.thin)


######################### Initialize the Parameters ##############################
source('initialize.R')
initialize()

###################### Start with a good configuration ###########################
#source('startSBC.R')
#startSBC()



########### Train the Model #########################################
source('burninDPMM.R')
burninDPMM()

source('gibbsDPMM.R')
gibbsDPMM()

########## Analyze the fit ##########################################
### Good feature selection from heatmap plus cindex plus randindex
source('MCMCanalyze.R')
MCMCanalyze()

###############################
### Some plots and analysis ####
#################################
#### Generating some plots and results ON Training data ###
logrank <- survdiff(smod ~ c.final)
1 - pchisq(unlist(logrank)$chisq,df =3)
surv.fit <- survfit(smod ~ c.final)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators \n Verhaak Cancer Data Set")
#+ ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('1', '2'))


############ Generating some Plots ##########################
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final))) + ggtitle(" SBC Clustering \n VerhaakCancer Data Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

