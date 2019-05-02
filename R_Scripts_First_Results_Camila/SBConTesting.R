####################################################################################################
################ This file prepares for testing the data on external validation data set #########
###################################################################################################

######## Input #########################
######## TRained Model ##################

#### Overall Survival Time Vector (N*1)
#Clinical_TestSet <- read_csv("/Volumes/GoogleDrive/My Drive/Documents/Life Science Informatics Master/Master thesis/Data/First_data_split/TestSet/Clinical_TestSet.csv")
time.new<-Clinical_TrainingSet[78:154,]$days_to_death
time.new<-log(time.new)
time.new
########### Event or Not Vector (N * 1) ####
censoring.new<-Clinical_TrainingSet[78:154,]$vital_status
censoring.new
######### mRNA (or miRNA expression values) (N*D) ###
#mRNAArray_TestSet <- read_csv("/Volumes/GoogleDrive/My Drive/Documents/Life Science Informatics Master/Master thesis/Data/First_data_split/TestSet/mRNAArray_TestSet.csv")
#mRNAArray_TestSet$X1<-NULL
Y.pre.test<-mRNAArray_TrainingSet

######### Signature of SBC ###############
signature.sbc


######## Verhaak Signature ###############
signature.vk

######## Verhaak Labels #################
labels.vk


########################################################################
###### We prepare the Data Structures Needed for the Running of the SBC ####
#############################################################################
Y <- data.matrix(Y.pre.test[,as.character(signature.sbc)][1:77,])
Y.new <- data.matrix(Y.pre.test[,as.character(signature.sbc)][78:154,])
smod.new <- Surv(time.new, censoring.new)
c.true.new<-labels.vk[78:154]
levels(c.true.new)<-c(1,2,3,4)
c.true.new<-as.factor(c.true.new)
c.true.new
################# Scale the data ##################
Y <- scale(Y, center = TRUE, scale = TRUE)
obj <- scale(Y, center = TRUE, scale = TRUE)
Y.new <- scale(Y.new, center = attr(obj,"scaled:center"), scale = (attr(obj,"scaled:scale")))

######## Predict on New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
source('predictCLASS.R')
#What does this function do?: Outputs a matrix that has dimension number of test samples * number of MCMC samples
#The output matrix is c.matrix.new
predictCLASS(Y.new)
## Check the predicted Rand Index

#### Choose that configuration which has the highest difference in survival curves ####################
lr <- c(0)
for (j in 1:Nps){
  lr[j] <-  1 - pchisq(unlist(survdiff(smod.new ~ c.matrix.new[,j]))$chisq,df=length(table(c.matrix.new[,j]))-1)
}
c.sbc.new <- c.matrix.new[,which.min(lr)]

###############################################################################
################################# Other ways of obtaining predictions #########
## If we fit cluster-specific models ###########################################

pre.sbc <- c(0)
for ( q in 1:F){
  ind <- which(c.sbc == q)
  ind.new <- which(c.sbc.new == q)

  time.tmp <- time[ind]

  Y.tmp <- Y[ind,]
  Y.tmp.new <- Y.new[ind.new,]

  reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
  pre.sbc[ind.new] <- predict(object = reg, newx = Y.tmp.new, s = "lambda.min")
}
predCIndex.sbc.aft  <<- as.numeric(concordance(smod.new ~ exp(-pre.sbc))[1])


#### Use the adhoc-prediction algorithm ####################
source('predictADHOCverhaak.R')

surv.fit <- survfit(smod.new ~ c.sbc.new)
p2 <- ggsurv(surv.fit, plot.cens = FALSE, main = " DPMM \n Kaplan Meier Estimators \n Verhaak Cancer Test Data Set") + ggplot2::guides(linetype = FALSE) 
############ Generating some Plots ##########################
pc <- prcomp(Y.new)
pc.pred <- predict(pc,newdata = Y.new)
p3 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.sbc.new))) + ggtitle(" SBC Clustering \n Verhaak Cancer Test Data Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 
#############
logrank.new <- survdiff(smod.new ~ c.sbc.new)
#df= Degrees of freedom should be number of clusters-1
1 - pchisq(unlist(logrank.new)$chisq,df =length(logrank.new$n)-1)

source('predictTIME.R')
predictchineseAFTtime(Y.new)
