####################################################################################################################################################
####################################################################################################################################################
###################### VISUALIZATION ###############################################################################################################
####################################################################################################################################################

##########################################################
########## Creating a Ranking for the Points ####################
###### Creating a Ratio for the Points ########################
rank <- matrix(0, nrow = N, ncol =Nps)

for (j in 1:Nps){
  
  Y.scaled <- matrix(0, nrow = N, ncol =D)
  for ( v in 1:4){
    clust <- which(c.list[[j]] == v)
    Y.scaled[clust,1:D] <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
  }
  
  
  for ( i in 1:N){
    rank[i,j] <- dMVN(as.vector(t(Y[i,1:D])), mean = mu.list[[j]][c.final[i],1:D], Q= S.list[[j]][c.final[i],1:D,1:D], log = TRUE) - dMVN(as.vector(t(Y[i,1:D])), mean = mu.list[[j]][1,1:D], Q= S.list[[j]][1,1:D,1:D], log = TRUE) +  dnorm(x = That[i], mean = beta0.list[[j]][c.final[i]] + betahat.list[[j]][c.final[i],1:D] %*% as.vector(t(Y.scaled[i,1:D])), sd = sqrt(sigma2.list[[j]][c.final[i]]), log =TRUE) -  dnorm(x = That[i], mean = beta0.list[[j]][1] + betahat.list[[j]][1,1:D] %*% as.vector(t(Y.scaled[i,1:D])), sd = sqrt(sigma2.list[[j]][1]), log =TRUE) 
  }
}

avg.rank <- apply(rank,1,mean)
order.zo <- range01(avg.rank)

order.train <- sort(order.zo,index.return = TRUE, decreasing = TRUE)
Y.order <- Y[order.train$ix,]
c.final.order <- c.final[order.train$ix]


######## Reordering Again ################
order.2 <- c(which(c.final.order==1),which(c.final.order==3),which(c.final.order==4),which(c.final.order==2))
Y.order.2 <- Y.order[order.2,]
c.final.order.2 <- c.final.order[order.2]




surv.ob <- Surv(exp(time),censoring)
Classes <- c.final
Classes <- as.factor(Classes)
surv.fit <- survfit(surv.ob ~ Classes)
p5 <- ggsurv(surv.fit, main = " Training Data Set\n Kaplan Meier Estimators \n pvalue 5e -05", surv.col = c("green","red","blue","orange"), cens.col ="Blue")  


############ Generating some Plots ##########################
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], color = as.factor(c.final) )) + ggtitle(" DPMM Clustering \n 47 Gene Signature Verhaak \n 42% Variance Explained") + geom_point(shape=19) + labs(y = "PC1", x = "PC2") + scale_color_manual(breaks = c("1", "2", "3","4"),values=  c("green","red","blue","orange"))



pdf('DPMMSignature.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y.order.2), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("green","red","blue","orange")[c.final.order.2], labRow = colnames(Y.order.2), labCol = NA, main = ' \n Training Set \n 47 Gene DPMM signature', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Best","Good Moderate","Bad moderate","Worst"),fill = c("Green","Blue","Orange","Red"), cex = 0.4)
p5
p1
dev.off()




############## Similar Plots for Testing Data Set ######################
########################################################################


rank <- matrix(0, nrow = N.new, ncol =Nps)

for (j in 1:Nps){
  
  Y.scaled <- matrix(0, nrow = N.new, ncol =D)
  for ( v in 1:4){
    clust <- which(c.matrix.new[,j] == v)
    Y.scaled[clust,1:D] <- scale(Y.new[clust,1:D], center = TRUE, scale = TRUE)
  }
  
  
  for ( i in 1:N.new){
    rank[i,j] <- dMVN(as.vector(t(Y.new[i,1:D])), mean = mu.list[[j]][c.sbc.new[i],1:D], Q= S.list[[j]][c.sbc.new[i],1:D,1:D], log = TRUE) - dMVN(as.vector(t(Y.new[i,1:D])), mean = mu.list[[j]][1,1:D], Q= S.list[[j]][1,1:D,1:D], log = TRUE) +  dnorm(x = time.new[i], mean = beta0.list[[j]][c.sbc.new[i]] + betahat.list[[j]][c.sbc.new[i],1:D] %*% as.vector(t(Y.scaled[i,1:D])), sd = sqrt(sigma2.list[[j]][c.sbc.new[i]]), log =TRUE) -  dnorm(x = time.new[i], mean = beta0.list[[j]][1] + betahat.list[[j]][1,1:D] %*% as.vector(t(Y.scaled[i,1:D])), sd = sqrt(sigma2.list[[j]][1]), log =TRUE) 
  }
}

avg.rank <- apply(rank,1,mean)
order.zo <- range01(avg.rank)

order.train <- sort(order.zo,index.return = TRUE, decreasing = TRUE)
Y.order.new <- Y.new[order.train$ix,]
c.final.order.new <- c.sbc.new[order.train$ix]


######## Reordering Again ################
order.2.new <- c(which(c.final.order.new==1),which(c.final.order.new==3),which(c.final.order.new==4),which(c.final.order.new==2))
Y.order.2.new <- Y.order.new[order.2.new,]
c.final.order.2.new <- c.final.order.new[order.2.new]


surv.ob <- Surv(exp(time.new),censoring.new)
Classes <- c.sbc.new
Classes <- as.factor(Classes)
surv.fit <- survfit(surv.ob ~ Classes)
p5 <- ggsurv(surv.fit, main = " Testing Data Set \n Kaplan Meier Estimators \n pvalue 3e -02", surv.col = c("green","red","blue","orange"), cens.col ="Blue")  


############ Generating some Plots ##########################
pc <- prcomp(Y.new)
pc.pred <- predict(pc,newdata = Y.new)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], color = as.factor(c.sbc.new) )) + ggtitle(" DPMM Clustering\n Testing Data \n 47 Gene Signature Verhaak \n 44% Variance Explained") + geom_point(shape=19) + labs(y = "PC1", x = "PC2") + scale_color_manual(breaks = c("1", "2", "3","4"),values=  c("green","red","blue","orange"))



pdf('TestingDPMMSignature.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y.order.2.new), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("green","red","blue","orange")[c.final.order.2.new], labRow = colnames(Y.order.2.new), labCol = NA, main = ' \n Testing Set \n 47 Gene DPMM signature', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Best","Good Moderate","Bad moderate","Worst"),fill = c("Green","Blue","Orange","Red"), cex = 0.4)
p5
p1
dev.off()

#########################################################################################################################################################
########################################################### Feature Selection #############################################################################
###########################################################################################################################################################

pdf("Heatmap_FeatureSelection_VerhaakDataSet.pdf")
hmcols<-colorRampPalette(c("white","black"))(128)
heatmap.2(t(as.matrix(heatmapdata[c(3,4),])),dendrogram="none", col =hmcols, margins=c(6,10), main = "Posterior prob.for Selection \n Verhaak Data Set ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE, labCol = c("Best","Worst","Moderate_Good","Moderate_Bad"))
dev.off()



####################################################################################################################################################
####################################################################################################################################################
########################### MAKING FUNCTIONAL /BIOLOGICAL SENSE OF THE DPMM SIGNATURE ##############################################################
####################################################################################################################################################
####################################################################################################################################################

###### Loading the Somatic Mutation Data Set ####################
load("/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/somaticMutations.rda")

######################################################################################################################
####### Apparently the Verhaak Samples DO NOT have Somatic Data associated with them ##################################
#####################################################################################################################

########## Look for Relationship between Our classes and Verhaak Classes ########################
matrix.assoc <- rbind(c(0,0,0,10), c(0,2,2,1), c(7,12,1,5), c(19,15,12,12))
fisher.test(matrix.assoc, workspace = 10000000,hybrid =TRUE)

################################################################################################################################################
##################################################### OVER REPRESENTED GO and KEGG terms in the sigature #######################################
#################################################################################################################################################
library(GOstats)
library(GO.db)
library(KEGG.db)
library(Category)

### Creating a GenNames to GO mapping
library("org.Hs.eg.db")

######## Getting Enriched GO and KEGG terms enriched #########################

################# Getting Entrez IDS from Gene Names ###################

xx <- as.list(org.Hs.egALIAS2EG)
entrez.signature <- list(0)

for ( i in 1:length(signature.dpmm)){
  entrez.signature[[i]] <- xx[[signature.dpmm[i]]]   
}
entre <- unlist(entrez.signature)

names.universe <- rownames(Y.train.prelim)
entrez.universe <- list(0)

for ( i in 1:length(names.universe)){
  entrez.universe[[i]] <- xx[[names.universe[i]]]   
}
uni <- unlist(entrez.universe)


params.bp <- new("GOHyperGParams", geneIds=entre, universeGeneIds=uni, annotation = "org.Hs.eg.db",ontology="BP", pvalueCutoff= 0.5, conditional= TRUE, testDirection="over")
params.mf <- new("GOHyperGParams", geneIds=entre, universeGeneIds= uni, annotation = "org.Hs.eg.db",ontology="MF", pvalueCutoff= 0.5, conditional= TRUE, testDirection="over")
params.cc <- new("GOHyperGParams", geneIds= entre, universeGeneIds= uni, annotation = "org.Hs.eg.db",ontology="CC", pvalueCutoff= 0.5, conditional= TRUE, testDirection="over")

hgOver.bp <- hyperGTest(params.bp)
hgOver.mf <- hyperGTest(params.mf)
hgOver.cc <- hyperGTest(params.cc)

#### KEGG pathways enrichment #########################
params.kegg <- new("KEGGHyperGParams", geneIds= entre, universeGeneIds=uni, annotation = "org.Hs.eg" , pvalueCutoff= 0.5,  testDirection="over")
hgKEGG <-  hyperGTest(params.kegg)


############################################
######### Doing multiple Test Correction #################
adjusted.pvalue.bh <- p.adjust(pvalues(hgOver.bp), method ="BH")[1:40]
table.vv <- cbind(summary(hgOver.bp)[1:40,],adjusted.pvalue.bh)
write.csv(table.vv, file = 'Verhaak.BP.csv')

######### Doing multiple Test Correction #################
adjusted.pvalue.bh <- p.adjust(pvalues(hgOver.mf), method ="BH")[1:20]
table.vv <- cbind(summary(hgOver.cc)[1:20,],adjusted.pvalue.bh)
write.csv(table.vv, file = 'Verhaak.MF.csv')

adjusted.pvalue.bh <- p.adjust(pvalues(hgKEGG), method ="BH")[1:20]
table.vv <- cbind(summary(hgKEGG)[1:20,],adjusted.pvalue.bh)
write.csv(table.vv, file = 'Verhaak.KEGG.csv')


##################################################################################################################
###############
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

##### Visualization of the Results ###############################
#### USE Top 35 Molecular Features OBTAINED FROM THE MODEL ############
sum.diff12  <- c(0)
for ( i in 1:Nps){
  sum.diff12 <- sum.diff12 + abs(mu.list[[i]][1,] - mu.list[[i]][2,])  
}
sum.diff12.sc <- (1/Nps) * (sum.diff12)/(diag(W))
cluste12 <- range01(sum.diff12.sc)

sum.diff23  <- c(0)
for ( i in 1:Nps){
  sum.diff23 <- sum.diff23 + abs(mu.list[[i]][2,] - mu.list[[i]][3,])  
}
sum.diff23.sc <- (1/Nps) * (sum.diff23)/(diag(W))

cluste23 <- range01(sum.diff23.sc)

sum.diff34  <- c(0)
for ( i in 1:Nps){
  sum.diff34 <- sum.diff34 + abs(mu.list[[i]][3,] - mu.list[[i]][4,])  
}
sum.diff34.sc <- (1/Nps) * (sum.diff34)/(diag(W))

cluste34 <- range01(sum.diff34.sc)

sum.diff41  <- c(0)
for ( i in 1:Nps){
  sum.diff41 <- sum.diff34 + abs(mu.list[[i]][4,] - mu.list[[i]][1,])  
}
sum.diff41.sc <- (1/Nps) * (sum.diff41)/(diag(W))

cluste41 <- range01(sum.diff41.sc)

sum.diff13  <- c(0)
for ( i in 1:Nps){
  sum.diff13 <- sum.diff13 + abs(mu.list[[i]][1,] - mu.list[[i]][3,])  
}
sum.diff13.sc <- (1/Nps) * (sum.diff13)/(diag(W))

cluste13 <- range01(sum.diff13.sc)

sum.diff24  <- c(0)
for ( i in 1:Nps){
  sum.diff24 <- sum.diff24 + abs(mu.list[[i]][2,] - mu.list[[i]][4,])  
}
sum.diff24.sc <- (1/Nps) * (sum.diff24)/(diag(W))
cluste24 <- range01(sum.diff24.sc)

####################### Feature Signficance matrix

cluste <- cbind(cluste12,cluste23,cluste34, cluste41, cluste13, cluste24)
heatmapdata <- cluste
rownames(heatmapdata) <-  colnames(Y)
colnames(heatmapdata) <- c("1 v s2","2 vs 3","3 v s4","4 vs 1","1 vs 3","2 vs 4")
hmcols<-colorRampPalette(c("white","black"))(128)
pdf("FeatureImportance.pdf")
heatmap.2(heatmapdata , margins=c(6,10),col = hmcols, main = "SBC signature \n Feature Importance \n GBM I ", cexCol = 0.85, cexRow = 0.7, Rowv = TRUE, Colv = TRUE, trace = "none", xlab = "Clusters") 
dev.off()


########################### New Code ADDED 1 AUGUST 2019 ###################################
######################### Code Added to visualize the miRNA data set  ###################
####################################### LIMMA  #############################################

##### Load the miRNA data for Verhaak Samples
mirna 

##### Load the SBC Labels
pheno  <- as.data.frame(c.sbc) 
names(pheno) <- c('labels')
##### Create the design matrix with these labels
mm <- model.matrix( ~ labels, pheno)
##### You can also define a contrast matrix which only looks at cluster specific contrasts
cont.matrix <- makeContrasts(cluster1vscluster2 = cluster1 - cluster2, levels= mm)
###### Fit limma model
fit <- lmFit(mirna, design = mm)
####### Fit particular contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
###### Use eBayes Normalization
fit3 <- eBayes(fit2)
###### Look for DE table
DEG_table <- topTable(fit3, adjust="BH", coef='cluster1vscluster2', number = Inf)


############################ Here for the difference between all the clusters ############
##### Load the miRNA data for Verhaak Samples
mirna 

##### Load the SBC Labels
pheno  <- as.data.frame(c.sbc) 
names(pheno) <- c('labels')
##### Create the design matrix with these labels
mm <- model.matrix( ~ labels, pheno)
##### You can also define a contrast matrix which only looks at cluster specific contrasts
cont.matrix <- makeContrasts(cluster1vscluster2 = cluster1 - cluster2, levels= mm)
###### Fit limma model
fit <- lmFit(mirna, design = mm)
####### Fit particular contrasts
#fit2 <- contrasts.fit(fit, cont.matrix)
###### Use eBayes Normalization
fit3 <- eBayes(fit)
###### Look for DE table
DEG_table <- topTable(fit3, adjust="BH", coef='cluster1vscluster2', number = Inf)


