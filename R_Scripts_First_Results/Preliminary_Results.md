
---

title: "Notebook for the SBC on the GBM dataset from the Firebrowse Database(TCGA)"
author: "Camila Duitama (Based on Ashar Ahmad's code for the SBC model)"
date: "May 2019"
output:
  html_document:
    df_print: paged
    toc: yes
    toc_depth: '2'
  pdf_document:
    highlight: zenburn
    df_print: kable
    fig_caption: yes
    fig_height: 6
    fig_width: 10
    number_sections: yes
    toc: yes
    toc_depth: '2'
    
---

#Loading of libraries

Libraries needed to run the coded must be loaded first

#Loading/Splitting training and testing sets

Data is separated into training and testing


##Data visualization

Initial visualization of censoring, survival times and mRNAArray Data of the training and testing set

```r
table(censoring) %>%  
  kable(caption = "Censoring frequency for the training set") %>%
  kable_styling(bootstrap_options = c("striped", "hover","responsive"),full_width = T, position = "center")
```

<table class="table table-striped table-hover table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>Censoring frequency for the training set</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> censoring </th>
   <th style="text-align:right;"> Freq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 77 </td>
  </tr>
</tbody>
</table>


```r
table(censoring.new) %>%
  kable(caption = "Censoring frequency for the testing set") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed","responsive"),full_width = T, position = "center")
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>Censoring frequency for the testing set</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> censoring.new </th>
   <th style="text-align:right;"> Freq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 77 </td>
  </tr>
</tbody>
</table>



```r
  Y.pre.train[1:5,1:5] %>%
  kable(caption = "mRNA MicroArray data for the Training set " ) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed","responsive"),full_width = T, position = "center")
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>mRNA MicroArray data for the Training set </caption>
 <thead>
  <tr>
   <th style="text-align:right;"> AACS </th>
   <th style="text-align:right;"> FSTL1 </th>
   <th style="text-align:right;"> ELMO2 </th>
   <th style="text-align:right;"> CREB3L1 </th>
   <th style="text-align:right;"> RPS11 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 6.539245 </td>
   <td style="text-align:right;"> 9.794400 </td>
   <td style="text-align:right;"> 6.213981 </td>
   <td style="text-align:right;"> 4.836276 </td>
   <td style="text-align:right;"> 10.81124 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 7.186891 </td>
   <td style="text-align:right;"> 4.945053 </td>
   <td style="text-align:right;"> 5.230444 </td>
   <td style="text-align:right;"> 5.818606 </td>
   <td style="text-align:right;"> 10.47730 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 7.675038 </td>
   <td style="text-align:right;"> 10.840095 </td>
   <td style="text-align:right;"> 6.620676 </td>
   <td style="text-align:right;"> 5.333213 </td>
   <td style="text-align:right;"> 10.63727 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 7.996010 </td>
   <td style="text-align:right;"> 8.931571 </td>
   <td style="text-align:right;"> 7.552416 </td>
   <td style="text-align:right;"> 6.087341 </td>
   <td style="text-align:right;"> 11.00153 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 8.355122 </td>
   <td style="text-align:right;"> 4.240622 </td>
   <td style="text-align:right;"> 6.707334 </td>
   <td style="text-align:right;"> 4.865492 </td>
   <td style="text-align:right;"> 10.68588 </td>
  </tr>
</tbody>
</table>


```r
  Y.pre.test[1:5,1:5] %>%
  kable(caption = "mRNA MicroArray data for the Testing set " ) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed","responsive"),full_width = F, position = "center")
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>mRNA MicroArray data for the Testing set </caption>
 <thead>
  <tr>
   <th style="text-align:right;"> AACS </th>
   <th style="text-align:right;"> FSTL1 </th>
   <th style="text-align:right;"> ELMO2 </th>
   <th style="text-align:right;"> CREB3L1 </th>
   <th style="text-align:right;"> RPS11 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 5.936346 </td>
   <td style="text-align:right;"> 10.623493 </td>
   <td style="text-align:right;"> 6.996271 </td>
   <td style="text-align:right;"> 4.883836 </td>
   <td style="text-align:right;"> 10.51334 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5.665690 </td>
   <td style="text-align:right;"> 10.109494 </td>
   <td style="text-align:right;"> 6.795128 </td>
   <td style="text-align:right;"> 4.862345 </td>
   <td style="text-align:right;"> 10.37863 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5.714858 </td>
   <td style="text-align:right;"> 9.320218 </td>
   <td style="text-align:right;"> 6.893214 </td>
   <td style="text-align:right;"> 4.227045 </td>
   <td style="text-align:right;"> 10.29740 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 6.586214 </td>
   <td style="text-align:right;"> 10.458615 </td>
   <td style="text-align:right;"> 7.077229 </td>
   <td style="text-align:right;"> 5.274177 </td>
   <td style="text-align:right;"> 11.35884 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 8.967289 </td>
   <td style="text-align:right;"> 7.380626 </td>
   <td style="text-align:right;"> 7.504833 </td>
   <td style="text-align:right;"> 4.365866 </td>
   <td style="text-align:right;"> 11.23741 </td>
  </tr>
</tbody>
</table>


#Signature calculation

This chunk calculates  SBC gene signature
It's based on the idea of Univariate testing of Survival Data features
NOTE: For now the rows where there is no censoring data will be removed from ALL datasets, until we find an imputation method.



The SBC signature on the dataset looks like this


```r
 # to_plot<-data.frame(Genes=signature.sbc)
 #  to_plot%>%
 #    kable(caption = "SBC Signature") %>%
 #    kable_styling(bootstrap_options = c("striped", "hover", "condensed","responsive"),full_width = F, position = "center")%>%
 #    scroll_box(width = "500px", height = "200px")
```




```r
Label<- c()
for (i in 1:length(Clinical_TrainingSet$days_to_death)){
  if (i<=n) {
  Label<-append(Label,"Training")
  }else{
  Label<-append(Label,"Testing")
  }}
timedf<-data.frame(logtime=log(Clinical_TrainingSet$days_to_death),set=Label)
mu <- ddply(timedf, "Label", summarise, grp.mean=mean(logtime))
ggplot(data =timedf,aes(x=logtime,color=Label,fill=Label)) +
  geom_histogram(bins = 15, alpha=0.5, position="dodge")+ 
  labs(y="Counts", x = "log(time)")+ ggtitle("Log(time) for training and testing data")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Label),linetype="dashed")
```

![plot of chunk Histogram of log(time) for both training and testing set](figure/Histogram of log(time) for both training and testing set-1.png)


#Train the model

##Preparation of the data needed for the training of the model

```r
########################################################################
###### We prepare the Data Structures Needed for the Running of the SBC ####
#############################################################################

######## Verhaak Signature ###############
Verhaak_gene_signature <-read_csv("/Volumes/GoogleDrive/My Drive/Documents/Life Science Informatics Master/Master thesis/Data/Verhaak_Data/Verhaak_GS_according_to_Training_Set_Split1.csv", 
                                  col_names = FALSE, col_types = cols(X1 = col_skip()))
Verhaak_gene_signature<-Verhaak_gene_signature[-indexes,]
signature.vk<-Verhaak_gene_signature$X2
######## Verhaak Labels #################
labels.vk<-Clinical_TrainingSet$Subtype

#######Getting signature matrix on training and testing######
Y <- data.matrix(Y.pre.train[,as.character(signature.sbc)])
Y.new <- data.matrix(Y.pre.test[,as.character(signature.sbc)])
smod.new <- Surv(time.new, censoring.new)
c.true <- as.factor(labels.vk[1:n])
c.true <- unclass(c.true)
c.true.new<-as.factor(labels.vk[-(1:n)])
c.true.new<-unclass(c.true.new)
```


##Training results




##Plots of results from training

```r
###############################
### Some plots and analysis ####
#################################
#### Generating some plots and results ON Training data ###
logrank <- survdiff(smod ~ c.final)
pval<-1 - pchisq(unlist(logrank)$chisq,df =3)
surv.fit <- survfit(smod ~ c.final)
p5 <- ggsurv(surv.fit, plot.cens=FALSE,main = " DPMM \n Kaplan Meier Estimators \n Verhaak Cancer Data Set")+ ggplot2::guides(linetype = FALSE)+labs(subtitle = paste("P-value:",toString(round(pval,digits = 5))))
p5
```

![plot of chunk Survival curves from training](figure/Survival curves from training-1.png)


```r
############ Generating some Plots ##########################
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final))) + ggtitle(" SBC Clustering \n VerhaakCancer Data Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 
p1
```

![plot of chunk PCA from training](figure/PCA from training-1.png)
#Testing

##Plot results from testing

```r
logrank.new <- survdiff(smod.new ~ c.sbc.new)
#df= Degrees of freedom should be number of clusters-1
pval.new<-1 - pchisq(unlist(logrank.new)$chisq,df =length(logrank.new$n)-1)
surv.fit <- survfit(smod.new ~ c.sbc.new)
p5.new <- ggsurv(surv.fit, plot.cens=FALSE,main = " DPMM \n Kaplan Meier Estimators \n Verhaak Cancer Data Set")+ ggplot2::guides(linetype = FALSE)+labs(subtitle = paste("P-value:",toString(round(pval.new,digits = 5))))
p5.new
```

![plot of chunk Survival curves from testing](figure/Survival curves from testing-1.png)


```r
pc <- prcomp(Y.new)
pc.pred <- predict(pc,newdata = Y.new)
p3.new <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.sbc.new))) + ggtitle(" SBC Clustering \n Verhaak Cancer Test Data Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 
p3.new
```

![plot of chunk PCA from testing](figure/PCA from testing-1.png)

```r
#############
```
