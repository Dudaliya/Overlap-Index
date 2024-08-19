# Overlap-Index
 With the t-test for two independent samples, we are not in the same territory as we are with the statistical induction from the experiment on the lady tasting tea. There we had a clear hypothesis and a clearer probabilistic approach to test the hypothesis. With the t-test, we depend heavily on the theoretical distributions. This may occasionally lead to some awkward intuitions.
https://www.researchgate.net/publication/359443514_Some_Practical_Considerations_about_the_t-Test_for_Two_Independent_Samples
R Code:
---
title: "Two Sample t-Test - Overlap Index"
author: "Darshi Arachige"
date: "18/03/2022"
output: html_document
---

The Overlap Index
Let us define a measure to assess the overlap between two sets of sample data and further assume that x1,…xn and y1…ym are two samples each of which ascendingly ordered. This index looks at the overlap of observations from two samples using hundred bins of equal size. Any grouping in bins can cause some loss of information. But binning helps us to make a form of probabilistic statement such as the probability of an observation falling in a bin is 0.01. In each bin, except for the first and last bins, there can be none of the observations or at least one observation.If the samples are from superpopulations that are infinite in size the empty bins simply represent the bins where there aren’t any observations in the realised samples.  


The algorithm used is as follows:  

1.	Find the minimum x1 and the maximum ym for the pooled sample.  

2.	Divide the difference between the maximum and minimum by 100.  

3.	Then using the interval calculated above determine the 100 bins that range from the minimum to the maximum.  

4.	Now, impute x1 + (ym - x1)/200 for x1 and ym - (ym - x1)/200   - x1)/200 for ym.  

5.	Assign each observation from each sample to one of those bins.  

6.	Calculate the Overlap Index as the range between the minimum and maximum bins that are overlapping.  

This can also be implemented as a faster calculation on the continuous scale:  

1.	The Overlap Index can then be defined as:  

(Minimum of (xn, ym) - Maximum of (x1, y1))/ (ym – x1)


```{r, warning=FALSE,error=FALSE,eval=FALSE,tidy=TRUE}
library(plyr)
library(dbplyr)
library(tidyverse)
library(readxl)
```


```
'an example of calling sampling function'
'sampling_out<-sampling(data,xlm=40,n=10,n1=10, n2=500,method=0,dname="Short Topic")'

 
sampling<-function(data,xlm,n,n1,n2,method,dname){
'data - a set of data with two columns- Column 1 with the sample identifier, 1 or 2 and Column 2 with the values'
'xlm - x axis limits for the histograms of simulated data'  
'n- size of the first sample - for simulations only'
'n1 - size of the second sample - for simulations only'
'n2 - Number of simulations required'
'method - 0 - for the continuous data 1 - for the binned data '
'dname- A short description of data'
'the output will provide the overlap for the observed data and the overlap at percentiles for the simulations '
  
  if(method == 0){
    original<-infer(data=data)
  } else
  {
    original<-inference(data=data)
  }
  m<-as.matrix(data)
  xval<-subset(m[,2],m[,1]==1)
  yval<-subset(m[,2],m[,1]==2)
  minm0<-apply(m, 2, min)[[2]]
  maxm0<-apply(m, 2, max)[[2]]
  hist(xval,freq=FALSE,xlim=c(minm0,maxm0),main = paste("Histogram of " , dname),breaks=20,xlab="Bin Value",col=rgb(0, 0, 1, 0.5))
  hist(yval,freq=FALSE,xlim=c(minm0,maxm0),add=T,breaks=20,col=rgb(0, 1, 0, 0.5))
  
  m1<-matrix(0,nrow=n2,ncol=2)
  m2<-matrix(0,nrow=n2,ncol=2)
  
  for(j in 1:n2){
    data<-as.matrix(data)
    s1<-rep(0,n)
    s2<-rep(0,n1)
    mn1<-mean(subset(data[,2],data[,1]==1))
    mn2<-mean(subset(data[,2],data[,1]==2))
    
    sd1<-((var(subset(data[,2],data[,1]==1)))^0.5)
    sd2<-((var(subset(data[,2],data[,1]==2)))^0.5)
    
    s1<-rnorm(n, mean= mn1, sd = sd1)
    s2<-rnorm(n1, mean= mn2, sd = sd2)
    
    dable<-rbind(cbind(1,(matrix(s1, ncol = 1))),cbind(2,(matrix(s2, ncol = 1))))
    if(method == 0){
      dble<-infer(data=dable)
    } else
    {
      dble<-inference(data=dable)
    }
    
    m1[j,1]<-mean(s1)
    m2[j,1]<-mean(s2)
    m1[j,2]<-dble$overlap
  }
  qtle<-quantile(m1[,2],probs=c(0.01,0.05,0.5,0.95,0.99))
  hist(m1[,1],freq=FALSE,xlim=c(0,xlm),main =paste( "Means - Simulated data \n", dname),breaks=20,xlab="Bin Value",col=rgb(0, 0, 1, 0.5))
  hist(m2[,1],freq=FALSE,xlim=c(0,xlm),add=T,breaks=20,col=rgb(0, 1, 0, 0.5))
  
  hist(m1[,2],freq=FALSE,main = paste("Overlaps - Simulated data \n",dname),breaks=10,xlab="Bin Value",col="red")
  list(original$overlap,sim_qtile=qtle)
 
}



infer <-function (data){
  m0<-as.matrix(data)
  
  
  {
    minm<-apply(m0, 2, min)[[2]]
    maxm<-apply(m0, 2, max)[[2]]
    intv<-(maxm-minm)
    
    
    g1<-subset(m0[,2],m0[,1]==1)
    g2<-subset(m0[,2],m0[,1]==2)
    ming<-max(min(g1),min(g2))
    maxg<-min(max(g1),max(g2))
    
    
    if(ming>=maxg){prob1<-0}else{ 
      
      prob1<-(maxg-ming)/intv
      
    }
    
  }
  list(overlap=prob1)
}


inference <-function (data){
  m0<-as.matrix(data)
  m<-m0
  
  {
    minm0<-apply(m0, 2, min)[[2]]
    maxm0<-apply(m0, 2, max)[[2]]
    intv<-(maxm0-minm0)/100
    
    rep1 <-rep(0:100)
    rep11<-rep1[-1]
    
    rep2<-minm0+rep1*intv
    rep21<-rep2[-101]
    
    rep3<-(rep2)
    rep31<-rep3[-1]
    
    
    mn<-length(m0[,1])
    row<-rep(1:mn)
    rown<-rep(row, each=mn)
    m011<-rep(m0[,1], mn)
    m0[,2][m0[,2]==minm0]<-minm0+intv/2
    m0[,2][m0[,2]==maxm0]<-maxm0-intv/2
    m012<-rep(m0[,2], mn)
    
    mdf_<-(data.frame(rown=rep(rep(1:100),each=mn),v0=rep(m011,100),v1=rep(m012,100)))
    mg1<-subset(mdf_,v0==1)
    mg2<-subset(mdf_,v0==2)
    mdf2<-data.frame(rown=rep11,v2=rep21,v3=rep31)
    bin_sync1<-left_join(mdf2,mg1, by= "rown")
    bin_sync2<-left_join(mdf2,mg2, by= "rown")
    
    mg11<-as.matrix(bin_sync1)
    mg22<-as.matrix(bin_sync2)
    
    
    myg1<-binning(mg11)
    myg2<-data.frame(myg1)
    myg3<-(subset(myg2,X6==1))
    myg12<-binning(mg22)
    myg22<-data.frame(myg12)
    myg32<-(subset(myg22,X6==1))
    counter<-data.frame(ddply(myg3,.(X1,X6),summarize,N=length(X1)/mn,minx5=min(X5),minx1=min(X1),minx2=min(X2),minx3=min(X3)))
    counter2<-data.frame(ddply(myg32,.(X1,X6),summarize,N=length(X1)/mn,minx5=min(X5),minx1=min(X1),minx2=min(X2),minx3=min(X3)))
    count_sync<-inner_join(counter,counter2, by= "X1")
    counter_min<-min(counter$X1)
    counter_max<-max(counter$X1)
    counter2_min<-min(counter2$X1)
    counter2_max<-max(counter2$X1)
    counter_dmin<-max(counter_min,counter2_min)
    counter_dmax<-min(counter_max,counter2_max)
    
    
    if(counter_dmax<=counter_dmin){prob1<-0}else{ 
      
      prob1<-(counter_dmax-counter_dmin+1)/100
      
    }
    
  }
  list(overlap=prob1)
}

binning<-function(m1){
  
  y<-matrix(0:0,length(m1[,1]),6)
  y[,1:5]<-  m1[,1:5]
  for (i in 1:length(m1[,1])){
    if( ( (m1[i,5]>=m1[i,2])&&(m1[i,5]< m1[i,3]) ) || ( (m1[i,5]>=m1[i,2])&&(m1[i,2]== m1[i,3]) )  ) {
      y[i,6]<-1
    }else{y[i,6]<-NA}
  }
  list(y)
} 

