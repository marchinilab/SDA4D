## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
# #Load the library SDA4D
# #requires devtools
#install.packages('devtools')
#devtools::install_github('https://github.com/marchinilab/SDA4D')
# #load library
library(SDA4D)

## ----generatedata--------------------------------------------------------
set.seed(42)
generatedData <- generate_example_data()
 
lapply(generatedData,dim)

## ----runmethod,results=FALSE---------------------------------------------
 res<-RunSDA4D(data_tensor = generatedData$dataTensor,
               dimn_vector = c(200,500,16,3),
               max_iters = 50,
               num_components = 8,
               stopping=FALSE)

## ----plotELBO,fig.width=8,fig.height=6,echo=FALSE------------------------
#res<-list(Neg_FE=exp(rnorm(100)))
 # plot(-log(-res$Neg_FE[-1]),col='red',pch=4,xlab='Iteration',main='log(ELBO) over 100 iterations')
 # lines(-log(-res$Neg_FE[-1]),col='red')

## ----plotELBOqplot,echo=FALSE,fig.width=6,fig.height=4-------------------
library(ggplot2)
startfrom=5
qplot(x=c(startfrom:length(res$ELBO)),
      y=res$ELBO[-c(1:(startfrom-1))],
      geom=c("point", "line")
      )+
    ylab('ELBO')+
    xlab('Iteration')+
    ggtitle('ELBO over 50 iterations')

## ----outputnames---------------------------------------------------------
names(res)
res$maximumiteration
length(res$ELBO)

## ----outputnamesanddims--------------------------------------------------
lapply(res$A,dim)

lapply(res$WS,dim)

identical(res$WS$mom1,res$WS$m*res$WS$gamma)

