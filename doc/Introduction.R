## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
#Load the library SDA4D
#library(SDA4D)

## ----generatedata--------------------------------------------------------
set.seed(42)
# generatedData <- generate_example_data()
# 
# dim(generatedData$dataTensor)
# 
# lapply(generateData,dim)

## ----runmethod-----------------------------------------------------------
#res<-RunSDA4D(generatedData$dataTensor,dimn_vector = c(),maxiters = 100,stopping=FALSE)

## ----plotELBO,fig.width=8,fig.height=6,echo=FALSE------------------------
res<-list(Neg_FE=exp(rnorm(100)))
# plot(log(res$Neg_FE[-1]),col='red',pch=4,xlab='Iteration',main='log(ELBO) over 100 iterations')
# lines(log(res$Neg_FE[-1]),col='red')

## ----plotELBOqplot,echo=FALSE,fig.width=6,fig.height=4-------------------
library(ggplot2)
qplot(x=c(1:length(res$Neg_FE[-1])),y=res$Neg_FE[-1],geom=c("point", "line"))+xlab('Iteration')+ggtitle('log(ELBO) over 100 iterations')

## ----outputnames---------------------------------------------------------
names(res)

