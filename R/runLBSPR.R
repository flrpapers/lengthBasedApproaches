library(FLCore)
library(ggplotFL)
library(FLife)
library(mydas)
library(popbio)
  
library(plyr)
library(dplyr)
library(reshape)
  
library(LBSPR)
  
library(foreach)
library(doParallel)

source("popdyn.R")

cl=makePSOCKcluster(6)
registerDoParallel(cl)

if ("Windows"%in%Sys.info()){
  setwd("p:/papers/inPrep/lengthMethods")}else{
  setwd("/home/laurie/pCloudDrive/papers/inPrep/lengthMethods")}
dirMy=file.path(getwd(),"lbm")

load("lbm/data/om/design.RData")

priors=llply(lhs,popdyn)

ran=system2("ls",args="lbm/data/runs",stdout=TRUE)
ran=ran[grep("lbspr",ran)]
ran=as.numeric(substr(ran,7,nchar(ran)-6))
toRun=seq(dim(design)[1])[!seq(dim(design)[1])%in%ran]

foreach(i=30:90, #seq(dim(design)[1]), #toRun, 
        .combine=rbind,
        .multicombine=TRUE,
        .export=c("lhs","priors","design","dirMy"),
        .packages=c("LBSPR","FLCore","mydas","FLife","plyr","reshape")) %dopar% {
          
  load(file.path(dirMy,paste("data/om/lfd",i,"RData",sep=".")))
  source("lbspr.R")          
  
  lh   =lhs[[   design[i,"Species"]]]
  prior=priors[[design[i,"Species"]]]
  
  lb=mdply(data.frame(year=seq(60,120,1)), function(year) {
    lb=lbspr(lfd[,ac(year)],prior)
    model.frame(lb)})[,c(2,6:10)]
    
  save(lb,file=file.path(dirMy,paste("data/runs/lbspr",i,"RData",sep=".")))}


## Correct Defaults ############################################################

ran=system2("ls",args="lbm/data/runs",stdout=TRUE)
ran=ran[grep("lbsprCorrectDefaults",ran)]
ran=as.numeric(substr(ran,22,nchar(ran)-6))
toRun=seq(dim(design)[1])[!seq(dim(design)[1])%in%ran]

foreach(i=seq(dim(design)[1]),  
        .combine=rbind,
        .multicombine=TRUE,
        .export=c("lhs","design","dirMy"),
        .packages=c("LBSPR","FLCore","mydas","FLife","plyr","reshape","popbio")) %dopar% {

          lh=lhs[[design[i,"Species"]]]
          lh["s"]   =design[i,"s"]
          lh["sel3"]=design[i,"sel3"]

          if (design[i,"m"]=="low m")
            prior=popdyn(lh,eq=lhEql(lh,m=function(x,params) {
            length=wt2len(stock.wt(x),params)
            #sets M= M at Linf
            length=length%=%params["linf"]
            exp(params["m1"]%+%(params["m2"]%*%log(length))%+%(params["m3"]%*%log(params["linf"]))%+%log(params["k"]))}))
          else
            prior=popdyn(lh)
              
          load(file.path(dirMy,paste("data/om/lfd",i,"RData",sep=".")))
          source("lbspr.R")          
          
          lb=mdply(data.frame(year=seq(60,120,1)), function(year) {
            lb=lbspr(lfd[,ac(year)],prior)
            model.frame(lb)})[,c(2,6:10)]
          
          save(lb,file=file.path(dirMy,paste("data/runs/lbsprUpdateDefaults",i,"RData",sep=".")))}
stopCluster(cl)

