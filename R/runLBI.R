#####################################################################################################
## Calculate LBIs                                                                                  ##
#####################################################################################################
library(ggplotFL)
library(FLBRP)
library(FLasher)
library(FLife)
library(mydas)
library(popbio)

library(plyr)
library(dplyr)
library(reshape)

library(spatstat)

library(doParallel)
library(foreach)

cl=makePSOCKcluster(7)
registerDoParallel(cl)

if ("Windows"%in%Sys.info()){
    setwd("p:/papers/inPrep/lengthMethods")}else{
    setwd("/home/laurie/pCloudDrive/papers/inPrep/lengthMethods")}
dirMy=file.path(getwd(),"lbm")

source("haupt.R")

load("lbm/data/om/design.RData")

ran=system2("ls",args="data/om",stdout=TRUE)
ran=ran[grep("lbi",ran)]
ran=as.numeric(substr(ran,5,nchar(ran)-6))
toRun=seq(dim(design)[1])[!seq(dim(design)[1])%in%ran]

foreach(i=seq(dim(design)[1]),             
          .combine=rbind,
          .multicombine=TRUE,
          .export=c("design","lhs","haupt","gislasonM","dirMy"),
          .packages=c("FLCore","FLife","plyr","mydas","reshape","spatstat")) %dopar% {
    
            
    source("oemLn.R")
            
    yrs=seq(55,120,1)
            
    par=lhs[[design[i,"Species"]]]
            
    load(file.path(dirMy,paste("data/om/om",i,"RData",sep=".")))
        
    ## lfd
    set.seed(6789)
    ak  =invAlk(par,cv=0.1)  
    lfd=lenSample(catch.n(om)[,ac(yrs)],ak,nsample=design[i,"nsample"])
    save(lfd,ak,file=file.path(dirMy,paste("data/om/lfd",i,"RData",sep=".")))
    
    ## Haupt Assessment
    z=haupt(lfd,par,par["sel1"]*0.9)
    m=gislasonM(par)
    save(z,m,file=file.path(dirMy,paste("data/runs/haupt",i,"RData",sep=".")))
   
    ## upper limit
    z=haupt(lfd,par,par["sel1"]*0.9,par["linf"])
    m=gislasonM(par)
    save(z,m,file=paste("/home/laurie/pCloudDrive/papers/inPrep/lengthMethods/lbm/data/runs/haupt.lmax",i,"RData",sep="."))
    
    
    ### LBI
    lbi =transform(subset(as.data.frame(lfd,drop=TRUE),data>0),
    wt  =c(par["a"])*len^c(par["b"]),
    lopt=c(2/3*par["linf"]))
    lbi =ddply(lbi, .(year,iter), with, lenInd(len,data,wt,lopt))
    lbi=cbind(lbi,linf=c(par["linf"]),k=c(par["k"]),l50=c(par["l50"]))
    save(lbi,file=paste("lbm/data/runs/lbi",i,"RData",sep="."))
    }

foreach(i=seq(dim(design)[1]),             
        .combine=rbind,
        .multicombine=TRUE,
        .export=c("design","lhs","haupt","gislasonM"),
        .packages=c("FLCore","FLife","plyr","mydas","reshape","spatstat")) %dopar% {
     
            yrs=seq(55,120,1)
            
            par=lhs[[design[i,"Species"]]]
            
            load(paste("/home/laurie/pCloudDrive/papers/inPrep/lengthMethods/lbm/data/om/lfd",i,"RData",sep="."))
            
            ## Haupt Assessment
            z=haupt(lfd,par,par["sel1"]*0.9,par["linf"]*.9)
            m=gislasonM(par)
            save(z,m,file=paste("/home/laurie/pCloudDrive/papers/inPrep/lengthMethods/lbm/data/runs/haupt.lmax",i,"RData",sep="."))
        }

foreach(i=seq(dim(design)[1]),             
        .combine=rbind,
        .multicombine=TRUE,
        .export=c("design","lhs","haupt","gislasonM"),
        .packages=c("FLCore","FLife","plyr","mydas","reshape","spatstat")) %dopar% {

          bhz<-function(lfd,par){
            lfd   =lfd[dimnames(ldf)[[1]]>=c(par["lc"]),]
            lmean=quantSums(lfd*an(ages(lfd)))%/%quantSums(lfd)
            
            par["k"]%*%(par["linf"]%-%lmean)%/%(lmean%-%par["lc"])}
          
          yrs=seq(55,120,1)
          
          load(paste("/home/laurie/pCloudDrive/papers/inPrep/lengthMethods/lbm/data/om/lfd",i,"RData",sep="."))
          
          par=lhs[[design[i,"Species"]]]
          par=rbind(FLPar("lc"=FLife:::vonB(par["sel1"],par)),lhs[[design[i,"Species"]]])
          
          z.bh=bhz(lfd,par)
          save(z.bh,file=paste("/home/laurie/pCloudDrive/papers/inPrep/lengthMethods/lbm/data/runs/bh.z",i,"RData",sep="."))
          
          ## Haupt Assessment
          z=haupt(lfd,par,par["sel1"]*0.9,par["linf"]*.9)
          m=gislasonM(par)
          save(z,m,file=paste("/home/laurie/pCloudDrive/papers/inPrep/lengthMethods/lbm/data/runs/haupt.lmax",i,"RData",sep="."))
        }
