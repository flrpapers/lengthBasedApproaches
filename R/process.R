library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(FLBRP)

library(pROC)

if ("Windows"%in%Sys.info()){
  setwd("p:/papers/inPrep/lengthMethods")}else{
    setwd("~/Desktop/inPrep/lbm")}


load("lbm/data/om/design.RData")

## LBIs #############################################################################################
## Check  that ran
ran=system2("ls",args="data/runs",stdout=TRUE)
ran=ran[grep("lbi",ran)]
ran=as.numeric(substr(ran,5,nchar(ran)-6))
ran=sort(ran[!is.na(ran)])

lbi=mdply(data.frame(.id=ran), function(.id){
  load(paste(file.path("data/om",  "om"),   .id,"RData",sep="."))
  load(paste(file.path("data/runs","lbi"),  .id,"RData",sep="."))
  load(paste(file.path("data/runs","haupt"),.id,"RData",sep="."))
  
  lbi=merge(lbi,model.frame(FLQuants("FFmsy"=fbar(om)/refpts(eq)["msy","harvest"],
                                     "BBmsy"=ssb( om)/refpts(eq)["msy","ssb"]),drop=T))
  lbi=merge(lbi,model.frame(FLQuants(z=z),drop=T))
 
  load(paste(file.path("data/runs","haupt.lmax"),.id,"RData",sep="."))
  lbi=merge(lbi,model.frame(FLQuants(z.lmax=z),drop=T))
 
  load(paste(file.path("data/runs","bh.z"),.id,"RData",sep="."))
  lbi=merge(lbi,model.frame(FLQuants(z.bh=z.bh),drop=T))
  
  lbi=transform(lbi,
                l95   = (l95/linf)/0.8,
                l25   = l25/l50,
                lmax5 = (lmax5/linf)/0.8,
                lmean = lmean/(2/3*linf),
                lbar  = lbar/l50,
                lmaxy = lmaxy/(linf*2/3),
                lc    = lc/l50,
                pmega = pmega/0.3)
  
  lbi})
save(lbi,file=file.path("data/results","lbi.RData"))

#### LBSPR #######################################################################################################
ran=system2("ls",args="lbm/data/runs",stdout=TRUE)
ran=ran[grep("lbspr",ran)]
ran=as.numeric(substr(ran,7,nchar(ran)-6))
ran=sort(ran[!is.na(ran)])

lbspr=mdply(data.frame(.id=ran), function(.id){
 load(paste(file.path("lbm/data/runs/","lbspr"),  .id,"RData",sep="."))
 lb})
save(lbspr,file=file.path("lbm/data/results","lbspr.RData"))

#### LBSPR Correct Defaults ######################################################################################
ran=system2("ls",args="lbm/data/runs",stdout=TRUE)
ran=ran[grep("lbsprUpdateDefaults",ran)]
ran=as.numeric(substr(ran,21,nchar(ran)-6))
ran=sort(ran[!is.na(ran)])

lbspr2=mdply(data.frame(.id=ran), function(.id){
  load(paste(file.path("lbm/data/runs","lbsprUpdateDefaults"),  .id,"RData",sep="."))
  lb})
save(lbspr2,file=file.path("lbm/data/results","lbspr2.RData"))


#### LIME ########################################################################################################
lime=mdply(expand.grid(.id=seq(90),iYr=seq(60,118,2)), function(.id,iYr){
  load(paste(file.path("lbm/data/runs","lime"),  .id,iYr,"RData",sep="."))
 
  dat=lime[as.numeric(lime$year)>=iYr,]
  dat})[,-2]
save(lime,file=file.path("lbm/data/results","lime.RData"))
