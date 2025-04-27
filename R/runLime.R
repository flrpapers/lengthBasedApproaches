library(LIME)
library(FLCore)
library(ggplotFL)
library(FLife)
library(mydas)
library(popbio)
  
library(plyr)
library(dplyr)
library(reshape)

library(ggpubr)
  
library(LIME)
  
library(foreach)
library(doParallel)
cl=makePSOCKcluster(7)
registerDoParallel(cl)

if ("Windows"%in%Sys.info()){
  setwd("p:/papers/inPrep/lengthMethods")}else{
  setwd("/home/laurie/pCloudDrive/papers/inPrep/lengthMethods")}
dirMy=file.path(getwd(),"lbm")

load("lbm/data/om/design.RData")


## Wrapper functions
lime_lh<-function(x,selex_type="flat",CVlen=0.2,m=NULL){
  
  create_lh_list(
    #list(
    # growth
    linf  =unlist(c(x["linf"])),
    vbk   =unlist(c(x["k"])),
    t0    =unlist(c(x["t0"])),
    lwa   =unlist(c(x["a"])),
    lwb   =unlist(c(x["b"])),
    AgeMax=unlist(c(vonB(params=x,length=x["linf"]*0.9))),
    
    # M
    M   =ifelse(is.null(m), mean(gislason(FLQuant(2:12),as(x,"FLPar"))),mean(gislason(FLQuant(40),as(x,"FLPar")))),
    
    # mat
    maturity_input="length",
    M50  =unlist(c(x["l50"])),
    
    # SRR
    h  =unlist(c(x["s"])),
    
    # selex
    selex_input="length",
    selex_type =selex_type,
    S50    =unlist(c(x["a50"])),
    S95    =unlist(c(x["a50"]+x["ato95"])),
    dome_sd=2,
    #
    CVlen=CVlen)}

lens<-function(x){
  fpp=cast(x,year~length,value="data")
  fpp=ddply(fpp, .(year), function(x) {x[is.na(x)]=0;x})
  
  fp=as.matrix(fpp[,-1],rownames=T)
  names(dimnames(fp))=c("year","length")
  dimnames(fp)$year=fpp[,1]
  fp}

runLime<-function(lfd,lh, mFlag=NULL){

  lh =lime_lh(lh,mFlag)
  
  res=foreach(iter=dimnames(lfd)$iter, 
        .combine=rbind,
        .multicombine=TRUE,
        .export=c("lh","lfd"),
        .packages=c("LIME","FLCore","data.table","FLife","plyr")) %dopar% {
                
    ln    =matrix(iter(lfd,iter),dim(lfd)[2],dim(lfd)[1],dimnames=dimnames(lfd)[2:1])
    inputs=create_inputs(lh=lh, list(years=dimnames(ln)$year, LF=ln))

    res=try(run_LIME(tempdir(), 
               input=inputs, data_avail=c("LC"), 
               vals_selex_ft=inputs$S_fl, 
               est_selex_f=FALSE, 
               C_type=0,
               derive_quants=TRUE,
               newtonsteps=FALSE)[c("Report","Derived")])

        if (!("try-error"%in%is(res))){
            res=try(data.frame("iter" =iter,
                       "year" =dimnames(ln)$year,
                       "F"    =res[["Report"]]$F_y,
                       "FFmsy"=res[["Derived"]]$FFmsy))
          
          if (("try-error"%in%is(res))) return(NULL) else return(res)
          }
        }
  res}

## Run
for (i in seq(dim(design)[1])){
   load(paste("lbm/data/om/lfd",i,"RData",sep="."))
  
   lh=lhs[[design[i,"Species"]]]
  
   for (year in seq(60,118,2)){
     print(year)
     lime=runLime(lfd[,ac(year+c(-1,0,1))],lh)
     save(lime,file=paste("lbm/data/runs/lime",i,year,"RData",sep="."))
     }

   #iSpp=rep(rep(1:5,each=6),3)[i]
   #lh=lhs[[iSpp]]
   #limeAll    =runLime(lfd[,ac(60:120)],lh)
   #save(limeAll,file=paste("/home/laurence-kell/Desktop/papers/COM3/R/lime/limeAll",i,"RData",sep="."))
}

design=cbind(.id=seq(90),design)

## Defaults ####################################################################
lime=mdply(expand.grid(.id=seq(90), year=seq(60,118,2)), function(.id,year){
   fl = paste("lbm/data/runs/lime",.id,year,"RData",sep=".")
   
   if (!file.exists(fl)) {
     load(paste("lbm/data/om/lfd",.id,"RData",sep="."))
     lh  =lhs[[design[.id,"Species"]]]
     lime=runLime(lfd[,ac(year+c(-1,0,1))],lh)
     save(lime,file=paste("lbm/data/runs/lime",.id,year,"RData",sep="."))
     return(NULL)}
   
   load(fl)
   lime$year=an(lime$year)
   cbind(design[.id,],lime[lime$year>=year,])
   })

save(lime,file="lbm/data/results/lime.RData")


## Correct Defaults ####################################################################
lime2=mdply(expand.grid(.id=seq(90), year=seq(60,118,2)), function(.id,year){
  fl = paste("lbm/data/runs/limeCorrectDefaults",.id,year,"RData",sep=".")
  
  if (!file.exists(fl)) {
    load(paste("lbm/data/om/lfd",.id,"RData",sep="."))
    lh  =lhs[[design[.id,"Species"]]]
    lh=lhs[[design[i,"Species"]]]
    lh["s"]   =design[i,"s"]
    lh["sel3"]=design[i,"sel3"]
      
    lime=runLime(lfd[,ac(year+c(-1,0,1))],lh, mFlag=design[i,"m"]=="low m")
    save(lime,file=paste("lbm/data/runs/limeCorrectDefaults",.id,year,"RData",sep="."))
    return(NULL)}
  
  load(fl)
  lime$year=an(lime$year)
  cbind(design[.id,],lime[lime$year>=year,])
})

save(lime,file="lbm/data/results/limeCorrectDefaults.RData")

ggboxplot(lime,y="F",x="year")+
 facet_grid(Species~.,scale="free")+
 scale_x_discrete(breaks=seq(60,120,10))

lm=as.FLQuant(transmute(subset(lime,as.numeric(as.character(year))>=as.numeric(as.character(yr))&.id==1),#
                        data=F,year=year,iter=iter))
plot(im)

ggboxplot(subset(lime,year==yr),y="FFmsy",x="year")+
 facet_grid(spp~.,scale="free")+
 geom_hline(aes(yintercept=1),col="red")+
 scale_x_discrete(breaks=seq(60,120,10))

l80=subset(lime,year==80)
l80=merge(lime,l80[,c(".id","iter","F","FFmsy")],by=c(".id","iter"))
l80=transform(l80,F=F.x/F.y,
                   FFmsy=FFmsy.x/FFmsy.y)
 
ggboxplot(subset(l80,year==yr),y="F",x="year")+
 facet_grid(spp~.,scale="free")+
 geom_hline(aes(yintercept=1),col="red")+
 scale_x_discrete(breaks=seq(60,120,10))
 
load("/home/laurence-kell/Desktop/papers/COM3/R/lime/limeAll.2.RData")
limeAll$year=factor(limeAll$year,levels=unique(sort(as.numeric(limeAll$year))))

ggplot(limeAll)+
   geom_boxplot(aes(year,F))+
   scale_x_discrete(breaks=seq(60,120,10))
 
