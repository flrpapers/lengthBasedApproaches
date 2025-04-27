#####################################################################################################
## Condition OMs                                                                                   ##
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

library(doParallel)
library(foreach)

if ("Windows"%in%Sys.info()){
  setwd("p:/papers/inPrep/lengthMethods")}else{
  setwd("/home/laurie/pCloudDrive/papers/inPrep/lengthMethods")}

dirMy="papers/inPrep/lengthMethods/lbm"


### Conditioning Operating Models on life histories ###########################################################
spp  =c("Pollachius pollachius","Psetta maxima","Raja clavata",
       "Sprattus sprattus sprattus","Thunnus obesus")
stock=c("Pollack","Turbot","Ray","Sprat","Bigeye")

## Read in parameters from fishbase 
load(url("https://github.com//fishnets//fishnets//blob//master//data//fishbase-web//fishbase-web.RData?raw=True"))

## lh
lh=subset(fb,species%in%spp) 

names(lh)[c(14:17)]=c("l50","l50min","l50max","a50")
lh=lh[,c("species","linf","k","t0","a","b","a50","l50","l50min","l50max")]

lh[is.na(lh["l50"]),"l50"]=(lh[is.na(lh["l50"]),"l50min"]+lh[is.na(lh["l50"]),"l50max"])/2
lh=lh[,-(9:10)]

# Add extra parameters for pollack
lh[lh$t0>=0,"t0"]=NA
lh=rbind(lh[1,],lh)
lh[1,-c(1,7:8)]=c(84.6,0.19,-0.94,0.01017,2.98)
lh[1,"l50"]=0.72*lh[1,"linf"]^0.93

lh=transform(lh,l50linf=l50/linf)

fctr=ddply(lh,.(species),with,mean(linf))
fctr=as.character(fctr[order(fctr$V1),"species"])

lh$species=factor(as.character(lh$species),levels=fctr)

lhs=dlply(lh,.(species), with, {
  
  res=lhPar(FLPar(
    linf =mean(linf,na.rm=T),
    k    =mean(k,   na.rm=T),
    t0   =mean(t0,  na.rm=T),
    l50  =mean(l50, na.rm=T),
    a    =mean(a,   na.rm=T),
    b    =mean(b,   na.rm=T)))
  
  res=rbind(res,FLPar(lopt =2/3*mean(res["linf"])))
  # also Linf*(3/(3+M/K)) 
  res=rbind(res,FLPar(lmega=1.1*mean(res["lopt"])))
  res})

cor=cor(model.frame(lh)[,c("linf","k","l50","t0","a","b","l50linf")],
        use="pairwise.complete.obs")


## OM Scenarios #########################################################################
design=data.frame(s      =c( 0.7,        0.9,        0.7,        0.7,        0.7,     0.7),
                  sel1   =c( 0.0,        0.0,        0.0,        0.0,        0.0,     1.0),
                  sel3   =c( 5000,       5000,        50,       5000,       5000,    5000),
                  nsample=c( 250,        250,        250,        100,        250,     250),
                  m      =c(rep("gislason",4),"low M","gislason"))
design=mdply(expand.grid(Species=spp,CV=c("0.3","0.5","AR"),stringsAsFactors=FALSE), 
             function(Species,CV) design)

save(design,lh,lhs,cor,file="lbm/data/om/design.RData")

## Condition OMs ########################################################################
nits=100
set.seed(234)
srDev=FLQuants(NULL)
srDev[["0.3"]]=rlnoise(nits,FLQuant(0,dimnames=list(year=1:130)),0.3,0.0);set.seed(234)
srDev[["0.5"]]=rlnoise(nits,FLQuant(0,dimnames=list(year=1:130)),0.5,0.0);set.seed(234)
srDev[["AR"]] =rlnoise(nits,FLQuant(0,dimnames=list(year=1:130)),0.3,0.7)

set.seed(1234)
uDev=rlnoise(nits,FLQuant(0,dimnames=list(year=1:130)),0.2,0.0)

## Condition OMs ####################################################################################
runOM<-function(i,design,lhs,srDev,  
                f=FLQuant(c(rep(0.1,60),seq(0.1,2.0,length.out=40)[-40],
                          seq(2.0,0.7,length.out=11),rep(0.7,20)))){
  
  #modify parameters based on design                                                
  par        =lhs[[design[i,"Species"]]]
  par["s"]   =design[i,"s"]
  par["sel1"]=max(1,par["sel1"]-design[i,"sel1"])
  par["sel3"]=design[i,"sel3"]
  
  ## simulated stock
  if (design[i,"m"]=="gislason")
    eq=lhEql(par)
  else{
    gislasonFn<-function(x,params) {
      
      length=wt2len(stock.wt(x),params)
      #sets M= M at Linf
      length=length%=%params["linf"]
      exp(params["m1"]%+%(params["m2"]%*%log(length))%+%(params["m3"]%*%log(params["linf"]))%+%log(params["k"]))}
    
    eq=lhEql(par,m=function(x,params) gislasonFn(x,params))}
  
    fapexAge<-function(object){
      tmp=harvest(object)
      tmp[]=fapex(object)
      tmp=FLQuant(ages(tmp)[harvest(object)==tmp])
      tmp}
    
  ## historical exploitation pattern
  fbar(eq)=f%*%fbar(eq)[catch(eq)==max(catch(eq))]
  fbar(eq)=qmin(fbar(eq),refpts(eq)["crash","harvest"]*0.9)
  
  ## stochastic simulations
  om      =window(as(eq,"FLStock"))
  stock.n(om)[stock.n(om)<0]=0.0001
  catch.n(om)[catch.n(om)<0]=0.0001
  landings.n(om)[landings.n(om)<0]=0.0001
  om      =FLasher::ffwd(propagate(om,dims(srDev[[1]])$iter),
               fbar=propagate(fbar(eq)[,-1],dims(srDev[[1]])$iter),
               sr=eq,deviances=srDev[[design[i,"CV"]]])
  
  scl=refpts(eq)["msy","harvest"]
  range(eq)[c("minfbar","maxfbar")]=mean(fapexAge(om))
  range(om)[c("minfbar","maxfbar")]=mean(fapexAge(om))
  
  eq=brp(eq)
  fbar(eq)=fbar(eq)%*%refpts(eq)["msy","harvest"]/scl
  
  ## save 
  save(om,eq,file=paste("lbm/data/om/om",i,"RData",sep="."))
  
  design[i,]}

cl=makePSOCKcluster(4)
registerDoParallel(cl)

foreach(i=seq(dim(design)[1]),#toRun, 
               .combine=list,
               .multicombine=TRUE,
               .export=c("runOM","srDev","lhs","design","haupt","ffwd"),
               .packages=c("FLCore","ggplotFL","FLBRP","FLasher",
                           "FLife","mydas","JABBA","mpb",
                           "popbio",
                           "plyr","dplyr","reshape","GGally")) %dopar% {
                             
      runOM(i,lhs=lhs,design=design,srDev)}

## Check and rerun as required
ran=system2("ls",args="inPrep/lengthMethods/lbm/data",stdout=TRUE)
ran=ran[grep("om",ran)]
ran=as.numeric(substr(ran,4,nchar(ran)-6))
ran=sort(ran[!is.na(ran)])
toRun=seq(dim(design)[1])[!seq(dim(design)[1])%in%ran]

