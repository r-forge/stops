
library(stops)
data(kinshipdelta)


diss <- as.matrix(kinshipdelta)
set.seed(210485)

#############################
# Test all losses dispatch
structures <- c("clinearity")
optstrain <- stops(diss,loss="strain",verbose=3,structures=structures,upper=5,lower=0) #works
optstress <- stops(diss,loss="stress",verbose=3,structures=structures,upper=5,lower=1) #not works
optsmacofSym <- stops(diss,loss="smacofSym",verbose=3,structures=structures,upper=5,lower=0) #not works
optsmacofSph <- stops(diss,loss="smacofSphere",verbose=3,structures=structures,upper=5,lower=0) #not works
optsammon <- stops(diss,loss="sammon",verbose=3,structures=structures,upper=5,lower=0) #works
optsammon2 <- stops(diss,loss="sammon2",verbose=3,structures=structures,upper=5,lower=0) #not works
optelastic <- stops(diss,loss="elastic",verbose=3,structures=structures,upper=5,lower=0) #not works
optsstress <- stops(diss,loss="sstress",verbose=3,structures=structures,upper=5,lower=0) #works
optrstress <- stops(diss,loss="rstress",verbose=3,structures=structures,upper=5,lower=0) #works
optpowermds <- stops(diss,loss="powermds",verbose=3,structures=structures,upper=c(5,5),lower=c(0,0)) #works
optpowerstress <- stops(diss,loss="powerstress",verbose=3,structures=structures,upper=c(5,5,5),lower=c(0,0,0)) #works
optpowersammon <- stops(diss,loss="powersammon",verbose=3,structures=structures,upper=c(5,5,5),lower=c(0,0,0)) #works
optpowerelastic <- stops(diss,loss="powerelastic",verbose=3,structures=structures,upper=c(5,5,5),lower=c(0,0,0)) #works
optpowerstrain <- stops(diss,loss="powerstrain",verbose=3,structures=structures,upper=5,lower=0) #works
optisomapk <- stops(diss,loss="isomap",verbose=3,structures=structures,upper=10,lower=2) #works
optisomapeps <- stops(diss,loss="isomapeps",verbose=3,structures=structures,lower=50,upper=100) #works
optbc <- stops(diss,loss="bcstress",verbose=3,structures=structures,lower=c(0,0,0),upper=c(5,5,5)) #works
optlmds <- stops(diss,loss="lmds",theta=1,verbose=3,structures=structures,lower=c(1,0),upper=c(20,2)) #works

# print method
optstrain
optstress
optsmacofSym
optsammon
optsammon2
optelastic
optsstress
optrstress
optpowermds
optpowerstress
optpowersammon
optpowerelastic
optpowerstrain
optelastic
optisomapk
optisomapeps
optbc
optlmds

# plot method
plot(optstrain)
plot(optstress)
plot(optsmacofSym)
plot(optsammon)
plot(optsammon2)
plot(optelastic)
plot(optsstress)
plot(optrstress)
plot(optpowermds)
plot(optpowerstress)
plot(optpowersammon)
plot(optpowerelastic)
plot(optpowerstrain)
plot(optelastic)
plot(optisomapk)
plot(optisomapeps)
plot(optbc)
plot(optlmds)


###########
#test stoploss combination
verbose <- 3
stressweight <- 1
structures <- c("clinearity","cdependence")
type <- c("additive")
s1a <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose)
s1a
#test stoploss with clinearity multiplicative
s1m <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="multiplicative",verbose=verbose)
s1m


#########
# Test optimizers
sannres <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose,optimmethod="SANN")
aljres <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose,optimmethod="ALJ")
psores <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose,optimmethod="pso")
krigres <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose,optimmethod="Kriging")
tgpres <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose,optimmethod="tgp")

#########
# Test different weights
structures <- c("clinearity","cdependence")
halfres <- stoploss(diss,stressweight=1,structures=structures,strucweight=c(-0.5,-0.5),verbose=verbose,optimmethod="ALJ")
fullres <- stoploss(diss,stressweight=0.5,structures=structures,strucweight=c(-1,-10),verbose=verbose,optimmethod="ALJ")
all.equal(halfres$conf,fullres$conf)


######
# Test  different theta lengths
structures <- c("clinearity","cdependence")
#theta scalar
thetas <- stoploss(diss,stressweight=1,loss="powerstress",theta=1,structures=structures,strucweight=c(-0.5,-0.5),verbose=verbose,optimmethod="ALJ",upper=5,lower=0)
#theta vector of 2
thetav <- stoploss(diss,stressweight=1,loss="powerstress",theta=c(2,2),structures=structures,strucweight=c(-0.5,-0.5),verbose=verbose,optimmethod="ALJ",upper=5,lower=0)
#theta vector of 3
thetav <- stoploss(diss,stressweight=1,loss="powerstress",theta=c(2,2,4),structures=structures,strucweight=c(-0.5,-0.5),verbose=verbose,optimmethod="ALJ",upper=5,lower=0)





#test stop_smacofSym
sres <- stop_smacofSym(dis,structures=structures,strucpars=strucpars)
sres <- stop_smacofSym(dis,structures=structures,strucpars=strucpars,type="multiplicative")
sres <- stop_smacofSym(dis,theta=c(1,2,1),structures=structures,strucpars=strucpars,type="multiplicative")
sres <- stop_smacofSym(dis,theta=c(2,2,2),structures=structures,strucpars=strucpars,type="multiplicative")
sres <- stop_smacofSym(dis,theta=2,structures=structures,strucpars=strucpars,type="multiplicative")

#test stop_powerstress
pres <- stop_powerstress(dis,structures=structures,strucpars=strucpars)
pres <- stop_powerstress(dis,structures=structures,strucpars=strucpars,type="multiplicative")
pres <- stop_powerstress(dis,theta=c(1,2,1),structures=structures,strucpars=strucpars,type="multiplicative")
pres <- stop_powerstress(dis,theta=c(2,2,2),structures=structures,strucpars=strucpars,type="additive")
pres <- stop_powerstress(dis,theta=2,structures=structures,strucpars=strucpars,type="multiplicative")
pres <- stop_powerstress(dis,theta=c(2,2,-2),weightmat=dis,structures=structures,strucpars=strucpars,type="additive")


library(stops)
data(kinshipdelta)
resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="additive",strucweight=c(-0.5,-0.5))

resa1<-stops(kinshipdelta,structures=c("cclusteredness","clinearity"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="additive")

resa2<-stops(kinshipdelta,structures=c("cclusteredness","clinearity"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="additive",strucweight=c(0.5,0.5))

resm<-stops(kinshipdelta,structures=c("cclusteredness","clinearity"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="multiplicative")

resm<-stops(kinshipdelta,structures=c("cclusteredness","clinearity"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="multiplicative",strucweight=c(-0.5,-0.5))

resm1<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cmanifoldness"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="multiplicative",strucweight=c(0.5,0.5))


#all stops losses
library(stops)
data(kinshipdelta)

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="stress",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,-0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cdependence"),loss="powerstress",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cassociation","cmanifoldness","clinearity"),loss="powermds",verbose=3,strucpars=list(NULL,NULL,NULL),type="additive",strucweight=c(-0.5,-0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="strain",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="powerstrain",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="powersammon",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="rstress",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="sstress",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="powerelastic",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="sammon",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="smacofSphere",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","ccomplexity"),loss="sammon2",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,theta=c(1,1,1),structures=c("cfunctionality","clinearity"),loss="elastic",verbose=3,strucpars=c(list(NULL),list(NULL)),type="additive",strucweight=c(-0.5,0.5),lower=c(1,0.5,-2),upper=c(3,4,-2))
resa




#test for kapp=lambda=nu=1
dis <- smacof::kinshipdelta
kappa <- 1
lambda <- 1
nu <- 1
eps <- 10
minp <- 2
rango <- c(0,1.3) 
fit <- powerStressMin(dis,kappa=kappa,lambda=lambda,nu=nu)
c1<-cordillera(fit$conf,epsilon=eps,minpts=minp,rang=rango)

expect_equal(coploss(fit,epsilon=eps,minpts=minp,rang=rango,stressweight=0,cordweight=1),stoploss(fit,strucpars=list(list(epsilon=eps,minpts=minp,rang=rango)),structures="cclusteredness",stressweight=0,strucweight=-1))


#Bug form using the wrong list - huge cordillera - becuae this way q is set to 10 as we only say "eps" - epsilon is no problem?
resa<-stops(kinshipdelta,structures=c("cclusteredness"),loss="powermds",verbose=3,strucpars=list(eps=10,minpts=2,rang=c(0,1.3)),type="additive",strucweight=c(-1),stressweight=0)
resa
resa<-stops(kinshipdelta,structures=c("cclusteredness"),loss="powermds",verbose=3,strucpars=list(epsilon=10,minpts=2,rang=c(0,1.3)),type="additive",strucweight=c(-1),stressweight=0)
resa





##### Trying out all losses in pcops with the correct for theta and upper and lower
library(stops)
dis <- smacof::kinshipdelta #df
dis <- cops::matchphi #matrix
dis <- smacof::morse #dist object

strucs <- c("cfunctionality","cstringiness")
strucpars <- list(list(alpha=0.9,C=13),list(NULL))

dis <- as.matrix(dis)

#strain
p1 <- stops(dis,loss="strain",type="ratio",theta=1,lower=0,upper=5,verbose=3,structures=strucs,strucpars=strucpars)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#stress 
p2 <- stops(dis,loss="stress",type="ratio",theta=1,lower=0,upper=5,verbose=3,structures=strucs,strucpars=strucpars)
p2
summary(p2)
plot(p2,"confplot")
plot(p2,"Shepard")
plot(p2,"transplot") 
#sammon 
p3 <- stops(dis,loss="sammon",type="ratio",theta=2,lower=0,upper=5,verbose=3,structures=strucs,strucpars=strucpars)
p3
summary(p3)
plot(p3,"confplot") #issues?
plot(p3,"Shepard")
plot(p3,"transplot")
#sammon2
p4 <- stops(dis,loss="sammon2",type="ratio",theta=1,lower=0,upper=5,verbose=3,structures=strucs,strucpars=strucpars)
p4
summary(p4)
plot(p4,"confplot")
plot(p4,"Shepard")
#rstress
p5 <- stops(dis,loss="rstress",type="ratio",theta=1,lower=0,upper=5,verbose=3,structures=strucs,strucpars=strucpars)
p5
summary(p5)
plot(p5,"confplot")
plot(p5,"Shepard") #Wrong plot
plot(p5,"Shepard")
#elastic 
p6 <- stops(dis,loss="elastic",type="ratio",theta=1,lower=0,upper=2,verbose=3,structures=strucs,strucpars=strucpars) 
p6
summary(p6)
plot(p6,"confplot")
plot(p6,"Shepard")
#powerstrain #same as strain
p7 <- stops(dis,loss="powerstrain",type="ratio",theta=1,lower=0,upper=5,verbose=3,structures=strucs,strucpars=strucpars) 
p7
summary(p7)
plot(p7,"confplot")
plot(p7,"Shepard")
#sstress
p8 <- stops(dis,loss="sstress",type="ratio",theta=1,lower=0,upper=5,verbose=3,structures=strucs,strucpars=strucpars) 
p8
summary(p8)
plot(p8,"confplot")
plot(p8,"Shepard")
#powermds
p9 <- stops(dis,loss="powermds",type="ratio",theta=c(0,1),lower=c(0,0),upper=c(5,5),verbose=3,structures=strucs,strucpars=strucpars) 
p9
summary(p9)
plot(p9,"confplot")
plot(p9,"Shepard")
#powersammon
p10 <- stops(dis,loss="powersammon",type="ratio",theta=c(1,1),lower=0,upper=5,verbose=3,structures=strucs,strucpars=strucpars) 
p10
summary(p10)
plot(p10,"confplot")
plot(p10,"Shepard")
#powerelastic 
p11 <- stops(dis,loss="powerelastic",type="ratio",theta=c(0.0963,1.2471),lower=c(0,0),upper=c(5,5),verbose=3,structures=strucs,strucpars=strucpars) 
p11
summary(p11)
plot(p11,"confplot")
plot(p11,"Shepard")
#rpowerstress
p12 <- stops(dis,loss="rpowerstress",type="ratio",theta=c(1,1),lower=c(0,0),upper=c(5,5),verbose=3,structures=strucs,strucpars=strucpars) 
p12
summary(p12)
plot(p12,"confplot")
plot(p12,"Shepard")
#apstress #Issues Check again or change it to ups/tau? 
p13 <- stops(dis,loss="apstress",type="ratio",theta=c(1,1,1),lower=c(1,1,1),upper=c(3,3,3),verbose=3,structures=strucs,strucpars=strucpars) 
p13
summary(p13)
plot(p13,"confplot")
plot(p13,"Shepard")
#powerstress
p14 <- stops(dis,loss="powerstress",type="ratio",theta=c(1,1,1),lower=0,upper=5,verbose=3,structures=strucs,strucpars=strucpars) 
p14
summary(p14)
plot(p14,"confplot")
plot(p14,"Shepard")
#bcmds
p15 <- stops(dis,loss="bcstress",type="ratio",theta=c(1,1,1),lower=0,upper=5,verbose=3,structures=strucs,strucpars=strucpars) 
p15
summary(p15)
plot(p15,"confplot")
plot(p15,"Shepard")
#lmds
p16 <- stops(dis,loss="lmds",type="ratio",theta=c(3,0),upper=c(5,1),lower=c(2,0),verbose=2,structures=strucs,strucpars=strucpars) 
p16
summary(p16)
plot(p16,"confplot")
plot(p16,"Shepard")

save(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,file="oldps.rda")




                                        #try all 2nd kinship 3rd matchphi
losses <- c("stress",#1 ok ok
            "smacofSym",#1 ok ok
            "strain",#1 ok ok
            "rstress", #1 ok ok
          #  "elastic",#1 ok ok
            "apstress", #3  ok
losses <- c(           "sammon2",#1 ok ok
      "powerstrain",#1 ok ok
      "sammon", #1 ok ok
      "powerelastic",#2 ok ok wrong
      "powerstress",#3 ok ok
      "sstress", #1 ok ok
      #)  
#works with 3 works with 1 but sometimes issues with 2
#why only one value by the wrong losses?  
#losses <-     c(
    "powermds", #2 wrong, worked, wrong, worked ok     
      "powersammon",#2 wrong, worked, worked ok
"rpowerstress" #2 wrong, worked ok
  )
      
 
  

dis <- as.matrix(dis)
for(i in losses)
{
 cat("Testing:",i,"\n")   
 p1 <- pcops(dis,loss=i,type="ratio",lower=c(0.1,0.1,0.1),upper=c(5,5,5),verbose=2,itmaxi=5000)
 print(p1)
 summary(p1)
 par(mfrow=c(2,1))
 plot(p1,"confplot",label.conf=list(pos=5))
 plot(p1,"reachplot") 
# plot(p1,"Shepard")
 #plot(p1,"transplot") transplot for stress
 par(mfrow=c(1,1))
}



################################
     data(stops::Swissroll)
     dissm<-as.matrix(dist(Swissroll[,1:3]))
     cols<-Swissroll[,4]
     

     #Setting up structurenedness parameters
     strucpars<-list(list(epsilon=10,minpts=2,scale=3),list(NULL))
     
     #STOPS with strain
     resstrain<-stops(dissm,loss="strain",theta=1,structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,itmax=10)
     resstrain
     #summary(resstrain)
     plot(resstrain,col=cols,label.conf=list(col=cols))
     #Fun fact: With strain clinearity must be 0 as the
     #two principal axes are orthogonal
     #and this can't be changed by taking
     #the dissimilarities to a power
     
     
     
     #STOPS with stress or smacofSym
     resstress<-stops(dissm,loss="stress",theta=1,
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
     resstress
     #summary(resstress)
     plot(resstress,col=cols,label.conf=list(col=cols))
     plot(resstress,"Shepard")
     
     #STOPS with smacofSphere
     ressph<-stops(dissm,loss="smacofSphere",theta=1,
     structures=c("cclusteredness","cdependence"),itmaxps=5000,
     strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,verbose=3) #takes long
     ressph
     #summary(ressph)
     plot(ressph,col=cols,label.conf=list(col=cols))
     
     #STOPS with sammon
     ressam<-stops(dissm,loss="sammon",theta=1,
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
     ressam
  #   summary(ressam)
     plot(ressam,col=cols,label.conf=list(col=cols))
     plot(ressam,"transplot")
     
     #STOPS with sammon2
     ressam<-stops(dissm,loss="sammon2",theta=1,
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,verbose=2)
     ressam
   #  summary(ressam)
     plot(ressam,col=cols,label.conf=list(col=cols))
     plot(ressam,"Shepard")
     
     #STOPS with elastic 
     ressam<-stops(dissm,loss="elastic",theta=1,
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,verbose=3)
     ressam
    # summary(ressam)
     plot(ressam,col=cols,label.conf=list(col=cols))
     plot(ressam,"transplot")
     
     #STOPS with sstress #Check
     resss<-stops(dissm,loss="sstress",type="ratio",theta=0.1,
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=2,verbose=3)
     resss
     #summary(resss)
     plot(resss,col=cols,label.conf=list(col=cols))
     plot(resss,"Shepard")
     plot(resss,"transplot")
     
     #STOPS with powerstress
     respstress<-stops(dissm,loss="powerstress",
     structures=c("cclusteredness","cdependence"),theta=c(1,1,1),
     strucpars=strucpars,weightmat=dissm,
     itmaxps=5000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5),verbose=3)
     respstress
     #summary(respstress)
     plot(respstress,col=cols,label.conf=list(col=cols))
     
     #STOPS with restricted powerstress
     respstressr<-stops(dissm,loss="rpowerstress",theta=c(1,1),
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,weightmat=dissm,
     itmaxps=5000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5),verbose=3)
     respstressr
     #summary(respstressr)
     plot(respstressr,col=cols,label.conf=list(col=cols))
     
     #STOPS with powermds
     respmds<-stops(dissm,loss="powermds",
     structures=c("cclusteredness","cdependence"),theta=c(1,1),
     strucpars=strucpars,weightmat=dissm,
     itmaxps=5000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5),verbose=3)
     respmds
     #summary(respmds)
     plot(respmds,col=cols,label.conf=list(col=cols))
     
     #STOPS with powersammon
     respmds<-stops(dissm,loss="powersammon",theta=c(1,1),
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,weightmat=dissm,
     itmaxps=5000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5),verbose=3)
     respmds
     #summary(respmds)
     plot(respmds,col=cols,label.conf=list(col=cols))
     
     #STOPS with powerelastic
     respmds<-stops(dissm,loss="powerelastic",theta=c(1,1),
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,weightmat=dissm,
     itmaxps=5000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5),verbose=3)
     respmds
     #summary(respmds)
     plot(respmds,col=cols,label.conf=list(col=cols))
     
     #STOPS with ordinal rstress 
     resr<-stops(dissm,loss="rstress",type="ordinal",theta=1,
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,weightmat=dissm,
     itmaxps=5000,optimmethod="ALJ",lower=0.5,upper=5,verbose=3)
     resr
     #summary(resr)
plot(resr,col=cols,label.conf=list(col=cols))
plot(resr,"Shepard")
     
     #STOPS with approximated powerstress
     respstressa<-stops(dissm,loss="powerstress",
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,weightmat=dissm,theta=c(1,1,1),
     itmaxps=5000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5),verbose=3)
     respstressa
     #summary(respstressa)
     plot(respstressa,col=cols,label.conf=list(col=cols))
     plot(respstressa,"transplot")


     #STOPS with bcmds
     resbcstress<-stops(dissm,loss="bcmds",
     structures=c("cclusteredness","cdependence"),
     strucpars=strucpars,optimmethod="ALJ",lower=c(0.5,1,0.5),upper=c(5,5,5),verbose=2)
     resbcstress
     #summary(resbcstress)
     plot(resbcstress,col=cols,label.conf=list(col=cols))
     
     #STOPS with lmds
     reslmds<-stops(dissm,loss="lmds",theta=c(1,1),
     structures=c("cclusteredness","clinearity"),
     strucpars=strucpars,optimmethod="ALJ",lower=c(2,0),upper=c(10,1),verbose=2)
     reslmds
     #summary(reslmds)
     plot(reslmds,col=cols,label.conf=list(col=cols))
     
     #STOPS with Isomap (the k version)
     resiso2<-stops(dissm,loss="isomap",theta=5,
     structures=c("cclusteredness","clinearity"),
     strucpars=strucpars,optimmethod="ALJ",lower=3,upper=10)
     resiso2
    # summary(resiso2)
     plot(resiso2,col=cols,label.conf=list(col=cols))
     
     #STOPS with Isomap (the eps version)
     resiso<-stops(dissm,loss="isomapeps",
     structures=c("cclusteredness","clinearity"),theta=0.5,
     strucpars=strucpars,optimmethod="ALJ",lower=0.1,upper=1)
     resiso
    # summary(resiso)
     plot(resiso,col=cols,label.conf=list(col=cols))
     #################
     
     data(kinshipdelta,package="smacof")
     strucpar<-list(NULL,NULL) #parameters for indices
     res1<-stops(kinshipdelta,loss="stress",
     structures=c("cclumpiness","cassociation"),strucpars=strucpar,
     lower=0,upper=10,itmax=10)
     res1
     
     
     data(BankingCrisesDistances)
     strucpar<-list(c(epsilon=10,minpts=2),NULL) #parameters for indices
     res1<-stops(BankingCrisesDistances[,1:69],loss="stress",verbose=0,
     structures=c("cclusteredness","clinearity"),strucpars=strucpar,
     lower=0,upper=10)
     res1
     
     strucpar<-list(list(alpha=0.6,C=15,var.thr=1e-5,zeta=NULL),
     list(alpha=0.6,C=15,var.thr=1e-5,zeta=NULL))
     res1<-stops(BankingCrisesDistances[,1:69],loss="stress",verbose=0,
     structures=c("cfunctionality","ccomplexity"),strucpars=strucpar,
     lower=0,upper=10)
     res1
     


########### Test methods for jackmds 
 data(kinshipdelta)
 dissm<-as.matrix(kinshipdelta)
 strucpar<-list(list(alpha=0.6,C=15,var.thr=1e-5,zeta=NULL),
 list(alpha=0.6,C=15,var.thr=1e-5,zeta=NULL))
##smacofB object in $fit works
res1<-stops(dissm,loss="stress",verbose=2,structures=c("cfunctionality","ccomplexity"),strucpars=strucpar,lower=0,upper=0.5)
res1$fit
jk1 <- jackmds(res1$fit)
jk2 <- jackmds(res1)
jk1
jk2


#smacofP object in $fit 
res2<-stops(dissm,loss="powerstress",theta=c(1,1,1),verbose=3,structures=c("cfunctionality","ccomplexity"),strucpars=strucpar,lower=c(0,0,0),upper=c(10,10,10))
res<-stops(dissm,loss="rstress",type="ordinal",theta=1,verbose=2,structures=c("cfunctionality","ccomplexity"),strucpars=strucpar,lower=0,upper=10)
res2$fit #we see itmaxi is not substituted in the call, so it won't work.
## we may fix it by either changing it in the stops call with substitute, I believe (but super hard as we have so many of these calls), or via hard code in jackmds
res3 <- powerstressMin(dissm,itmax=100000)
res3$call
res3$call$itmax

jk1 <- jackmds(res2$fit)
jk2 <- jackmds(res)
jk1
jk2
plot(jk2)

system.time(jk2 <- jackmds(res2))
system.time(jk3 <- jackmds(res3))

r <- 0.4
fit <- rStressMin(dissm,type="ordinal",r=r,itmax=1000) ## 2D ordinal MDS
jackmds(fit)
#works, but if we rm(r) then it no longer does
fit$pars
#works

########

dats <- na.omit(PVQ40[,1:5])
dissm <- as.matrix(dist(t(dats)))   ## Euclidean distances 
 
strucpar<-list(list(alpha=0.6,C=15,var.thr=1e-5,zeta=NULL),
list(alpha=0.6,C=15,var.thr=1e-5,zeta=NULL))
#smacofB object in $fit works
res1<-stops(dissm,loss="stress",verbose=2,structures=c("cfunctionality","ccomplexity"),strucpars=strucpar,lower=0,upper=10)
res1
jk2 <- bootmds(res1,data=dats)


#smacofP in $fit works
res2<-stops(dissm,loss="powerstress",theta=c(1,1,1),verbose=3,structures=c("cfunctionality","ccomplexity"),strucpars=strucpar,lower=c(0,0,0),upper=c(10,10,10))
bs2 <- bootmds(res2,data=dats) #works

bs2 <- bootmds(res2$fit,data=dats) #doesnt work; would need to change bootmds.smacofP

res<-stops(dissm,loss="rstress",type="ratio",theta=1,verbose=2,structures=c("cfunctionality","ccomplexity"),strucpars=strucpar,lower=0,upper=10)

bs3 <- bootmds(res,data=dats) #works


#t1 <- stop_rstress(dissm,theta=1,structures="clinearity")
#t1$fit
#r
#rm(r)
#eval(t1$fit$call)

#t2 <- stop_rstress2(dissm,theta=1,structures="clinearity")
#t2$fit
#itmaxi <- 10000
#eval(t2$fit$call)
