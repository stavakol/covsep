### parameters

## # sim.name<-"test-gaussian-boot"        # name of the simulation
## # #simulation.seed=2       # simulation seed
## # N<-20                    # number of observations
## # gam<-0.1      # violation of the null hypothesis [reasonable values between 0 and 0.1]
## # B<-100                    # bootstrap replicates used in the procedures
## # nrep<- 10               # number of simulations replication for the power analysis
## # studentize=TRUE           # boolean wich control if the studentized version of the test is used
## # sim.norm=TRUE             # boolean. If TRUE, simulated data come from a Gaussian process, otherwise they are generated from a multivariate t distribution on the discretized grid
## # pval.method="HS.GAUSSBOOT"   ## or "GAUSS.BOOT" or "EMP.BOOT" or "HS.EMPBOOT", "HS.GAUSSBOOT", "CLT"
## # ##  GAUSS.BOOT=FALSE         # boolean. If TRUE, the test based on the parametric gaussian bootstrap is used.
## # ##  EMP.BOOT=TRUE            # boolean. If TRUE, the test based on the empirical bootstrap is used.
## # alpha<-0.05               # level of the test
## # #l1=c(1,2,3,3,4,5,5,32)                #  basis in the first direction
## # #l2=c(1,2,2,3,3,3,4,7)                     #  basis in the second direction
## # #l1=c(1,2,10)                
## # #l2=c(1,2,4)           # eigendirection to be used in the simulation ## DO NOT LEAVE ANY SPACES
## # l1=1
## # l2=1



#####################################################################
#                                                                   #  
#  Script for the power simulation study of separability tests      #
#                                                                   #  
#####################################################################

# Get info for the simulation
system("hostname")
date()
sessionInfo()

library(mvtnorm)


## DEBUG
#options(error=recover)

# get arguments from command line
args <- commandArgs(TRUE)
print(ls())
print(args)
for(arg in args) eval(parse(text = arg))
rm(arg, args)

## load list of variables
source("source/variable-list.R")
##check if all variables are specified
for(i in 1:length(var.names)) if (!exists(var.names[i])){ stop(paste("Variable", var.names[i], "is not defined, but is needed for the script")) }


# import functions for computation of the test statistic and bootstrap tests.

#source("source/load-functions.R")

# this need to be substituted in the load-functions.R. The name of the
# functions for the test statistics and the test are not changed, but the output is different.
source("source/bootstrap_tests_rectangles.R")
source("source/test-statistic.R")
source("source/HS_test.R")
#source("source/test-statistic-empboot-rectangles.R")


## check that the parameters are not conflicting
if( !(pval.method %in% c("CLT", "GAUSS.BOOT", "EMP.BOOT", "HS.EMPBOOT", "HS.GAUSSBOOT")) ){
  stop(paste("pvalue method '", pval.method, "' undefined"))
}

## Upload the covariance function used in the simulations

load("source/Cov_for_simulations.RData")

# the vector for recording the pvalues
#pvalue<-NULL
if(pval.method %in% c("GAUSS.BOOT","EMP.BOOT")){
    pvalue=array(NA, c(nrep, length(l1)))
} else if (pval.method %in% c("CLT","HS.EMPBOOT", "HS.GAUSSBOOT")){
    pvalue=array(NA, c(nrep))
}



# Vectorization map to project 4d cov to 2d

I<-matrix(0,dim(C1)[1]*dim(C2)[1],2)

for (i in 1:dim(I)[1]){
  I[i,]<-c(ceiling(i/dim(C2)[1]),i-(ceiling(i/dim(C2)[1])-1)*dim(C2)[1])
}



# Define the covariance operator for the simulation, based on gam

CL<-matrix(0,dim(C1)[1]*dim(C2)[1],dim(C1)[1]*dim(C2)[1])
x<-seq(0,1,len=dim(C1)[1])
y<-seq(0,1,len=dim(C2)[1])

for (p in 1:(dim(C1)[1]*dim(C2)[1])){
  for (q in 1:(dim(C1)[1]*dim(C2)[1])){
    CL[p,q]<-gam*(1/((I[p,2]-I[q,2])^2+1)*exp(-(I[p,1]-I[q,1])^2/((I[p,2]-I[q,2])^2+1)))+(1-gam)*C1[I[p,1],I[q,1]]*C2[I[p,2],I[q,2]]
  }
}


# correlation operator if needed for the multivariate t simulations

if(!sim.norm){   
  CLC<-cov2cor(CL)
}



# set simulation seed
set.seed(simulation.seed)

# variable to measure running time 
ptm=proc.time()

# simulation replicates

for (it in 1:nrep){ 
  
  ## print simulation progress
  if((it %% floor(1+nrep/200))==0){ cat("\n**********  Rep", round(100*it/nrep), "%  ***********\n", sep="") }
  
  # data simulation
  
  if(sim.norm){ 
    X<-rmvnorm(N,mean=rep(0,dim(CL)[1]),sigma=CL)
  }else{   
    X<-rmvt(N,sigma=CLC,df=6)
  }
  
  # inverse of the vectorization:
  
  Y<-array(0,c(N,dim(C1)[1],dim(C2)[1]))
  
  for (i in 1:dim(X)[2]){
    Y[,ceiling(i/dim(Y)[3]),i-(ceiling(i/dim(Y)[3])-1)*dim(Y)[3]]<-X[,i]
  }
  
  
  # pvalue 
  
  if(pval.method == "GAUSS.BOOT"){
    #pvalue<-rbind(pvalue,gauss.boot.test(Y,l1,l2,studentize=studentize,B))
    pvalue[it,]=gaussian_bootstrap_test(Y,l1,l2,studentize=studentize,B)
  } else if(pval.method == "EMP.BOOT"){
    pvalue[it,]=empirical_bootstrap_test(Y,l1,l2,studentize=studentize,B)
    #pvalue<-rbind(pvalue,emp.boot.test(Y,l1,l2,studentize=studentize,B))
  } else if (pval.method == "HS.EMPBOOT"){
      pvalue[it]=HS_empirical_bootstrap_test(Y, B)
  } else if (pval.method == "HS.GAUSSBOOT"){
      pvalue[it]=HS_gaussian_bootstrap_test(Y, B)
  } else if (pval.method == "CLT"){
      pvalue[it] = clt_test(Y,l1,l2)
  }
}

# compute total running time
tot.time=proc.time() - ptm
print(tot.time)

print(paste("finished at ", date() ) )

# power calculation
  if(pval.method %in% c("GAUSS.BOOT","EMP.BOOT")){
P=apply(pvalue,2, function(x){return( sum(x <= alpha)/length(x))})
  } else if (pval.method %in% c("HS.EMPBOOT", "HS.GAUSSBOOT", "CLT")){
P= sum(pvalue <= alpha)/length(pvalue)
  }

# name of variables to save
to.save=c(var.names, 
          "P",
          "pvalue",
          "N",
          "gam",
          "B",
          "nrep",
          "alpha",
          "studentize",
          "sim.norm",
          "pval.method",
          "l1",
          "l2",
          "tot.time"
)

# create a list with the variables to save
results=list()
for(i in 1:length(to.save)){
  eval(parse(text=paste("results$", to.save[i], "=", to.save[i], sep="")))
}


#save(P,pvalue,N,gam,B,nrep,alpha,studentize,sim.norm,CLT,GAUSS.BOOT,EMP.BOOT,l1,l2,file=paste("results/results_sim-", sim.name, "-",Sys.time(),".RData",sep=""))
save(results, file=paste("results/", sim.name, ".RData",sep=""))



##############
