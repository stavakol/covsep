library(tools)
library(parallel)
cores <- detectCores()
print(cores)

library(covsep)
N  <- 100

## testing 'rmtnorm'
Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)))
str(Data)
Data.meanzero = rmtnorm(N, C1, C2)
str(Data)

## testing 'renormalize_mtnorm'
ans  <-  renormalize_mtnorm(Data.meanzero[1,,], C1, C2, type='diag')
str(ans)
ans  <-  renormalize_mtnorm(Data.meanzero[1,,], C1, C2, type='full')
str(ans)
ans=sapply(1:N, function(x){ sum( ( renormalize_mtnorm(Data.meanzero[x,,], C1,C2, type='full')[1:2,1:2])^2) })
qqplot(qchisq(ppoints(N), df=4), ans)
abline(0,1)

## testing 'marginal_covariances'
marginal.covariances  <-  marginal_covariances(Data)
str(marginal.covariances)
print( sum((C1 - marginal.covariances$C1)^2) )
print( sum((C2 - marginal.covariances$C2)^2) )

## testing 'projected_differences'
ans <- projected_differences(Data, l1=1, l2=1)
str(ans)
ans <- projected_differences(Data, l1=3, l2=2)
str(ans)
ans <- projected_differences(Data, l1=1, l2=1, with.asymptotic.variances=FALSE)
str(ans)
ans <- projected_differences(Data, l1=3, l2=2, with.asymptotic.variances=FALSE)
str(ans)

## testing 'clt_test'
# data with non-zero mean
print(clt_test(Data, L1=1, L2=1) )
print(clt_test(Data, L1=1, L2=3) )
print(clt_test(Data, L1=3, L2=1) )
print(clt_test(Data, L1=c(1,2), L2=c(1,4)) )
nrep  <- 100
pvals <- replicate(nrep, {
          Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
          clt_test(Data, L1=1:3, L2=1:3)
})
op <- par(mfrow=c(3,1), oma=c(0,0,3,0))
{
hist(pvals[1,])
hist(pvals[2,])
hist(pvals[3,])
    title(main='Pvalues of CLT test (Gaussian data with non-zero mean)', outer=TRUE)
}

# data with mean zero
print(clt_test(Data.meanzero, L1=1, L2=1) )
print(clt_test(Data.meanzero, L1=1, L2=3) )
print(clt_test(Data.meanzero, L1=3, L2=1) )
print(clt_test(Data.meanzero, L1=c(1,2), L2=c(1,4)) )
nrep  <- 100
pvals.meanzero <- replicate(nrep, {
          Data.meanzero = rmtnorm(N, C1, C2) 
          clt_test(Data.meanzero, L1=1:3, L2=1:3)
})
op <- par(mfrow=c(3,1), oma=c(0,0,3,0))
{
hist(pvals.meanzero[1,], xlim=0:1)
hist(pvals.meanzero[2,], xlim=0:1)
hist(pvals.meanzero[3,], xlim=0:1)
    title(main='Pvalues of CLT test (Gaussian data with mean zero)', outer=TRUE)
}



### testing 'gaussian_bootstrap_test'
nrep  <- 100
B  <- 1000


system.time( gaussian_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='no', B=B) )

system.time({
    replicate(4, gaussian_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='no', B=100, parallel=FALSE) ) 
})
gaussian_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='diag', B=B)
gaussian_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='full', B=B)

## verif distribution of p-values
# setup cluster
#cluster <- makePSOCKcluster(cores)
cluster <- makePSOCKcluster(2)
clusterExport(cluster, c(ls(), 'C1', 'C2', 'psnice', ls('package:covsep')))

nrep=10
B=100
print('--time without parallel--')
system.time(
            {
            replicate(nrep, 
                      {
                 Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
                gaussian_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='no', B=B);
                      }
            )
}
            )
print('--time with parallel--')
system.time({
    ans=simplify2array(mclapply(1:nrep, function(j){ 
                 Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
                 return(gaussian_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='no', B=B))
}, mc.cores = 3 , mc.silent = TRUE ) )
})

system.time(parLapply(cluster, 1:4, function(x) {
                     psnice(value = 19)
                                      Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
                     gaussian_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='no', B=B)
}
                     ))

ptm=proc.time()
pvals.nostud=simplify2array(parLapply(cluster, 1:nrep, function(x) {
                     psnice(value = 19)
                                      Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
                                      gaussian_bootstrap_test(Data,
                                                              L1=1:3,
                                                              L2=1:3,
                                                              studentize='no',
                                                              B=B)
}
)
)
tot.time= proc.time()-ptm
print(tot.time)


ptm=proc.time()
pvals.diag=simplify2array(parLapply(cluster, 1:nrep, function(x) {
                     psnice(value = 19)
                                      Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
                                      gaussian_bootstrap_test(Data,
                                                              L1=1:3,
                                                              L2=1:3,
                                                              studentize='diag',
                                                              B=B)
}
)
)
tot.time= proc.time()-ptm
print(tot.time)

ptm=proc.time()
pvals.full=simplify2array(parLapply(cluster, 1:nrep, function(x) {
                     psnice(value = 19)
                                      Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
                                      gaussian_bootstrap_test(Data,
                                                              L1=1:3,
                                                              L2=1:3,
                                                              studentize='full',
                                                              B=B)
}
)
)
tot.time= proc.time()-ptm
print(tot.time)

#shut down cluster
stopCluster(cluster)

# plot results
op <- par(mfrow=c(3,1), oma=c(0,0,3,0))
{
    for(i in 1:3){
        hist(pvals.nostud[i,], xlim=c(0,1))
    }
    title(main='Gaussian Bootstrap, no stud', outer=TRUE)
    for(i in 1:3){
        hist(pvals.diag[i,], xlim=c(0,1))
    }
    title(main='Gaussian Bootstrap, diag stud', outer=TRUE)
    for(i in 1:3){
        hist(pvals.full[i,], xlim=c(0,1))
    }
    title(main='Gaussian Bootstrap, full stud', outer=TRUE)
}
par(op)

### testing 'empirical_bootstrap_test'
nrep  <- 100
B  <- 100

empirical_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='no', B=B)
empirical_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='diag', B=B)
empirical_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='full', B=B)

## verif distribution of p-values
# setup cluster
cluster <- makePSOCKcluster(cores)
clusterExport(cluster, c(ls(), 'C1', 'C2', 'psnice'))
system.time(parLapply(cluster, 1:4, function(x) {
                     psnice(value = 19)
                                      Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
                     empirical_bootstrap_test(Data, L1=1:3, L2=1:3, studentize='no', B=B)
}
                     ))

ptm=proc.time()
pvals.nostud=simplify2array(parLapply(cluster, 1:nrep, function(x) {
                     psnice(value = 19)
                                      Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
                                      empirical_bootstrap_test(Data,
                                                              L1=1:3,
                                                              L2=1:3,
                                                              studentize='no',
                                                              B=B)
}
)
)
tot.time= proc.time()-ptm
print(tot.time)


ptm=proc.time()
pvals.diag=simplify2array(parLapply(cluster, 1:nrep, function(x) {
                     psnice(value = 19)
                                      Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
                                      empirical_bootstrap_test(Data,
                                                              L1=1:3,
                                                              L2=1:3,
                                                              studentize='diag',
                                                              B=B)
}
)
)
tot.time= proc.time()-ptm
print(tot.time)

ptm=proc.time()
pvals.full=simplify2array(parLapply(cluster, 1:nrep, function(x) {
                     psnice(value = 19)
                                      Data = rmtnorm(N, C1, C2, matrix(2, nrow(C1), nrow(C2)) )
                                      empirical_bootstrap_test(Data,
                                                              L1=1:3,
                                                              L2=1:3,
                                                              studentize='full',
                                                              B=B)
}
)
)
tot.time= proc.time()-ptm
print(tot.time)

#shut down cluster
stopCluster(cluster)


# plot results
op <- par(mfrow=c(3,1), oma=c(0,0,3,0))
{
    for(i in 1:3){
        hist(pvals.nostud[i,], xlim=c(0,1))
    }
    title(main='Empirical Bootstrap, no stud', outer=TRUE)
    for(i in 1:3){
        hist(pvals.diag[i,], xlim=c(0,1))
    }
    title(main='Empirical Bootstrap, diag stud', outer=TRUE)
    for(i in 1:3){
        hist(pvals.full[i,], xlim=c(0,1))
    }
    title(main='Empirical Bootstrap, full stud', outer=TRUE)
}
par(op)


## testing 'difference_fullcov'

system.time( ans <- difference_fullcov(Data))
str(ans)

## testing 'HS_gaussian_bootstrap_test'

system.time(  ans <- HS_gaussian_bootstrap_test(Data, B=100))
str(ans)

## testing 'HS_empirical_bootstrap_test'

system.time(  ans <- HS_empirical_bootstrap_test(Data, B=100))
str(ans)

## testing 'generate_surface_data'

set.seed(1)
ans1 <- rmtnorm(N, C1, C2) 
set.seed(1)
ans2 <- generate_surface_data(N, C1, C2, gamma=0, distribution='gaussian') 

system.time(replicate(3, rmtnorm(N, C1, C2) ))
system.time(replicate(3, generate_surface_data(N, C1, C2, gamma=0, distribution='gaussian') ))


######### automatic simulations

system.time({
tmp=lapply(simulation.settings[1:200], function(args){
                do.call(simulation_covsep, c(list(outdir="/home/st624/tmp/"), args))
})
})


ans=simulation_covsep(outdir='.', sim.name, simulation.seed, N, gam, B, nrep,
                  studentize, distribution, pval.method, alpha, l1, l2)

##

l=list()
for (i in 1:10){
    l[[i]] <- list(a=i, b=0)
}

