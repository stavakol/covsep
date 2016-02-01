simulation_covsep  <- function(outdir='.', sim.name,
                               simulation.seed,
                               N,
                               gam,
                               B,
                               nrep,
                               studentize,
                               distribution,
                               pval.method,
                               alpha,
                               l1,
                               l2, return.results=FALSE){

    #sink(file=paste0(outdir, sim.name, '.Rout'))

    hostname=system("hostname", intern=TRUE)
    cat(paste(hostname, Sys.time(), Sys.getpid(), toString(sessionInfo()), sep='\n'),
        file=paste0(outdir, sim.name, '.Rout'), append=TRUE)

    ## load package
    library(covsep)

    # the vector for recording the pvalues
    #pvalue<-NULL
    if(pval.method %in% c("CLT","GAUSS.BOOT","EMP.BOOT")){
        pvalue=array(NA, c(nrep, length(l1)))
    } else if (pval.method %in% c("HS.EMPBOOT", "HS.GAUSSBOOT")){
        pvalue=array(NA, c(nrep))
    }

    # set simulation seed
    set.seed(simulation.seed)

    # variable to measure running time 
    ptm=proc.time()

    # simulation replicates

    for (it in 1:nrep){ 
        save(.Random.seed, file=paste0(outdir, sim.name, '.RData.dump'))

        ## print simulation progress and hearbeat
        #if((it %% floor(1+nrep/200))==0){ 
            cat(paste(Sys.time(), Sys.info()['nodename'], 'PID:', Sys.getpid()), '\n')
            cat(paste(Sys.time(), Sys.info()['nodename'], 'PID:', Sys.getpid()), '\n', file=paste0(outdir, sim.name, '.Rout'), append=TRUE)  
            cat("\n**********  Rep", round(100*it/nrep, 2), "%  ***********\n", sep="", file=paste0(outdir, sim.name, '.Rout'), append=TRUE) 
    #}

        # data simulation
        Y <- generate_surface_data(N, C1, C2, gamma=gam, distribution=distribution) 


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
            pvalue[it, ] = clt_test(Y,l1,l2)
        }
    }; rm(it)

    # compute total running time
    tot.time=proc.time() - ptm
    print(tot.time)

    print(paste(hostname, " finished at ", date() ) )

    # power calculation
    if(pval.method %in% c("GAUSS.BOOT","EMP.BOOT", "CLT")){
        P=apply(pvalue,2, function(x){return( sum(x <= alpha)/length(x))})
    } else if (pval.method %in% c("HS.EMPBOOT", "HS.GAUSSBOOT")){
        P= sum(pvalue <= alpha)/length(pvalue)
    }

    # name of variables to save
    to.save=c(c("sim.name",
              "simulation.seed",
              "N",
              "gam",
              "B",
              "nrep",
              "studentize",
              "pval.method",
              "alpha",
              "l1",
              "l2"), 
              "P",
              "pvalue",
              "N",
              "gam",
              "B",
              "nrep",
              "alpha",
              "studentize",
              "pval.method",
              "l1",
              "l2",
              "tot.time"
              )

#    # create a list with the variables to save
#    results=list()
#    for(i in 1:length(to.save)){
    #results[[to.save[i]]] = 
#        eval(parse(text=paste("results$", to.save[i], "=", to.save[i], sep="")))
#    }
    #results=sapply(ls(), get)
    results=mget(ls())


    #save(P,pvalue,N,gam,B,nrep,alpha,studentize,CLT,GAUSS.BOOT,EMP.BOOT,l1,l2,file=paste("results/results_sim-", sim.name, "-",Sys.time(),".RData",sep=""))
    save(results, file=paste(outdir, sim.name, ".RData",sep=""))
    if(return.results)    return(invisible(results))
    else return(0)

}

##############
