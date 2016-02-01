library(covsep)

arguments.names <- c("sim.name", "simulation.seed", "N", "gam", "B",
                     "nrep", "studentize", "distribution", "pval.method", "alpha",
                     "l1", "l2");
##############################

#outdir="/home/st624/tmp/"
outdir="/home/st624/simulations/JA-DP-ST/cov-separability/results/"
prepend="SIM2016a"  ## same as 'SIM19jan', but with projection on "rectangles", and HS norm tests
first.seed=2000 # so it doesn't interfere with SIM12feb

## Parameters that change
N.seq<-c(10,25,50,100)                     # number of observations
gam.seq<-seq(0,.1,by=.01)                    # violation of the null hypothesis [reasonable values between 0 and 0.1]
pval.method.seq=c("CLT","EMP.BOOT","GAUSS.BOOT")
#pval.method.seq=c("GAUSS.BOOT")
sim.norm.seq=c(TRUE, FALSE)             # boolean. If TRUE, simulated data come from a Gaussian process, otherwise they are generated from a multivariate t distribution on the discretized grid
studentize.seq=c('full')           # boolean which control if the studentized version of the test is used

## Parameters that are fixed
nrep=1000
B=1e3
l1=c(1,2,10)                
l2=c(1,2,4)           # eigendirection to be used in the simulation ## DO NOT LEAVE ANY SPACES
alpha<-0.05               # level of the test
##########################################################################################

## Length of each changing parameter
N.seq.len=          length(N.seq)          
gam.seq.len=        length(gam.seq)
pval.method.seq.len=length(pval.method.seq)
sim.norm.seq.len=   length(sim.norm.seq)
studentize.seq.len =length(studentize.seq)

# dimensions of paremeters
param.dims=c(N.seq.len,          
             gam.seq.len,        
             pval.method.seq.len,
             sim.norm.seq.len,
             studentize.seq.len)

# seeds for all the simulations
seeds=array(seq(from=first.seed, by=1, len=prod(param.dims)), dim=param.dims)

## ans=foreach(noise.j = 1:noise.seq.len) %dopar% {
simulation.settings = list()
j=0 #setting number

for(N.i in 1:N.seq.len){
    N=N.seq[N.i]
    for(gam.i in 1:gam.seq.len){
        gam=gam.seq[gam.i]
        for(pval.method.i in 1:pval.method.seq.len){
            pval.method=pval.method.seq[pval.method.i]
            for(sim.norm.i in 1:sim.norm.seq.len){
                sim.norm=sim.norm.seq[sim.norm.i]
                for(studentize.i in 1:studentize.seq.len){
                    studentize=studentize.seq[studentize.i]

                    # set seed for the simulation
                    simulation.seed=seeds[N.i,gam.i,pval.method.i,sim.norm.i, studentize.i]
                    distribution  <- ifelse(sim.norm, 'gaussian', 'student')

                    sim.name=paste(prepend, "-N", N, "-gam", gam, "-pval.method", pval.method, "-sim.norm", sim.norm, "-stud", studentize,  sep="") 

                    ## only generate simulation setting if .RData file doesn't exist
                    if(!file.exists( paste(outdir, sim.name, ".RData",sep="")) && !file.exists( paste(outdir, sim.name, ".Rout",sep="")) ){
                       j  <- j+1; # update setting number
                       simulation.settings[[j]] = sapply(arguments.names, get)
             }
                }
            }
        }
    }
}
rm(simulation.seed)

print(warnings())

##############################
## wake "crow" "sloth" "garnish" "otter" "buffoon" "nordic" "glencoe" "yogi" "platypus" "solstice"

#juno 
COMPUTERS="buffoon cave equinox garnish glencoe goldrush grumpy haystack hiccup nordic norinori nosey otter oval pentopia platypus primrose rainbow saratoga senate sloth snake solstice stadium suraromu tomtom tapa"

FASTCOMPUTERS="gargoyle orwell snake garnish chiffchaff heraclius jigsaw norinori eris numberlink pluto janus minesweeper chopin haydn shikaku himalia metis anice statuepark rossini mozart socrates elara thebe sinope brahms paaliaq persaeus mustard magnets wagner scramble beans wryneck basil sun jupiter rubik democrates suguru hashi porridge tangram siarnaq verdi epicurus cayenne"

FASTCOMPUTERS.tmp="gargoyle orwell chiffchaff heraclius jigsaw eris numberlink pluto janus minesweeper chopin haydn shikaku himalia metis anice statuepark rossini mozart socrates elara thebe sinope brahms paaliaq persaeus mustard magnets wagner scramble beans wryneck basil sun jupiter rubik democrates suguru hashi porridge tangram siarnaq verdi epicurus cayenne"

ALLCOMPUTERS="equinox parrot bonxie buzzard blackcap parakeet gull kite plover whinchat skua sanderling primrose saratoga glencoe treecreeper minesweeper porridge hashi hippo democrates chopin yogi haydn vulture tacitus tapa platypus himalia metis magnets verdi grumpy jigsaw siarnaq nosey norinori hiccup otter eris heraclius numberlink thales pluto tangram janus wagner buffoon acrion mustard thebe persaeus paaliaq neverland sinope brahms rubik anice sloth statuepark rossini tapir jupiter sun basil mozart socrates fulmar epicurus elara inglorious beans scramble diogenes cayenne garnish gargoyle suguru shikaku wryneck stonechat magpie waxwing bunting snake chiffchaff senate cave nordic stadium pentopia haystack oval tomtom quail harrier orwell diver"

USED.COMP='garnish'

#system2("/alt/bin/wakeup", args=USED.COMP)
computers = unlist(strsplit(USED.COMP, ' ') )

all.cpus = c()
for(comp in computers)
{
    #system2("/usr/bin/ssh",  args=paste0(comp, " -f 'killall R'"))
    ncpu = system(paste0("/usr/bin/ssh ",  comp, " -f 'nproc'"), intern=TRUE)
    if( is.null(attr(ncpu, 'status')) ) # add computer only if it's responding
        all.cpus  <- c(all.cpus, rep(comp, 1)) # just one cpu per computer
            #c(all.cpus, rep(comp, as.integer(ncpu)-1)) # leave one cpu free
}

##cl : set up parallel cluster on 10 machines
library(tools)
library(parallel)
source("/local/scratch/public/st624/Dropbox/work/JA-DP-ST/Report Separability/new_package/covsep/simulations/sim_function.R")
#c("crow","sloth","garnish","otter","buffoon","nordic","glencoe","yogi","platypus","solstice")
cl <- makePSOCKcluster(names = all.cpus,  outfile=' /home/st624/sep-paper/new_package/covsep/simulations/cluster.out')
                       
#cluster export: send information to cluster
clusterExport(cl, c('simulation_covsep', 'psnice', 'outdir'))

system.time(
            parLapplyLB(cl, simulation.settings, function(args){
                        psnice(value = 19)
                do.call(simulation_covsep, c(list(outdir=outdir), args))
})
)
#                do.call(simulation_covsep, c(list(outdir=outdir), simulation.settings[[1]]))

stopCluster(cl)

