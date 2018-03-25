###############################################
#           SESIM - microcosm model           #
###############################################

library(vegan)
library(igraph)
library(reshape)
library(ggplot2)
library(gplots)
library(plotrix)

## Setting up metacommunity

# number of patches
numCom <- 5

### Connectivity

# absent
d.a <- matrix(c(Inf,Inf,Inf,Inf,Inf), nrow=5, ncol=5)
diag(d.a) <- 0
d.a
# linear
d.l <- matrix(c(0,1.0,Inf,Inf,Inf,1.0), nrow=5, ncol=5)
d.l
# circular
d.c <- d.l
d.c[5,1] <- d.c[1,5] <- 1.0
d.c
# global
d.g <- matrix(c(1.0,1.0,1.0,1.0,1.0), nrow=5, ncol=5)
diag(d.g) <- 0
d.g
connec <- c("d.a", "d.l", "d.c", "d.g")

# dispersal rate
rate <- 0.001

# species-specific dispersal ability
disp <- c(0, 0.13, 0.13, 0.13, 0.64, 0.62, 0.53, 0.90, 0.62, 0.44, 0.33, 0.71, 1.00)

# dispersal kernel
k <- 1

### Patch quality

# growth rates
r.source <- c(0,8.4 ,  16.8,   20.3, 1.32, 0.92, 1.64,  0.11, 2.63,  0.15, 1.78, 0.48, 1.31)
r.sink   <- c(0,0.084 ,  0.168,   0.203, 1.82, 0.16, 1.45,  0.00, 1.25,  0.77, 2.04, 1.30, 1.31)

# environmental effect per patch quality
e.source <- 1
e.sink <- 0.1


### Species

# food web
FW <- read.csv(file="pc-fw.csv", header=T)
# number of species
nSp <- nrow(FW)

# cell mass used to calculate productivity
mass <- c(4.3e-08,1.96e-07,3.76e-06,1.52e-08,4.77e-09,6.9e-08,9.68e-08,8.05e-08,2.27e-07)


### Environment

# number of environmental variables
nEnv <- 1
# environmental niche amplitude
eA <- 1
# strength of environmental effect
enveff <- 1
# environmental fluctuation period
eP <- 1


### Simulation parameters

# number of replicates
reps <- 10

# seeds
xseed <- as.integer(runif(reps)*100000)

# simulation length
Tmax <- 500

# sampling frequency
sampfreq <- 10
sampleV <- seq(0, Tmax, by=sampfreq)

# empty array for saving sampled data
X_save <- array(data=NA, dim=c(numCom, nSp, length(sampleV)))
light  <- array(data=NA, dim=c(numCom, 1  , length(sampleV)))

### Experimental design

# experiment factorial design
experiment <- data.frame(connectivity = rep(NA, each=length(connec)*length(ss.prop)*
                                              length(rate)*length(k)*length(eP)), 
                         patchquality = rep(NA),
                         dispersal=NA, kernel=NA, fluctuation=NA)

# data frames for results
results <- data.frame(loop = rep(1:(dim(experiment)[1]*numCom*reps)),
                      treatment = rep(1:(dim(experiment)[1]), each=5*reps),
                      r = rep(1:reps, each=5),
                      patch = rep(1:5))



## Model

treatment_id = 1
loop_id = 1

for(a in rate){
  for(dd in k){
    for(ePeriod in eP){
      for(eAMP in eA){
        for(conn in seq_along(connec)){
          for(ss in ss.prop){
            repeat{
              for(r in 1:reps){
                
                sampling = 1
                
                # setting replicate seed
                set.seed(xseed[r])
                
                # patch quality
                C <- matrix(c(rep(r.source, times=ss), rep(r.sink, times=(5-ss))), nrow=nSp, ncol=5)
                
                # dispersal ability
                cc <- get(connec[conn])
                d_exp <- exp(-dd*cc) - diag(nrow(cc))
                dispersal_m <- apply(d_exp, 1, function(x) x/sum(x))
                dispersal_m[is.na(dispersal_m)] <- 0
                
                # species interactions strenghts
                BB  <- FW
                matriz  <- BB
                novamatriz <- matriz
                for(j in 1:ncol(matriz))
                  for(i in 1:nrow(matriz))
                    novamatriz[i,j] <- runif(1, min=min(0, matriz[i,j]), max=max(0, matriz[i,j]))
                BB <- novamatriz/nSp
                
                # species initial abundances
                low  <- c(100000 ,100000, 100000 ,100000 ,20,100,20,100,100,20,20,20,100)
                high <- c(5000000,5000000,5000000,5000000,20,100,20,100,100,20,20,20,100)
                X <- matrix(c(rep(high, times=ss), rep(low, times=(5-ss))), nrow=nSp, ncol=5)
                Xd <- t(X)
                X_I <- array(data=NA, dim=c(numCom, nSp, reps))
                X_I[,,] <- Xd
                
                # initial environment
                if(ePeriod != 1){
                  envt.v <- c(1,0,1,0,1)
                } else {
                  envt.v <- c(1,1,1,1,1)
                }
                
                # species environmental optimum
                Env_Opt <- matrix(rep(0.5), nSp, nEnv)
                
                # and
                # empty data frames for food wed results
                fw.results   <- matrix(NA, nrow=numCom, ncol=4)
                fw.results.r <- rep(NA, times=4)
                
                # simulation of time steps
                for(l in 1:(Tmax)){
                  
                  # for starting conditions at 1
                  if(ePeriod != 1){
                    env.new <- rep(NA, 5)
                    for(e in seq_along(envt.v)){
                      if((e==1 | e==3 | e==5) & envt.v[e] != 0) env.new[e] <- envt.v[e] - (1/ePeriod)
                      if((e==1 | e==3 | e==5) & envt.v[e] == 0) env.new[e] <- 1 - abs(envt.v[e] - (1/ePeriod))
                      if((e==2 | e==4)        & envt.v[e] != 1) env.new[e] <- envt.v[e] + (1/ePeriod)
                      if((e==2 | e==4)        & envt.v[e] == 1) env.new[e] <- 1 - abs(envt.v[e] - (1/ePeriod))
                    }
                    envt.v <- round(env.new, digits=nchar(ePeriod))
                  }
                  
                  # environmental condition in each patch per species
                  Env <- matrix(rep(envt.v, each=nSp), nSp, numCom)
                  # environmental match
                  enveff <- c(rep(e.source, times=ss), rep(e.sink, times=(5-ss)))
                  enviro <- t(apply( abs((Env-Env_Opt[,1])), 1, function(x) x * enveff))
                  enviro[is.infinite(enviro)] <- 0
                  diff <- C/enviro
                  diff[is.na(diff)] <- enviro[is.na(diff)]
                  diff[is.infinite(diff)] <- enviro[is.infinite(diff)]
                  diff[diff==0] <- enviro[,]
                  
                  # interactions
                  interactions <- ((BB) %*% (X*(diff)))
                  interactions[is.na(interactions)] <- 0
                  interactions[X == 0] <- 0
                  
                  # dispersal
                  cMigrants <- matrix(NA, nSp, numCom)
                  for(j in 1:numCom){
                    cMigrants[,j] <- X[,j]*(1/(length(dispersal_m[j,][dispersal_m[j,] > 0 & dispersal_m[j,] != Inf])))
                  }
                  Migrants <- apply(cMigrants, 2, function(x) x * (a/disp))
                  Migrants[is.na(Migrants)] <- 0
                  Migrants[is.infinite(Migrants)] <- 0
                  if(conn == 1) Migrants[] <- 0
                  cImmigrants <- matrix(NA, nSp, numCom)
                  for(i in 1:nSp){
                    for(j in 1:numCom){
                      cImmigrants[i,j] <- sum((Migrants[i,])*dispersal_m[j,])
                    }
                  }
                  Immigrants <- cImmigrants
                  Immigrants[is.na(Immigrants)] <- 0
                  Immigrants[is.infinite(Immigrants)] <- 0
                  
                  # Lotka-Volterra model
                  Xt <- X + interactions + Immigrants - Migrants
                  
                  Xhold <- Xt
                  Xhold[(Xhold < 0.1)] <- 0
                  Xhold[is.na(Xhold)] <- 0
                  Xhold[is.infinite(Xhold)] <- 0
                  X <- Xhold
                  Xd <- t(X)
                  if(all(is.nan(X))) break
                  
                  if(l==sampleV[l/sampfreq+1] && l<=Tmax){
                    X_save    [,,sampling] <- Xd
                    light     [,,sampling] <- envt.v
                    sampling <- sampling + 1
                  } # end sampling
                  
                } # end all time steps
                
                # saving sampled data for one replicate and one treatment level
                treat <- data.frame(time=(rep(seq_len(dim(X_save)[3]),each=dim(X_save)[1]))*sampfreq, 
                                    dispersal = rep(rate), 
                                    kernel = rep(k), 
                                    replicate = rep(r), 
                                    patch = rep(1:5),
                                    resource    = c(rep("high", times=ss), rep("low", times=(5-ss))),
                                    fluctuation = rep(ePeriod),
                                    connectivity = rep(conn),
                                    quality = rep(ss))
                X_saved <- as.data.frame(apply(X_save, 2, cbind))
                assign(paste("rep_sampled_", loop_id, "_", r, sep=""), cbind(treat, X_saved))
                
                treat.I <-
                  data.frame (time        = (rep(1, each=dim(X_I)[3])),
                              dispersal   = rep(a),
                              kernel      = rep(dd),
                              replicate   = rep(r),
                              patch       = rep(1:5),
                              resource    = c(rep("high", times=ss), rep("low", times=(5-ss))),
                              fluctuation = rep(ePeriod),
                              connectivity = rep(conn),
                              quality      = rep(ss),
                              light        = rep("1"))
                X_initial <- as.data.frame(apply(X_I, 2, cbind))
                X_saved <- as.data.frame(apply(X_save, 2, cbind))
                light_saved <-  as.data.frame(apply(light, 2, cbind))
                names(light_saved)[names(light_saved)=="V1"] <- "light"
                assign(paste("rep_sampled_", loop_id, "_", r, sep=""), rbind(cbind(treat.I, X_initial),cbind(treat, light_saved, X_saved)))
                
                # adding treatment level to experiment data frame
                experiment$connectivity [treatment_id] <- connec[conn]
                experiment$patchquality [treatment_id] <- ss
                experiment$dispersal  [treatment_id] <- a
                experiment$kernel     [treatment_id] <- dd
                experiment$fluctuation[treatment_id] <- ePeriod
                
                # adding treatment level to results data frame
                results$connectivity   [(((loop_id-1)*5)+1):(loop_id*5)] <- connec[conn]
                results$patchquality   [(((loop_id-1)*5)+1):(loop_id*5)] <- ss
                results$dispersal   [(((loop_id-1)*5)+1):(loop_id*5)] <- a
                results$kernel      [(((loop_id-1)*5)+1):(loop_id*5)] <- dd
                results$fluctuation [(((loop_id-1)*5)+1):(loop_id*5)] <- ePeriod
                
                # adding results
                # diversity
                results$div_l [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity(Xd[,-1:-4], index="simpson")
                results$div_r [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-4], 2, sum)), index="simpson")
                # beta diversity
                results$div_b [(((loop_id-1)*5)+1):(loop_id*5)] <- (vegan::diversity((apply(Xd[,-1:-4], 2, sum)), index="simpson") / (vegan::diversity(Xd[,-1:-4], index="simpson"))) / (numCom - 1)
                # richness
                results$rich_l [(((loop_id-1)*5)+1):(loop_id*5)] <- specnumber(Xd[,-1:-4])
                results$rich_r [(((loop_id-1)*5)+1):(loop_id*5)] <- specnumber((apply(Xd[,-1:-4], 2, sum)))
                # evenness
                results$eve_l [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity(Xd[,-1:-4], index="simpson") / log(specnumber(Xd[,-1:-4])) 
                results$eve_r [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-4], 2, sum)), index="simpson") / log(specnumber((apply(Xd[,-1:-4], 2, sum)))) 
                results$eve_l[is.infinite(results$eve_l)] <- 0
                results$eve_r[is.infinite(results$eve_r)] <- 0
                # productivity
                productivity <- Xd[,-1:-4] * mass
                results$prod_l [(((loop_id-1)*5)+1):(loop_id*5)] <- apply(productivity, 1, sum)
                results$prod_r [(((loop_id-1)*5)+1):(loop_id*5)] <- sum(productivity)
                
                # food webs
                # local
                interactionm <- array(NA, dim=c(nSp,nSp,numCom))
                for(j in 1:numCom)
                  for (i in 1:nSp)
                    for(m in 1:nSp)
                      interactionm[i,m,j] <- X[m,j]*BB[i,m]
                # regional
                X_r  <- apply(X, 1, sum)
                intr <- apply(interactionm,1:2,sum)
                # properties
                for(i in 1:numCom){ fw.results[i,] <- FW.results(X[,i], interactionm[,,i]) }
                fw.results.r <- FW.results(X_r, intr)
                # number of links
                results$l_l [(((loop_id-1)*5)+1):(loop_id*5)] <- fw.results[,1]
                results$l_r [(((loop_id-1)*5)+1):(loop_id*5)] <- fw.results.r[1] 
                # link density
                results$ld_l [(((loop_id-1)*5)+1):(loop_id*5)] <- fw.results[,2]
                results$ld_r [(((loop_id-1)*5)+1):(loop_id*5)] <- fw.results.r[2]
                # connectance
                results$c_l [(((loop_id-1)*5)+1):(loop_id*5)] <- fw.results[,3]
                results$c_r [(((loop_id-1)*5)+1):(loop_id*5)] <- fw.results.r[3] 
                # number of trophic levels
                results$t_l [(((loop_id-1)*5)+1):(loop_id*5)] <- fw.results[,4]
                results$t_r [(((loop_id-1)*5)+1):(loop_id*5)] <- fw.results.r[4] 
                
                loop_id <- loop_id + 1
                
              } # finished all replicates for one treatment level
              
              # species abundance time series
              visual <- melt(X_save)
              mainplotitle <- paste("means", r, treatment_id, ".png", sep="")	
              png(filename=mainplotitle, width=29, height=21, units="cm", res=300)
              print(ggplot(visual, aes(x = X3, y = log(value), colour=factor(X2))) +
                stat_summary(fun.y=sum, na.rm=T, geom="line", size=1.2) +
                labs(x = "Time", y = "Abundance") +
                theme_bw())
              dev.off()
              
              id <- treatment_id
              treatment_id <- id+1
              
              break
              
            } # end repeat
          }
        }
      } 
    } 
  }
} # end

