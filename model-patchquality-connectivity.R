###############################################
#           SESIM - microcosm model           #
# for patch quality x connectivity experiment #
###############################################

library(vegan)
library(igraph)
library(reshape)
library(ggplot2)

## Setting up metacommunity

### Connectivity

# number of patches
numCom <- 5
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

# dispersal function
calc.immigration <- function(Nd,a,dispersal_matrix)	{
  immigrants <- dispersal_matrix%*%Nd*a
  return(immigrants)
}

# dispersal rate
rate <- 0.001

# species-specific dispersal ability
disp <- c(0.13, 0.13, 0.13, 0.64, 0.62, 0.53, 0.90, 0.62, 0.44, 0.33, 0.71, 1.00)

# dispersal kernel
k <- 1

### Patch quality

# growth rates
r.source <- c(8.4 ,  16.8,   20.3,  0, 0.92, 0,  0, 2.63,  0,  0, 1.78, 0)
r.sink   <- c(0.08,  0.16,   0.20,  0, 0.16, 0,  0, 1.25,  0,  0, 2.04, 0)

# patch quality
ss.prop <- c(0,1,2,3,4,5)

### Species

# food web
FW <- read.csv(file="foodweb.csv", header=T)
# number of species
nSp <- nrow(FW)

mass <- c(4.3e-08,1.96e-07,3.76e-06,1.52e-08,4.77e-09,6.9e-08,9.68e-08,8.05e-08,2.27e-07)

# number of environmental variables
nEnv <- 1
# environmental niche amplitude
eA <- 1
# strength of environmental effect
enveff <- 1

# environmental fluctuation period
eP <- 0

### Experimental design

# experiment factorial design
experiment <- data.frame(connectivity = rep(NA, each=length(connec)*length(ss.prop)), 
                         patchquality = rep(NA))

### Simulation parameters

# number of replicates
reps <- 10

# simulation length
Tmax <- 1000

# sampling frequency
sampfreq <- 100
sampleV <- seq(0, Tmax, by=sampfreq)
# empty array for saving sampled data
X_save <- array(data=NA, dim=c(numCom, nSp, Tmax/sampfreq))

# data frames for results
results <- data.frame(treatment = rep(1:(dim(experiment)[1]), each=numCom*reps),
                      r = rep(1:reps, each=numCom),
                      patch = rep(1:numCom))

# Model

treatment_id = 1
loop_id = 1

for(conn in seq_along(connec)){
  for(ss in ss.prop){
    repeat{
      for(r in 1:reps){
        cc <- get(connec[conn])
        sampling = 1
        
        ## patch quality
        C <- matrix(c(rep(r.source, times=ss), rep(r.sink, times=(5-ss))), nrow=12, ncol=5)
        
        # matrix of species interaction
        BB  <- data.matrix(as.data.frame(FW[,-1]))
        matriz  <- BB
        novamatriz <- matriz
        for(j in 1:ncol(matriz))
          for(i in 1:nrow(matriz))
            novamatriz[i,j] <- runif(1, min=min(0, matriz[i,j]), max=max(0, matriz[i,j]))
        BB <- novamatriz
        
        # dispersal ability
        d_exp <- exp(-k*cc) - diag(nrow(cc))
        dispersal_m <- apply(d_exp, 1, function(x) x/sum(x))
        dispersal_m[is.na(dispersal_m)] <- 0
        
        # species environmental optimum
        Env_Opt <- matrix(runif(nSp, min=0, max=eA), nSp, nEnv)
        # environmental condition in each patch per species
        Env <- matrix(rep(0.5), nrow=nSp, ncol=numCom)
        # environmental match
        enviro <- (Env-Env_Opt[,1])*enveff
        
        # species initial abundances
        X   <- data.matrix(read.csv(file="initialabundX.csv", header=T) [,-1])
        Xd  <- t(X)
        
        for(l in 1:(Tmax)){
          
          # interactions
          interactions <- ((BB + t(BB)) %*% (X*(C+enviro))) / ((BB + t(BB)) %*% (X)) 
          interactions[X == 0] <- 0
          # number of migrants in each patch per species
          Immigrants <- t(calc.immigration(Xd, rate, dispersal_m))
          Migrants   <- apply(X, 2, function(x) x * (rate*disp))
          # Lotka-Volterra model
          Xt <- X + interactions + Immigrants - Migrants
          
          Xhold <- Xt
          Xhold[(Xhold < 0.1)] <- 0
          X <- Xhold
          Xd <- t(X)
          if(all(is.nan(X))) break
          
          if(l==sampleV[l/sampfreq+1] && l<=Tmax){
            X_save[,,sampling] <- Xd
            sampling <- sampling + 1
            
          } # end sampling
        } # end all time steps
        
        # saving sampled data for one replicate and one treatment level
        treat <- data.frame(time=(rep(seq_len(dim(X_save)[3]),each=dim(X_save)[1]))*sampfreq, 
                            dispersal = rep(rate), 
                            kernel = rep(k), 
                            replicate = rep(r), 
                            patch = rep(1:5),
                            connectivity = rep(conn),
                            quality = rep(ss))
        X_saved <- as.data.frame(apply(X_save, 2, cbind))
        assign(paste("rep_sampled_", loop_id, "_", r, sep=""), cbind(treat, X_saved))
        
        # adding treatment level to experiment 
        experiment$connectivity [treatment_id] <- connec[conn]
        experiment$patchquality [treatment_id] <- ss
        
        # adding results
        results$connectivity [(((loop_id-1)*5)+1):(loop_id*5)] <- connec[conn]
        results$patchquality [(((loop_id-1)*5)+1):(loop_id*5)] <- ss
        # diversity
        results$div_l [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity(Xd[,-1:-3], index="simpson")
        results$div_r [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-3], 2, sum)),
                                                                            index="simpson")
        # beta diversity
        results$div_b [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-3], 2, sum)), 
                                                                            index="simpson") -
          mean(vegan::diversity(Xd[,-1:-3], index="simpson"))
        # richness
        results$rich_l [(((loop_id-1)*5)+1):(loop_id*5)] <- specnumber(Xd[,-1:-3])
        results$rich_r [(((loop_id-1)*5)+1):(loop_id*5)] <- specnumber((apply(Xd[,-1:-3], 2, sum)))
        # evenness
        results$eve_l [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity(Xd[,-1:-3], index="simpson") /
          log(specnumber(Xd[,-1:-3])) 
        results$eve_r [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-3], 2, sum)), 
                                                                            index="simpson") /
          specnumber((apply(Xd[,-1:-3], 2, sum))) 
        # productivity
        productivity <- Xd[,-1:-3] * mass
        results$prod_l [(((loop_id-1)*5)+1):(loop_id*5)] <- apply(productivity, 1, sum)
        results$prod_r [(((loop_id-1)*5)+1):(loop_id*5)] <- sum(productivity)
        
        loop_id <- loop_id + 1
        
      } # finished all replicates for one treatment level
      
      # regional species abundance time series
      visual <- melt(X_save)
      mainplotitle <- paste("abund_replicate", r, treatment_id, ".png", sep="")	
      png(filename=mainplotitle, width=29, height=21, units="cm", res=300)
      ggplot(visual, aes(x = X3, y = log(value), colour=factor(X2))) +
        geom_point() + 
        labs(x = "Time", y = "Abundance") +
        theme_bw()
      dev.off()
      
      id <- treatment_id
      treatment_id <- id+1
      
      break
      
    } # end repeat
  } 
} # end

write.csv(results, file="cp-results.csv")

names_samples <- grep("rep_sampled_", x=ls(), value=T)
rep_samples <- do.call(rbind, mget(names_samples))
str(rep_samples)
write.csv(rep_samples, file="cp-timeseries.csv")


