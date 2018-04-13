
# effect size
library(effsize);library(vegan);library(dplyr);#library(ggplot2);detach(package:plyr)

### patchwise differences 

rep_samples <- na.omit(rep_samples) # clean from NAs
rep_samples <- unique(rep_samples)
tran_abundance <- rep_samples[c(1,4:5,8:9,15:23)] # columns we want
tran_abundance[-(1:5)] <- decostand(tran_abundance[-(1:5)], "hellinger", na.rm = TRUE) # standardise
pca_all <- prcomp(tran_abundance[-(1:5)]) #PCA
df.pca <- as.data.frame(pca_all$x)
df.pca <- df.pca[1:2] # only PC1 & PC2
df.pca[3:7] <- tran_abundance[1:5] # identifiers - treatments
df.pca <- df.pca %>% # find centroids
  group_by(time, connectivity, quality) %>%
  mutate(centr_PC1 = mean(PC1), centr_PC2 = mean(PC2))
df.pca <- df.pca %>% # find distance
  mutate(dist1 = sqrt((PC1 - centr_PC1)^2 + (PC2 - centr_PC2)^2)) ## METHOD 1 - within group differences
df.pca$control1 <- rep(df.pca[df.pca$connectivity==1 & df.pca$quality==0,]$centr_PC1, times = 24)
df.pca$control2 <- rep(df.pca[df.pca$connectivity==1 & df.pca$quality==0,]$centr_PC2, times = 24)
df.pca <- df.pca %>% # find distance
  mutate(dist2 = sqrt((PC1 - control1)^2 + (PC2 - control2)^2)) ## METHOD 2 - between group differences

## Within group homogeneity

dists <- data.frame(matrix(NA,nrow=18,ncol=5))
colnames(dists) <- c("treatment","estimate","low","high","magnitude")
id <- 1
for (ss in ss.prop){
  for(conn in 2:4){
    d = as.numeric(cbind(df.pca[df.pca$connectivity == conn & df.pca$quality == ss,]$dist1, # distances
                         df.pca[df.pca$connectivity == 1 & df.pca$quality == 0,]$dist1))
    f = rep(c("Treatment","Control"), each = 550) # treatment vs control (connectivity absent & all low quality)
    co <- cohen.d(d ~ f, paired = TRUE) # pairwise using replicates and time
    dists[id,1] <- paste0(connec[conn],ss) # store
    dists[id,2] <- co$estimate
    dists[id,3:4] <- co$conf.int
    dists[id,5] <- co$magnitude
    id = id + 1
  }
}


## plotting connectivities on different plots
dists[2:5] <- sapply(dists[2:5], as.numeric)
for (conn in 2:4){ # plot
  pdf(paste0("rawplots/effsize/",connec[conn] ,"_withingroup.pdf"))
  dat <- dists[grep(paste0(connec[conn]),x=dists$treatment),]
  plot(0:5, dat$estimate, ylim=range(c(dat$low, dat$high)),
       pch=19, xlab="patch quality", ylab="effect size",
       main=paste0(connec[conn]))
  arrows(0:5, dat$low, 0:5, dat$high,length=0) # confidence interval
  for (p in 0:5){ # magnitude of difference
    if (dat[p+1,]$magnitude==2){points(p + 0.2, dat[p+1,]$estimate,pch=3)}
    else if (dat[p+1,]$magnitude==3){points(p + 0.2, dat[p+1,]$estimate,pch=4)}
    else if (dat[p+1,]$magnitude==4){points(p + 0.2, dat[p+1,]$estimate,pch=8)}
  }
  abline(h=0, lty=2, col = "gray30") # line for 0
  dev.off()
}

## plotting connectivities on same plot
pdf("rawplots/effsize/withingroup.pdf")
plot(NA, NA, xlim=range(c(dists$low, dists$high)),ylim = c(0,17), 
     pch=19, ylab="patch quality", xlab="effect size", yaxt = 'n', main = "Within group effect size")
abline(v = 0, lty=2, col = "gray30") # line for 0
abline(h = c(5.5, 11.5)) # breaking up groups
axis(2, at = 0:17, labels = (rep(0:5, 3))) # correctly label axis
axis(4, at = c(2.5, 8.5, 14.5), labels = c("linear","circular","global"), tick = FALSE) # connectivity
shift <- 0
for (conn in 2:4){ # plot
  dat <- dists[grep(paste0(connec[conn]),x=dists$treatment),]
  points(dat$estimate, 0:5 + shift, xlim=range(c(dat$low, dat$high)),ylim = c(0,17), 
       pch=19, ylab="patch quality", xlab="effect size", yaxt = 'n')
  arrows(dat$low, 0:5 + shift, dat$high, 0:5 + shift, length=0) # confidence interval
  for (p in 0:5){ # magnitude of difference
    if (dat[p+1,]$magnitude==2){points(dat[p+1,]$estimate + 0.01, p - 0.5 + shift, pch=3)}
    else if (dat[p+1,]$magnitude==3){points(dat[p+1,]$estimate + 0.01, p - 0.5 + shift, pch=4)}
    else if (dat[p+1,]$magnitude==4){points(dat[p+1,]$estimate + 0.01, p - 0.5 + shift, pch=8)}
  }
  shift <- shift + 6
}
dev.off()










## Between group homogeneity

dists <- data.frame(matrix(NA,nrow=18,ncol=5))
colnames(dists) <- c("treatment","estimate","low","high","magnitude")
id <- 1
for (ss in ss.prop){
  for(conn in 2:4){
    d = as.numeric(cbind(df.pca[df.pca$connectivity == conn & df.pca$quality == ss,]$dist2, # distances
                         df.pca[df.pca$connectivity == 1 & df.pca$quality == 0,]$dist2))
    f = rep(c("Treatment","Control"), each = 550) # treatment is treatment in question, control is no connectivity & all low quality
    co <- cohen.d(d ~ f, paired = TRUE) # pairwise using replicates & times
    dists[id,1] <- paste0(connec[conn],ss)
    dists[id,2] <- co$estimate
    dists[id,3:4] <- co$conf.int
    dists[id,5] <- co$magnitude
    id = id + 1
  }
}

## plotting connectivities on different plots
dists[2:5] <- sapply(dists[2:5], as.numeric)
for (conn in 2:4){ # plot
  pdf(paste0("rawplots/effsize/",connec[conn] ,"_betweengroup.pdf"))
  dat <- dists[grep(paste0(connec[conn]),x=dists$treatment),]
  plot(0:5, dat$estimate, ylim=range(c(dat$low, dat$high)),
       pch=19, xlab="patch quality", ylab="effect size",
       main=paste0(connec[conn]))
  arrows(0:5, dat$low, 0:5, dat$high,length=0) # confidence interval
  for (p in 0:5){ # magnitude
    if (dat[p+1,]$magnitude==2){points(p + 0.2, dat[p+1,]$estimate,pch=3)}
    else if (dat[p+1,]$magnitude==3){points(p + 0.2, dat[p+1,]$estimate,pch=4)}
    else if (dat[p+1,]$magnitude==4){points(p + 0.2, dat[p+1,]$estimate,pch=8)}
  }
  abline(h=0, lty=2, col = "gray30")
  dev.off()
}

## plotting connectivities on same plot
pdf("rawplots/effsize/betweengroup.pdf")
plot(NA, NA, xlim=range(c(dists$low, dists$high)),ylim = c(0,17), 
     pch=19, ylab="patch quality", xlab="effect size", yaxt = 'n', main = "Between group effect size")
abline(v = 0, lty=2, col = "gray30") # line for 0
abline(h = c(5.5, 11.5)) # breaking up groups
axis(2, at = 0:17, labels = (rep(0:5, 3))) # correctly label axis
axis(4, at = c(2.5, 8.5, 14.5), labels = c("linear","circular","global"), tick = FALSE) # connectivity
shift <- 0
for (conn in 2:4){ # plot
  dat <- dists[grep(paste0(connec[conn]),x=dists$treatment),]
  points(dat$estimate, 0:5 + shift, xlim=range(c(dat$low, dat$high)),ylim = c(0,17), 
         pch=19, ylab="patch quality", xlab="effect size", yaxt = 'n')
  arrows(dat$low, 0:5 + shift, dat$high, 0:5 + shift, length=0) # confidence interval
  for (p in 0:5){ # magnitude of difference
    if (dat[p+1,]$magnitude==2){points(dat[p+1,]$estimate + 0.01, p - 0.5 + shift, pch=3)}
    else if (dat[p+1,]$magnitude==3){points(dat[p+1,]$estimate + 0.01, p - 0.5 + shift, pch=4)}
    else if (dat[p+1,]$magnitude==4){points(dat[p+1,]$estimate + 0.01, p - 0.5 + shift, pch=8)}
  }
  shift <- shift + 6
}
dev.off()

