
# effect size
library(effsize);library(vegan);library(dplyr);library(ggplot2)#;detach(package:plyr)

### patchwise differences 

rep_samples <- na.omit(rep_samples) # clean from NAs
rep_samples <- unique(rep_samples)
tran_abundance <- rep_samples[c(1,4:5,8:9,15:23)] # columns we want
tran_abundance[-(1:5)] <- decostand(tran_abundance[-(1:5)], "hellinger", na.rm = TRUE) # standardise
pca_all <- prcomp(tran_abundance[-(1:5)]) #PCA
df.pca <- as.data.frame(pca_all$x)
df.pca <- df.pca[1:2] # only PC1 & PC2
df.pca[3:7] <- tran_abundance[1:5] # identifiers
df.pca <- df.pca %>% # find centroids
  group_by(time, connectivity, quality) %>%
  mutate(centr_PC1 = mean(PC1), centr_PC2 = mean(PC2))
df.pca <- df.pca %>% # find distance
  mutate(dist1 = sqrt((PC1 - centr_PC1)^2 + (PC2 - centr_PC2)^2)) ## METHOD 1
#df.pca <- df.pca[df.pca$time == 500,]
df.pca$control1 <- rep(df.pca[df.pca$connectivity==1 & df.pca$quality==0,]$centr_PC1, times = 24)
df.pca$control2 <- rep(df.pca[df.pca$connectivity==1 & df.pca$quality==0,]$centr_PC2, times = 24)
df.pca <- df.pca %>% # find distance
  mutate(dist2 = sqrt((PC1 - control1)^2 + (PC2 - control2)^2)) ## METHOD 2

## Within group homogeneity

dists <- data.frame(matrix(NA,nrow=18,ncol=5))
colnames(dists) <- c("treatment","estimate","low","high","magnitude")
id <- 1
for (ss in ss.prop){
  for(conn in 2:4){
    d = as.numeric(cbind(df.pca[df.pca$connectivity == conn & df.pca$quality == ss,]$dist1,
                         df.pca[df.pca$connectivity == 1 & df.pca$quality == 0,]$dist1))
    f = rep(c("Treatment","Control"), each = 550)
    co <- cohen.d(d ~ f, paired = TRUE)
    dists[id,1] <- paste0(connec[conn],ss)
    dists[id,2] <- co$estimate
    dists[id,3:4] <- co$conf.int
    dists[id,5] <- co$magnitude
    id = id + 1
  }
}

dists[2:5] <- sapply(dists[2:5], as.numeric)
for (conn in 2:4){
  pdf(paste0("rawplots/effsize/",connec[conn] ,"_withingroup.pdf"))
  dat <- dists[grep(paste0(connec[conn]),x=dists$treatment),]
  plot(0:5, dat$estimate, ylim=range(c(dat$low, dat$high)),
       pch=19, xlab="patch quality", ylab="effect size",
       main=paste0(connec[conn]))
  arrows(0:5, dat$low, 0:5, dat$high,length=0)
  for (p in 0:5){
    if (dat[p+1,]$magnitude==2){points(p + 0.2, dat[p+1,]$estimate,pch=3)}
    else if (dat[p+1,]$magnitude==3){points(p + 0.2, dat[p+1,]$estimate,pch=4)}
    else if (dat[p+1,]$magnitude==4){points(p + 0.2, dat[p+1,]$estimate,pch=8)}
  }
  abline(h=0, lty=2, col = "gray30")
  dev.off()
}




## Between group homogeneity

dists <- data.frame(matrix(NA,nrow=18,ncol=5))
colnames(dists) <- c("treatment","estimate","low","high","magnitude")
id <- 1
for (ss in ss.prop){
  for(conn in 2:4){
    d = as.numeric(cbind(df.pca[df.pca$connectivity == conn & df.pca$quality == ss,]$dist2,
                         df.pca[df.pca$connectivity == 1 & df.pca$quality == 0,]$dist2))
    f = rep(c("Treatment","Control"), each = 550)
    co <- cohen.d(d ~ f, paired = TRUE)
    dists[id,1] <- paste0(connec[conn],ss)
    dists[id,2] <- co$estimate
    dists[id,3:4] <- co$conf.int
    dists[id,5] <- co$magnitude
    id = id + 1
  }
}

dists[2:5] <- sapply(dists[2:5], as.numeric)
for (conn in 2:4){
  pdf(paste0("rawplots/effsize/",connec[conn] ,"_betweengroup.pdf"))
  dat <- dists[grep(paste0(connec[conn]),x=dists$treatment),]
  plot(0:5, dat$estimate, ylim=range(c(dat$low, dat$high)),
       pch=19, xlab="patch quality", ylab="effect size",
       main=paste0(connec[conn]))
  arrows(0:5, dat$low, 0:5, dat$high,length=0)
  for (p in 0:5){
    if (dat[p+1,]$magnitude==2){points(p + 0.2, dat[p+1,]$estimate,pch=3)}
    else if (dat[p+1,]$magnitude==3){points(p + 0.2, dat[p+1,]$estimate,pch=4)}
    else if (dat[p+1,]$magnitude==4){points(p + 0.2, dat[p+1,]$estimate,pch=8)}
  }
  abline(h=0, lty=2, col = "gray30")
  dev.off()
}


