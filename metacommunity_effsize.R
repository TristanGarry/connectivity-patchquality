## effect size & coefficient of variation

# effect size
library(effsize);library(vegan);library(dplyr);library(ggplot2)

rep_samples <- na.omit(rep_samples) # clean from NAs
tran_abundance <- rep_samples[c(1,4:5,8:9,15:23)] # columns we want
tran_abundance <- aggregate(tran_abundance[6:14], by = tran_abundance[(c(1:2,4:5))], sum) # effect size is 
                                                                        #by metacommunity rather than patch
tran_abundance[-(1:4)] <- decostand(tran_abundance[-(1:4)], "hellinger", na.rm = TRUE) # standardise
pca_all <- prcomp(tran_abundance[-(1:4)]) #PCA
df.pca <- as.data.frame(pca_all$x)
df.pca <- df.pca[1:2] # only PC1 & PC2
df.pca[3:6] <- tran_abundance[1:4] # identifiers
df.pca <- df.pca %>% # find centroids
  group_by(time, connectivity, quality) %>%
  mutate(centr_PC1 = mean(PC1), centr_PC2 = mean(PC2))
df.pca <- df.pca %>% # find distance
  mutate(dist = sqrt((PC1 - centr_PC1)^2 + (PC2 - centr_PC2)^2))
#df.pca <- df.pca[df.pca$time == 500,]

# dists <- data.frame(matrix(NA,nrow=18,ncol=5))
# colnames(dists) <- c("treatment","estimate","low","high","magnitude")
# id <- 1
# for (ss in ss.prop){
#   for(conn in 2:4){
#     dat <- df.pca[df.pca$connectivity == conn & df.pca$quality == ss |
#                     df.pca$connectivity == 1 & df.pca$quality == 0,]
#     c1 <- dat[1,]$centr_PC1;c2 <- dat[1,]$centr_PC2
#     dat <- dat %>%
#       mutate(dist = sqrt((PC1 - c1)^2 + (PC2 - c2)^2))
#     d = as.numeric(cbind(df.pca[df.pca$connectivity == conn & df.pca$quality == ss,]$dist
#                          ,df.pca[df.pca$connectivity == 1 & df.pca$quality == 0,]$dist))
#     f = rep(c("Treatment","Control"), each = 10)
#     co <- cohen.d(d ~ f)
#     dists[id,1] <- paste0(connec[conn],ss)
#     dists[id,2] <- co$estimate
#     dists[id,3:4] <- co$conf.int
#     dists[id,5] <- co$magnitude
#     id <- id + 1
#   }
# }


## effsize by timestep

dists <- data.frame(matrix(NA,nrow=198,ncol=6)); colnames(dists) <- c("treatment","estimate","low","high","magnitude","time")
id <- 1
for (conn in 2:4){
  dat <- df.pca[df.pca$connectivity == conn | df.pca$connectivity == 1,] # subset
  for (t in unique(dat$time)[-1]){
    dat2 <- dat[dat$time==t,] # subset
    for (ss in ss.prop){
      dat3 <- dat2[dat2$quality == ss | dat2$quality == 0,] # subset
      d = as.numeric(cbind(dat3[dat3$connectivity == conn & dat3$quality == ss,]$dist # isolate distances
                           ,dat3[dat3$connectivity == 1 & dat3$quality == 0,]$dist))
      f = rep(c("Treatment","Control"), each = 10)
      co <- cohen.d(d ~ f)
      dists[id,1] <- paste0(connec[conn],ss)
      dists[id,2] <- co$estimate
      dists[id,3:4] <- co$conf.int
      dists[id,5] <- co$magnitude
      dists[id,6] <- t
      id <- id + 1
    }
  }
}

## plotting is for timestep 500 only - ignore for now

is.na(dists)<-sapply(dists, is.nan); dists[is.na(dists)]<-0; dists[dists$magnitude==0,]$magnitude <- 1
dists[2:5] <- sapply(dists[2:5], as.numeric)
for (conn in 2:4){
  pdf(paste0("rawplots/effsize/",connec[conn] ,"test2.pdf"))
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


