## effect size & coefficient of variation

# effect s9ze
library(effsize)

# comps possible - comparing between connectivities, controls as comparison, using diversity ????
regional <- aggregate(div_r ~ connectivity * patchquality * r, data = results, mean)
dists <- data.frame(matrix(NA,nrow=18,ncol=5))
colnames(dists) <- c("treatment","estimate","low","high","magnitude")
id <- 1
for (ss in ss.prop){
  for(conn in 2:4){
    d = as.numeric(cbind(regional[regional$connectivity == connec[conn] & regional$patchquality == ss,]$div_r
                         ,regional[regional$connectivity == connec[conn] & regional$patchquality == 0,]$div_r))
    f = rep(c("Treatment","Control"), each = 10)
    co <- cohen.d(d ~ f)
    dists[id,1] <- paste0(connec[conn],ss)
    dists[id,2] <- co$estimate
    dists[id,3:4] <- co$conf.int
    dists[id,5] <- co$magnitude
    id <- id + 1
  }
}
is.na(dists)<-sapply(dists, is.nan); dists[is.na(dists)]<-0; dists[dists$magnitude==0,]$magnitude <- 1
dists[2:5] <- sapply(dists[2:5], as.numeric)
for (conn in 2:4){
  pdf(paste0("rawplots/effsize/",connec[conn] ,"3.pdf"))
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


