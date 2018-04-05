## metacommunity synchrony

meta.calculate.synchrony <- function(df1){
  per.treatment <- function(df,conn,prop){
    require(codyn);require(reshape2)
    data <- df[df$connectivity==conn & df$quality==prop,]; data <- data[-(2:3)];data <- data[-(3:5)]; data <- data[-(6:10)]
    data <- aggregate(. ~ connectivity * quality * replicate * time, data, sum)
    data <- melt(data[-(1:2)],id.vars=c("time","replicate"),na.rm=TRUE)
    return(synchrony(data,time.var="time",species.var="variable",abundance.var="value",metric="Gross",replicate.var="replicate"))
  }
  synch <- data.frame(matrix(nrow=240, ncol=2))
  rownum = 1
  for (c in 1:length(connec)){
    for (p in ss.prop){
      synch[rownum:(rownum+9),] <- per.treatment(df1,c,p)
      names <- unlist(synch[rownum:(rownum+9),1]) 
      for (i in 1:length(names)){names[i] <- paste0(connec[c],p)}
      synch[rownum:(rownum+9),1] <- names
      rownum = rownum + 10
    }
  }
  return(synch)
}

metacommunity_synchrony <- meta.calculate.synchrony(rep_samples)

for (conn in 1:length(connec)){
  dat <- metacommunity_synchrony[grep(paste0(connec[conn]),x=metacommunity_synchrony$X1),]
  dat$X1 <- factor(dat$X1, levels = paste0(connec[conn],c(0,1,2,3,4,5)))
  pdf(paste0("rawplots/synchrony/", connec[conn],".pdf"))
  boxplot(X2 ~ X1, data = dat, ylab="synchrony", xlab="treatment",main=paste0(connec[conn]))
  dev.off()
}

