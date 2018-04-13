## metacommunity synchrony - overall

meta.calculate.synchrony <- function(df1){
  per.treatment <- function(df,conn,prop){ # function to be applied to each treatment
    require(codyn);require(reshape2)
    data <- df[df$connectivity==conn & df$quality==prop,]; data <- data[c(1,4,8:9,15:23)] # subset to treatment
    data <- aggregate(. ~ connectivity * quality * replicate * time, data, sum) # sum metacommunity abundance
    data <- melt(data[-(1:2)],id.vars=c("time","replicate"),na.rm=TRUE) # make data usable for synchrony
    return(synchrony(data,time.var="time",species.var="variable",abundance.var="value", # within metacommunity synchrony
                     metric="Gross",replicate.var="replicate")) # with replicates as replicates
  }
  synch <- data.frame(matrix(nrow=240, ncol=2))
  rownum = 1
  for (c in 1:length(connec)){
    for (p in ss.prop){
      synch[rownum:(rownum+9),] <- per.treatment(df1,c,p) # apply synchrony and store
      names <- unlist(synch[rownum:(rownum+9),1]) 
      for (i in 1:length(names)){names[i] <- paste0(connec[c],p)}
      synch[rownum:(rownum+9),1] <- names
      rownum = rownum + 10
    }
  }
  return(synch)
}

metacommunity_synchrony <- meta.calculate.synchrony(rep_samples) # apply function to abundances

for (conn in 1:length(connec)){ # plot
  dat <- metacommunity_synchrony[grep(paste0(connec[conn]),x=metacommunity_synchrony$X1),]
  dat$X1 <- factor(dat$X1, levels = paste0(connec[conn],c(0,1,2,3,4,5)))
  pdf(paste0("rawplots/synchrony/", connec[conn],".pdf"))
  boxplot(X2 ~ X1, data = dat, ylab="synchrony", xlab="treatment",main=paste0(connec[conn]))
  dev.off()
}

## metacommunity synchrony - between patches

meta.patch.calculate.synchrony <- function(df1){
  per.treatment.patch <- function(df,conn,prop,r){ # function to be applied to each treatment & replicate
    require(codyn);require(reshape2)
    data <- df[df$connectivity==conn & df$quality==prop & df$replicate==r,] # subset by treatment & replicate
    #data <- unique(data[c(1,5,8:9,15:23)]) # not necessary anymore - there was an issue in double entries
    data <- melt(data[-(3:4)],id.vars=c("time","patch"),na.rm=TRUE) # data usable for synchrony
    return(synchrony(data,time.var="time",species.var="variable",abundance.var="value", # within patch synchrony
                     metric="Gross",replicate.var="patch")) # with patches as replicates
  }
  synch <- data.frame(matrix(nrow=2400, ncol=2))
  rownum = 1
  for (c in 1:length(connec)){
    for (p in ss.prop){
      for (r in 1:10){ # apply over each treatment & replicate and store
        synch[rownum:(rownum+9),] <- per.treatment.patch(df1,c,p,r)
        names <- unlist(synch[rownum:(rownum+9),1]) 
        for (i in 1:length(names)){names[i] <- paste0(connec[c],p)}
        synch[rownum:(rownum+9),1] <- names
        rownum = rownum + 10
      }
    }
  }
  return(synch)
}
metacommunity_patch_synchrony <- meta.patch.calculate.synchrony(rep_samples) # apply 
#metacommunity_patch_synchrony[!is.finite(unlist(metacommunity_patch_synchrony[2])),2] <- NA

for (conn in 1:length(connec)){ # plot
  dat <- metacommunity_patch_synchrony[grep(paste0(connec[conn]),x=metacommunity_patch_synchrony$X1),]
  dat$X1 <- factor(dat$X1, levels = paste0(connec[conn],c(0,1,2,3,4,5)))
  pdf(paste0("rawplots/synchrony/patch_", connec[conn],".pdf"))
  boxplot(X2 ~ X1, data = dat, ylab="patch synchrony", xlab="treatment",main=paste0(connec[conn]),
          ylim = c(-1.0,1.0))
  dev.off()
}








