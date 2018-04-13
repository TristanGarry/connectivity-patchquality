
is.na(results)<-sapply(results, is.infinite)
results[is.na(results)]<-NA

results$connectivity<-factor(results$connectivity, levels=c("d.a","d.l","d.c","d.g")) # properly order

metrics <- colnames(results[,10:26])

for (i in 1:length(metrics)){
  pdf(paste0("rawplots/metrics/",metrics[i] ,"_env_",enveff,".pdf"))
  boxplot(as.formula(paste(metrics[i], " ~ patchquality * connectivity")), data=results,col=3:8
          ,xlab="patch quality x connectivity",ylab=paste(metrics[i]),outline=FALSE)
  abline(v=6.5);abline(v=12.5);abline(v=18.5)
  dev.off()
}


