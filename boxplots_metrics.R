
is.na(results)<-sapply(results, is.infinite)
results[is.na(results)]<-NA

results$connectivity<-factor(results$connectivity, levels=c("d.a","d.l","d.c","d.g"))

metrics <- colnames(results[,6:14])

for (i in 1:length(metrics)){
  png(paste0("metric_graphs/",metrics[i] ,"_env_",enveff,".png"))
  boxplot(as.formula(paste(metrics[i], " ~ patchquality * connectivity")), data=results,col=3:8
          ,xlab="patch quality x connectivity",ylab=paste(metrics[i]),outline=FALSE)
  abline(v=6.5);abline(v=12.5);abline(v=18.5)
  dev.off()
}


