
## averages and sd
library(dplyr); library(ggplot2); library(reshape)
is.na(results)<-sapply(results, is.infinite) # clean results
results[is.na(results)]<-NA

## avg and sd of metrics
summary_stats <- results %>% 
  group_by(.dots=c("connectivity","patchquality")) %>% # by treatment
  summarise_at(vars(colnames(results[,-1:-8])),funs(mean(.,na.rm=TRUE),sd(.,na.rm=TRUE)))

## avg and sd of species composition
composition <- rep_samples %>% 
  rowwise() %>%
  transform(total = rowSums(rep_samples[15:23])) # total abundance
composition[15:23] <- composition[15:23]/composition[24][row(composition[15:23])] # scale to find relative abundance
composition <- composition[,-24]; composition <- composition[,-(11:14)] # remove extra columns
composition_local <- composition %>% 
  group_by(.dots=c("connectivity","quality","time")) %>% # find avgs and sd of compositions
  summarise_at(vars(colnames(composition[,-(1:10)])),funs(mean(.,na.rm=TRUE),sd(.,na.rm=TRUE)))


# plotting community composition

for (c in 1:length(connec)){ # plot
  for (q in ss.prop){
    png(paste0("rawplots/communitycomp/composition_box_",connec[c] , q,".png"))
    boxplot(x = composition[composition$connectivity==c & composition$quality==q & composition$time == max(composition$time - sampfreq),-(1:10)]
            ,xlab="species",ylab="composition",outline=FALSE)
    dev.off()
  }
}

