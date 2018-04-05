
## averages and sd
library(dplyr)
library(ggplot2)
library(reshape)
is.na(results)<-sapply(results, is.infinite)
results[is.na(results)]<-NA

## avg and sd of metrics
summary_stats <- results %>% 
  group_by(.dots=c("connectivity","patchquality")) %>% 
  summarise_at(vars(colnames(results[,-1:-8])),funs(mean(.,na.rm=TRUE),sd(.,na.rm=TRUE)))

## avg and sd of species composition
composition <- rep_samples %>% 
  rowwise() %>%
  transform(total = rowSums(rep_samples[15:23]))
composition[15:23] <- composition[15:23]/composition[24][row(composition[15:23])]
composition <- composition[,-24]; composition <- composition[,-(11:14)]
composition_local <- composition %>% 
  group_by(.dots=c("connectivity","quality","time")) %>% 
  summarise_at(vars(colnames(composition[,-(1:10)])),funs(mean(.,na.rm=TRUE),sd(.,na.rm=TRUE)))


# plotting community composition

for (c in 1:length(connec)){
  for (q in ss.prop){
    png(paste0("rawplots/communitycomp/composition_box_",connec[c] , q,".png"))
    boxplot(x = composition[composition$connectivity==c & composition$quality==q & composition$time == max(composition$time - sampfreq),-(1:10)]
            ,xlab="species",ylab="composition",outline=FALSE)
    dev.off()
    #png(paste0("tests/newmodel/composition_bar_",connec[c] , q,".png"))
    #barplot(height = t(matrix(composition_local[composition_local$connectivity == c & composition_local$quality == q & composition_local$time == max(composition$time - sampfreq),11:19]))
    #        ,names = colnames(composition[,-(1:10)]))
    #dev.off()
  }
}

