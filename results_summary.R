## averages and sd
library(dplyr)
is.na(results)<-sapply(results, is.infinite)
results[is.na(results)]<-NA

## avg and sd of metrics
summary_stats <- results %>% 
  group_by(.dots=c("connectivity","patchquality")) %>% 
  summarise_at(vars(colnames(results[,-1:-5])),funs(mean(.,na.rm=TRUE),sd(.,na.rm=TRUE)))

## avg and sd of species composition
composition <- rep_samples %>% 
  rowwise() %>%
  transform(total = rowSums(rep_samples[12:20]))
composition[12:20] <- composition[12:20]/composition[21][row(composition[12:20])]
composition <- composition[,-21]; composition <- composition[,-(8:11)]
composition_local <- composition %>% 
  group_by(.dots=c("connectivity","quality","time")) %>% 
  summarise_at(vars(colnames(composition[,-(1:7)])),funs(mean(.,na.rm=TRUE),sd(.,na.rm=TRUE)))

