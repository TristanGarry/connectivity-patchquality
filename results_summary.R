
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


for (c in 1:length(connec)){
  for (q in ss.prop){
    png(paste0("composition_graphs/composition_box_",connec[c] ,q ,".png"))
    boxplot(x = composition[composition$connectivity==c & composition$quality==q & composition$time == 1000,-(1:7)]
            ,xlab="species",ylab="composition",outline=FALSE)
    dev.off()
    png(paste0("composition_graphs/composition_bar_",c ,q ,".png"))
    barplot()
      
    dev.off()
  }
}
