

## per patch effects - * not presentable yet * 

is.na(results)<-sapply(results, is.infinite)
results[is.na(results)]<-NA

results$connectivity<-factor(results$connectivity, levels=c("d.a","d.l","d.c","d.g"))

metrics <- colnames(results[,10:26])

for (i in 1:length(metrics)){
  pdf(paste0("rawplots/metrics/",metrics[i] ,".pdf"))
  boxplot(as.formula(paste(metrics[i], " ~ patchquality * connectivity * patch")), data=results #,col=3:8
          ,xlab="patch quality x connectivity",ylab=paste(metrics[i]),outline=FALSE)
  dev.off()
}

## coefficient of variation

library(ggplot2)

# found online http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

library(dplyr); library(reshape)

temp <- rep_samples[-(2:3)];temp <- temp[-(4:5)]; temp <- temp[-(6:10)]
coef_var <- temp %>%
  group_by(.dots=c("connectivity","quality","replicate","time","patch")) %>% 
  transform(mean = rowMeans(temp[6:14])) %>%
  transform(sd = apply(temp[6:14], 1, sd))
coef_var <- coef_var %>%  
  transform(coefvar = coef_var[16]/coef_var[15])
coef_var <- coef_var[coef_var$time == 500,]; coef_var <- coef_var[-1]; coef_var <- coef_var[-(5:15)]
coef_var_sum <- summarySE(data = coef_var, measurevar = "sd.1"
                          , groupvars = c("patch","connectivity","quality"), na.rm = TRUE)
pd <- position_dodge(0.1)
for (conn in 1:length(connec)){
  pdf(paste0("rawplots/perpatch_coefvar/",connec[conn] ,".pdf"))
  print(ggplot(coef_var_sum[coef_var_sum$connectivity==conn,], aes(x=patch, y=sd.1, group=factor(quality), colour=factor(quality))) + 
    geom_errorbar(aes(ymin=sd.1-se, ymax=sd.1+se), width=.3, position=pd) +
    geom_line(position=pd, size=1) +
    geom_point(position=pd, size=1) + 
    ylab("coefficient of variation") + 
    labs(title = paste0(connec[conn])))
  dev.off()
}


## plots to do
# CCA 
# synchrony 
# CCA distances/differences from centroid







