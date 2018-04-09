

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

## coefficient of variation - abundances

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

temp <- rep_samples[-c(2:3,6:7,10:14)] # clean and tranform data
temp[-(1:5)] <- decostand(temp[-(1:5)], "hellinger", na.rm = TRUE)
coef_var <- temp %>% # mean and sd for each patch
  group_by(.dots=c("connectivity","quality","replicate","time","patch")) %>% 
  transform(mean = rowMeans(temp[6:14])) %>%
  transform(sd = apply(temp[6:14], 1, sd))
coef_var <- coef_var %>%  # coefficient of variation for each patch
  transform(coefvar = coef_var[16]/coef_var[15])
coef_var <- coef_var[coef_var$time == 500,-c(1,6:16)] # clean data
coef_var_sum <- summarySE(data = coef_var, measurevar = "sd.1" # variation in coefficient of variation for plotting
                          , groupvars = c("patch","connectivity","quality"), na.rm = TRUE)
pd <- position_dodge(width = 0.5)
for (conn in 1:length(connec)){ # plot by connectivity & save
  pdf(paste0("rawplots/perpatch_coefvar/",connec[conn] ,".pdf"))
  print(ggplot(coef_var_sum[coef_var_sum$connectivity==conn,], aes(x=patch, y=sd.1, group=factor(quality), colour=factor(quality))) + 
    geom_errorbar(aes(ymin=sd.1-se, ymax=sd.1+se, width=1), position=pd) +
    geom_line(position=pd, size=1.2) +
    geom_point(aes(x=patch, y=sd.1, group=factor(quality), colour=factor(quality)), position=pd, size = 3) + 
    ylab("coefficient of variation of abundances (sd/mean)") + 
    labs(title = paste0(connec[conn])) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(xintercept=c(1.5,2.5,3.5,4.5), colour='black', size = .1))
  dev.off()
}

for (ss in ss.prop){ # plot by patch quality & save
  pdf(paste0("rawplots/perpatch_coefvar/pc_",ss ,".pdf"))
  print(ggplot(coef_var_sum[coef_var_sum$quality==ss,], aes(x=patch, y=sd.1, group=factor(connectivity), colour=factor(connectivity))) + 
          geom_errorbar(aes(ymin=sd.1-se, ymax=sd.1+se, width=1), position=pd) +
          geom_line(position=pd, size=1.2) +
          geom_point(aes(x=patch, y=sd.1, group=factor(connectivity), colour=factor(connectivity)), position=pd, size = 3) + 
          ylab("coefficient of variation of abundances (sd/mean)") + 
          labs(title = paste0("patchquality=",ss)) + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
          geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5), colour='black', size = .1))
  dev.off()
}
## coefficient of variation - diversity - not finished

temp <- results[c(3:6,10)]
temp <- aggregate(temp[5], by = temp[1:4], mean)
coef_var <- temp %>%
  group_by(connectivity,patchquality,r,patch) %>% 
  transform(mean = mean(div_l)) %>%
  transform(sd = sd(div_l))

coef_var <- temp %>%
  group_by(.dots=c("connectivity","quality","replicate","time","patch")) %>% 
  transform(mean = rowMeans(temp[6:14])) %>%
  transform(sd = apply(temp[6:14], 1, sd))

coef_var <- coef_var %>%  
  transform(coefvar = coef_var[6]/coef_var[7])
coef_var_sum <- summarySE(data = coef_var, measurevar = "mean.1"
                          , groupvars = c("connectivity","patchquality","patch"), na.rm = TRUE)
pd <- position_dodge(width = 0.5)
for (conn in 1:length(connec)){
  pdf(paste0("rawplots/perpatch_coefvar/div_",connec[conn] ,".pdf"))
  print(ggplot(coef_var_sum[coef_var_sum$connectivity==connec[conn],], aes(x=patch, y=mean.1, group=factor(patchquality), colour=factor(patchquality))) + 
          geom_errorbar(aes(ymin=mean.1-se, ymax=mean.1+se, width=1), position=pd) +
          geom_line(position=pd, size=1.2) +
          geom_point(aes(x=patch, y=mean.1, group=factor(patchquality), colour=factor(patchquality)), position=pd, size = 3) + 
          ylab("coefficient of variation of diversity (sd/mean)") + 
          labs(title = paste0(connec[conn])) + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
          geom_vline(xintercept=c(1.5,2.5,3.5,4.5), colour='black', size = .1))
  dev.off()
}








