## assembly trajectory - metacommunity * needs more computing power *


require(vegan); require(ggplot2)

test <- na.omit(rep_samples)
test$interaction <- paste0(test$time,"-",test$connectivity,"-",test$quality)
test <- test[,15:24]

dist <- vegdist(test[,1:9],method="bray")
disp <- betadisper(dist, group = test[,]$interaction, sqrt.dist=F, type = "centroid")
temp <- cbind(na.omit(rep_samples[,1:10])[,],scores(disp)$sites)
test <- cbind(test, na.omit(rep_samples[,1:10]))
ggplot(data = temp[temp$connectivity==1 & temp$quality==0 & temp$patch==1,],aes(x=PCoA1,y=PCoA2,group=replicate,by=time)) +
  geom_path(colour=col,size=1) + geom_point(data=temp[temp$treatment==level & temp$time==30,],size=4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



