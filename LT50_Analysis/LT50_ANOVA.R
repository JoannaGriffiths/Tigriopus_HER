library(nlme)

setwd("~/Documents/Tigriopus_whole_genome_seq/Phys_data")

Data=read.delim("LT50_data")

attach(Data)
##########################################################################

hist(Data$LT50, main="Histogram")
shapiro.test(Data$LT50) #if p-value below 0.05 then data is not normally distributed
bartlett.test(Data$LT50, TREAT)
bartlett.test(Data$LT50, SEX)


model2 = lme(LT50 ~ TREAT*SEX, random = ~ 1|LINE, data = Data)
summary(model2)
aov_LT50 = aov(model2)
summary(aov_LT50)

library(ggplot2)
library(scales)

dat <- data.frame(
  Treatment = factor(c("S1", "S2", "S3", "S4", "S5", "C1", "C2", "C3", "C4", "C5", "BR", "SD"), levels=c("HS1", "HS2", "HS3", "HS4", "HS5", "CT1", "CT2", "CT3", "CT4", "CT5", "BR", "SD")),
  LT50 = c(37.44411, 37.26742, 36.80642, 36.70757, 36.84582, 35.61022, 35.2774, 35.95806, 35.84852, 36.0938, 35.22171, 37.308),   
  sterror = c(0.1026294, 0.1079124, 0.09447984, 0.0739102, 0.08145478, 0.08617707, 0.04119953, 0.07198648, 0.09950195, 0.08289513, 0.1071527, 0.08668574)    
)
dat

quartz()
ggplot(data=dat, aes(x=Treatment, y=LT50, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.8) +
  geom_errorbar(data=dat, aes(ymin=LT50-sterror, ymax=LT50+sterror), width=.2, position=position_dodge(.8)) +
  scale_fill_manual(values=c("#FF9900", "#FF9900", "#FF9900", "#FF9900", "#FF9900", "#3399FF", "#3399FF", "#3399FF", "#3399FF", "#3399FF", "#00CC66", "#FFFF66")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12)) +
  labs(y=expression("LT50 ("^o*"C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=12)) +
  scale_y_continuous(limits=c(34,38),oob = rescale_none)

##graph for females
dat <- data.frame(
  Treatment = factor(c("S1", "S2", "S3", "S4", "S5", "C1", "C2", "C3", "C4", "C5", "BR", "SD", "BR-SD", "SD-BR"), levels=c("S1", "S2", "S3", "S4", "S5", "C1", "C2", "C3", "C4", "C5", "BR", "SD", "BR-SD", "SD-BR")),
  groups = factor(c("S", "S", "S", "S", "S", "C", "C", "C", "C", "C", "P", "P", "F1", "F1"), levels = c("S", "C", "P", "F1")),
  LT50 = c(36.93656, 36.99241, 36.70432, 36.4117, 36.6736, 35.50809, 35.20537, 35.74767, 35.71135, 35.73916, 35.22171, 37.308, 36.59413, 36.48214),   
  sterror = c(0.09493073, 0.08696999, 0.1165887, 0.06719275, 0.07972128, 0.09062819, 0.107618, 0.07572416, 0.08200082, 0.1100629, 0.1071527, 0.08668574, 0.05062267, 0.0576367)    
)
dat

#graphs with 4 colors: "#FF9900", "#FF9900", "#FF9900", "#FF9900", "#FF9900", "#3399FF", "#3399FF", "#3399FF", "#3399FF", "#3399FF", "#00CC66", "#FFFF66", "#9933CC", "#9933CC"
#grey68", "grey68", "grey68", "grey68", "grey68", "grey33", "grey33", "grey33", "grey33", "grey33", "grey68", "grey68", "grey33", "grey33"
quartz()
ggplot(data=dat, aes(x=Treatment, y=LT50, fill=groups)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.8) +
  geom_errorbar(data=dat, aes(ymin=LT50-sterror, ymax=LT50+sterror), width=.2, position=position_dodge(.8)) +
  scale_fill_manual(values=c("grey68", "grey33", "grey68", "grey33")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12)) +
  labs(y=expression("LT50 ("^o*"C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.position = "none", 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=12)) +
  scale_y_continuous(limits=c(34,38),oob = rescale_none)


##graph for males
dat <- data.frame(
  Treatment = factor(c("S1", "S2", "S3", "S4", "S5", "C1", "C2", "C3", "C4", "C5", "BR", "SD", "BR-SD", "SD-BR"), levels=c("S1", "S2", "S3", "S4", "S5", "C1", "C2", "C3", "C4", "C5", "BR", "SD", "BR-SD", "SD-BR")),
  groups = factor(c("S", "S", "S", "S", "S", "C", "C", "C", "C", "C", "P", "P", "F1", "F1"), levels = c("S", "C", "P", "F1")),
  LT50 = c(37.44411, 37.26742, 36.80642, 36.70757, 36.84582, 35.61022, 35.2774, 35.95806, 35.84852, 36.0938, 35.22171, 37.308, 36.59413, 36.48214),   
  sterror = c(0.1026294, 0.1079124, 0.09447984, 0.0739102, 0.08145478, 0.08617707, 0.04119953, 0.07198648, 0.09950195, 0.08289513, 0.1071527, 0.08668574, 0.05062267, 0.0576367)    
)
dat

quartz()
ggplot(data=dat, aes(x=Treatment, y=LT50, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.8) +
  geom_errorbar(data=dat, aes(ymin=LT50-sterror, ymax=LT50+sterror), width=.2, position=position_dodge(.8)) +
  scale_fill_manual(values=c("#FF9900", "#FF9900", "#FF9900", "#FF9900", "#FF9900", "#3399FF", "#3399FF", "#3399FF", "#3399FF", "#3399FF", "#00CC66", "#FFFF66", "#9933CC", "#9933CC")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12)) +
  labs(y=expression("LT50 ("^o*"C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=12)) +
  scale_y_continuous(limits=c(34,38),oob = rescale_none)



