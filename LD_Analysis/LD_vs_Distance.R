######This code has been adapted from the paper: Brennan RS, Garrett AD, Huber KE, Hargarten H, Pespeni MH, Brennan RS, Pespeni MH. 2019. Rare genetic variation and balanced polymorphisms are important for survival in global change conditions. Proc. R. Soc. B 286.
#https://github.com/PespeniLab/urchin_single_gen_selection/blob/master/scripts/Fig_03.R

library(scales)
library(plyr)
library(dplyr)
library(tidyr)
library("DataCombine")


setwd("~/Documents/Tigriopus_whole_genome_seq/LD/LDx_results")
setwd("~/LSU/Research/Tigriopus_selection/LD_analysis")
###############################read in 10,000bp window cmh results 
cmh <-read.table("10000window_SD_cmh_overlap.txt", header=TRUE)
cmh <- cmh[c(2,3,6)]
cmh <- separate(cmh, FeatureID, c("BP_start", "BP_end"))
Replaces <- data.frame(from=c("Chromosome_"), to=c(""))
cmh <- FindReplace(data = cmh, Var = "Chromosome", replaceData = Replaces,
                    from = "from", to = "to", exact = FALSE)


cmh$id_wind <- paste(cmh$Chromosome, cmh$BP_start, sep=".")
#create column with chromosome and bp merged together
cmh$id <- paste(cmh$Chromosome, cmh$Start, sep=".")

cmh_sig <- read.delim("cmh_window_geomean_sigonly_SD_round2")
cmh_sig_id_list <- cmh_sig[c("ID")]
colnames(cmh_sig_id_list) <- c("id_wind")

colnames(cmh_sig_id_list) <- c("id_wind")
cmh_neut <- cmh[!cmh$id_wind %in% cmh_sig_id_list$id_wind,]

cmh_top <- merge(cmh_sig_id_list, cmh, by="id_wind")


#read in LD results for selected replicate 1S
LD_1S=read.csv("LD_1S", header = T)
#drop column 1
LD_1S <- LD_1S[-c(1)]
#add column with replicate number
LD_1S$Rep <- rep(1,nrow(LD_1S))
#create column with chromosome and bp merged together
LD_1S$id1 <- paste(LD_1S$Chr, LD_1S$snp1, sep=".")
LD_1S$id2 <- paste(LD_1S$Chr, LD_1S$snp2, sep=".")

#read in replicate 2S
LD_2S=read.csv("LD_2S", header = T)
LD_2S <- LD_2S[-c(1)]
LD_2S$Rep <- rep(2,nrow(LD_2S))
LD_2S$id1 <- paste(LD_2S$Chr, LD_2S$snp1, sep=".")
LD_2S$id2 <- paste(LD_2S$Chr, LD_2S$snp2, sep=".")

#read in replicate 3S
LD_3S=read.csv("LD_3S", header = T)
LD_3S <- LD_3S[-c(1)]
LD_3S$Rep <- rep(3,nrow(LD_3S))
LD_3S$id1 <- paste(LD_3S$Chr, LD_3S$snp1, sep=".")
LD_3S$id2 <- paste(LD_3S$Chr, LD_3S$snp2, sep=".")

#read in replicate 4S
LD_4S=read.csv("LD_4S", header = T)
LD_4S <- LD_4S[-c(1)]
LD_4S$Rep <- rep(4,nrow(LD_4S))
LD_4S$id1 <- paste(LD_4S$Chr, LD_4S$snp1, sep=".")
LD_4S$id2 <- paste(LD_4S$Chr, LD_4S$snp2, sep=".")

#read in replicate 5S
LD_5S=read.csv("LD_5S", header = T)
LD_5S <- LD_5S[-c(1)]
LD_5S$Rep <- rep(5,nrow(LD_5S))
LD_5S$id1 <- paste(LD_5S$Chr, LD_5S$snp1, sep=".")
LD_5S$id2 <- paste(LD_5S$Chr, LD_5S$snp2, sep=".")

#read in replicate 1U
LD_1U=read.csv("LD_1U", header = T)
LD_1U <- LD_1U[-c(1)]
LD_1U$Rep <- rep(1,nrow(LD_1U))
LD_1U$id1 <- paste(LD_1U$Chr, LD_1U$snp1, sep=".")
LD_1U$id2 <- paste(LD_1U$Chr, LD_1U$snp2, sep=".")

#read in replicate 2U
LD_2U=read.csv("LD_2U", header = T)
LD_2U <- LD_2U[-c(1)]
LD_2U$Rep <- rep(2,nrow(LD_2U))
LD_2U$id1 <- paste(LD_2U$Chr, LD_2U$snp1, sep=".")
LD_2U$id2 <- paste(LD_2U$Chr, LD_2U$snp2, sep=".")

#read in replicate 4U
LD_4U=read.csv("LD_1U", header = T)
LD_4U <- LD_4U[-c(1)]
LD_4U$Rep <- rep(3,nrow(LD_4U))
LD_4U$id1 <- paste(LD_4U$Chr, LD_1U$snp1, sep=".")
LD_4U$id2 <- paste(LD_4U$Chr, LD_1U$snp2, sep=".")

#read in replicate 5U
LD_5U=read.csv("LD_5U", header = T)
LD_5U <- LD_5U[-c(1)]
LD_5U$Rep <- rep(4,nrow(LD_5U))
LD_5U$id1 <- paste(LD_5U$Chr, LD_5U$snp1, sep=".")
LD_5U$id2 <- paste(LD_5U$Chr, LD_5U$snp2, sep=".")

#read in replicate 6U
LD_6U=read.csv("LD_6U", header = T)
LD_6U <- LD_6U[-c(1)]
LD_6U$Rep <- rep(5,nrow(LD_6U))
LD_6U$id1 <- paste(LD_6U$Chr, LD_6U$snp1, sep=".")
LD_6U$id2 <- paste(LD_6U$Chr, LD_6U$snp2, sep=".")



#merge LD results and cmh results by chr:bp column first for selected sites then neutral sites
df_sel1_1S <- merge(LD_1S, cmh_top, by.x="id1", by.y="id")
df_sel2_1S <- merge(LD_1S, cmh_top, by.x="id2", by.y="id")
df_sel_1S <- rbind(df_sel1_1S, df_sel2_1S)

df_neut1_1S <- merge(LD_1S, cmh_neut, by.x="id1", by.y="id")
df_neut2_1S <- merge(LD_1S, cmh_neut, by.x="id2", by.y="id")
df_neut_1S <- rbind(df_neut1_1S, df_neut2_1S)


df_sel1_2S <- merge(LD_2S, cmh_top, by.x="id1", by.y="id")
df_sel2_2S <- merge(LD_2S, cmh_top, by.x="id2", by.y="id")
df_sel_2S <- rbind(df_sel1_2S, df_sel2_2S)

df_neut1_2S <- merge(LD_2S, cmh_neut, by.x="id1", by.y="id")
df_neut2_2S <- merge(LD_2S, cmh_neut, by.x="id2", by.y="id")
df_neut_2S <- rbind(df_neut1_2S, df_neut2_2S)


df_sel1_3S <- merge(LD_3S, cmh_top, by.x="id1", by.y="id")
df_sel2_3S <- merge(LD_3S, cmh_top, by.x="id2", by.y="id")
df_sel_3S <- rbind(df_sel1_3S, df_sel2_3S)

df_neut1_3S <- merge(LD_3S, cmh_neut, by.x="id1", by.y="id")
df_neut2_3S <- merge(LD_3S, cmh_neut, by.x="id2", by.y="id")
df_neut_3S <- rbind(df_neut1_3S, df_neut2_3S)


df_sel1_4S <- merge(LD_4S, cmh_top, by.x="id1", by.y="id")
df_sel2_4S <- merge(LD_4S, cmh_top, by.x="id2", by.y="id")
df_sel_4S <- rbind(df_sel1_4S, df_sel2_4S)

df_neut1_4S <- merge(LD_4S, cmh_neut, by.x="id1", by.y="id")
df_neut2_4S <- merge(LD_4S, cmh_neut, by.x="id2", by.y="id")
df_neut_4S <- rbind(df_neut1_4S, df_neut2_4S)


df_sel1_5S <- merge(LD_5S, cmh_top, by.x="id1", by.y="id")
df_sel2_5S <- merge(LD_5S, cmh_top, by.x="id2", by.y="id")
df_sel_5S <- rbind(df_sel1_5S, df_sel2_5S)

df_neut1_5S <- merge(LD_5S, cmh_neut, by.x="id1", by.y="id")
df_neut2_5S <- merge(LD_5S, cmh_neut, by.x="id2", by.y="id")
df_neut_5S <- rbind(df_neut1_5S, df_neut2_5S)

df_neut1_1U <- merge(LD_1U, cmh_top, by.x="id1", by.y="id")
df_neut2_1U <- merge(LD_1U, cmh_top, by.x="id2", by.y="id")
df_neut_1U <- rbind(df_neut1_1U, df_neut2_1U)

df_neut1_2U <- merge(LD_2U, cmh_top, by.x="id1", by.y="id")
df_neut2_2U <- merge(LD_2U, cmh_top, by.x="id2", by.y="id")
df_neut_2U <- rbind(df_neut1_2U, df_neut2_2U)

df_neut1_4U <- merge(LD_4U, cmh_top, by.x="id1", by.y="id")
df_neut2_4U <- merge(LD_4U, cmh_top, by.x="id2", by.y="id")
df_neut_4U <- rbind(df_neut1_4U, df_neut2_4U)

df_neut1_5U <- merge(LD_5U, cmh_top, by.x="id1", by.y="id")
df_neut2_5U <- merge(LD_5U, cmh_top, by.x="id2", by.y="id")
df_neut_5U <- rbind(df_neut1_5U, df_neut2_5U)

df_neut1_6U <- merge(LD_6U, cmh_top, by.x="id1", by.y="id")
df_neut2_6U <- merge(LD_6U, cmh_top, by.x="id2", by.y="id")
df_neut_6U <- rbind(df_neut1_6U, df_neut2_6U)

#calculate distance between SNPs from LD results and make new column
df_sel_1S$distance <- abs(df_sel_1S$snp1 - df_sel_1S$snp2)
df_neut_1S$distance <- abs(df_neut_1S$snp1 - df_neut_1S$snp2)

df_sel_2S$distance <- abs(df_sel_2S$snp1 - df_sel_2S$snp2)
df_neut_2S$distance <- abs(df_neut_2S$snp1 - df_neut_2S$snp2)

df_sel_3S$distance <- abs(df_sel_3S$snp1 - df_sel_3S$snp2)
df_neut_3S$distance <- abs(df_neut_3S$snp1 - df_neut_3S$snp2)

df_sel_4S$distance <- abs(df_sel_4S$snp1 - df_sel_4S$snp2)
df_neut_4S$distance <- abs(df_neut_4S$snp1 - df_neut_4S$snp2)

df_sel_5S$distance <- abs(df_sel_5S$snp1 - df_sel_5S$snp2)
df_neut_5S$distance <- abs(df_neut_5S$snp1 - df_neut_5S$snp2)


df_neut_1U$distance <- abs(df_neut_1U$snp1 - df_neut_1U$snp2)
df_neut_2U$distance <- abs(df_neut_2U$snp1 - df_neut_2U$snp2)
df_neut_4U$distance <- abs(df_neut_4U$snp1 - df_neut_4U$snp2)
df_neut_5U$distance <- abs(df_neut_5U$snp1 - df_neut_5U$snp2)
df_neut_6U$distance <- abs(df_neut_6U$snp1 - df_neut_6U$snp2)



#merge all replicate files for selected and neutral SNPs separately (vertically)
df_sel_all <- rbind(df_sel_1S, df_sel_2S, df_sel_3S, df_sel_4S, df_sel_5S)
df_neut_all <- rbind(df_neut_1S, df_neut_2S, df_neut_3S, df_neut_4S, df_neut_5S)
df_unsel_neut_all <- rbind(df_neut_1U, df_neut_2U, df_neut_4U, df_neut_5U, df_neut_6U)

#df_neut_all <- rbind(df_neut_2S, df_neut_3S, df_neut_4S, df_neut_5S)

#perform stats/model
library(nlme)
m.sel.1S <- lm(mle_est ~ log10(distance), data=df_sel_1S)
fitted(m.sel.1S)
anova(m.sel.1S)
summary(m.sel.1S)
# intercept: 0.65521; distance: -0.08294

m.neut.1S <- lm(mle_est ~ log10(distance), data=df_neut_1S)
anova(m.neut.1S)
summary(m.neut.1S)
# intercept: 0.569853; distance: -0.103162


m.sel.2S <- lm(mle_est ~ log10(distance), data=df_sel_2S)
anova(m.sel.2S)
summary(m.sel.2S)
# intercept: 0.807265; distance: -0.134851

m.neut.2S <- lm(mle_est ~ log10(distance), data=df_neut_2S)
anova(m.neut.2S)
summary(m.neut.2S)
# intercept: 0.778915; distance: -0.103162


m.sel.3S <- lm(mle_est ~ log10(distance), data=df_sel_3S)
anova(m.sel.3S)
summary(m.sel.3S)
# intercept: 0.811825; distance: -0.108724

m.neut.3S <- lm(mle_est ~ log10(distance), data=df_neut_3S)
anova(m.neut.3S)
summary(m.neut.3S)
# intercept: 0.781190; distance: -0.107251


m.sel.4S <- lm(mle_est ~ log10(distance), data=df_sel_4S)
anova(m.sel.4S)
summary(m.sel.4S)
# intercept: 0.788265; distance: -0.119399

m.neut.4S <- lm(mle_est ~ log10(distance), data=df_neut_4S)
anova(m.neut.4S)
summary(m.neut.4S)
# intercept: 0.789024; distance: -0.117879


m.sel.5S <- lm(mle_est ~ log10(distance), data=df_sel_5S)
anova(m.sel.5S)
summary(m.sel.5S)
# intercept: 0.774598; distance: -0.164261

m.neut.5S <- lm(mle_est ~ log10(distance), data=df_neut_5S)
anova(m.neut.5S)
summary(m.neut.5S)
# intercept: 0.794923; distance: -0.152991


m.sel.all <- lme(mle_est ~ log10(distance), random=~1|Rep, data=df_sel_all)
#m.sel.all <- lme(mle_est ~ distance, random=~1|Rep, data=df_sel_all)
anova(m.sel.all)
summary(m.sel.all)
# intercept: 0.7762289; distance (slope): -0.122

m.neut.all <- lme(mle_est ~ log10(distance), random=~1|Rep, data=df_neut_all)
anova(m.neut.all)
summary(m.neut.all)
# intercept: 0.7401241; distance (slope): -0.120

m.unsel.neut.all <- lme(mle_est ~ log10(distance), random=~1|Rep, data=df_unsel_neut_all)
anova(m.unsel.neut.all)
summary(m.unsel.neut.all)



############

#figure out 95% confidence intervals for neutral SNPs 

############

df_sel_neut_all <- rbind(df_sel_all, df_neut_all)
perm_length <- length(which(df_sel_all == TRUE))

perm.out2 <- list()
perm.raw2 <- list()

for (i in 1:500){
  #first to distinguish between our neutral and selected SNPs--first call everything false (neutral)
  rand.sel <- rep("FALSE", nrow(cmh))
  #now identify the random "selected" values by calling them TRUE. The number selected is the same as the length of the dataset specificied above
  rand.sel[sample(1:nrow(cmh), perm_length)]<- TRUE
  selected.perm <- data.frame( cmh$Chromosome, cmh$Start, cmh$id, rand.sel )
  colnames(selected.perm) <- c("chr", "pos", "id", "selected")
  #merge whole LD results with the "selected" cmh results, all.x: if TRUE, then extra rows will be added to the output, one for each row in x that has no matching row in y. These rows will have NAs in those columns that are usually filled with values from y. The default is FALSE, so that only rows with data from both x and y are included in the output.
  dat.perm <- merge(df_sel_neut_all, selected.perm, by.x="id1", by.y="id", all.x=TRUE)
  
  #dat.perm3 <- merge(df_sel_neut_all, selected.perm, by.x="id2", by.y="id", all.x=TRUE)
  #this makes sure we get the rows back that didn't merge, so all SNPs remain in the dataset (almost all?)
  dat.perm3 <- merge(dat.perm, selected.perm, by.x="id2", by.y="id", all.x=TRUE)
  
  dat.perm3$distance <- abs(dat.perm3$snp1-dat.perm3$snp2)
  #keep rows where distance is greater than 100bp
  dat.perm3 <- dat.perm3[which(dat.perm3$distance <=100),]
  #keep rows where either column selected.x or selected.y has a value of True (ie not False)
  dat.perm4 <- dat.perm3[which(dat.perm3$selected.y == TRUE | dat.perm3$selected.x == TRUE),]
  mod.perm <- lm(mle_est ~ log10(distance), data = dat.perm4)
  #if dividing by 25, the remainder is 0, then print which loop its on
  if (i%%25 == 0){print(paste("current i loop",i))} # printing progress
  perm.out2[[(length(perm.out2)+1)]] <- cbind(sort(dat.perm4$distance, decreasing=FALSE), sort(fitted(mod.perm), decreasing=TRUE))
  perm.raw2[[(length(perm.out2)+1)]] <- dat.perm4
}

# for each permutation, take the mean of the bin, then combine all perms.

out <- list()
for (i in 1: length(perm.out2)){
  
  tmp.df <- data.frame(x=perm.out2[[i]][,1],
                       y= perm.out2[[i]][,2],
                       bin = cut(
                         perm.out2[[i]][,1],
                         breaks=seq(from=0, to=100, by=1),
                         labels=seq(from=1, to=100, by=1)))
  out[[i]] <- data.frame(bin=seq(from=1, to=100, by=1), ld=tapply(tmp.df$y,tmp.df$bin , mean))
}

out.new <- do.call(rbind, out)

densities.qtiles <- out.new %>%
  group_by(bin) %>%
  summarise(q05 = quantile(ld, 0.025, na.rm=TRUE),
            q50 = quantile(ld, 0.5, na.rm=TRUE),
            q95 = quantile(ld, 0.975, na.rm=TRUE))

ld.qtiles <- densities.qtiles

#compare out.new to selected and non-selected data
## estimate intercept and slope of all perms

intercept <- c()
slope <- c()

for (i in 2:length(perm.raw2)){
  m1 <- lm(mle_est ~ log10(distance),data=perm.raw2[[i]])
  intercept[i] <- summary(m1)$coefficients[1,1]
  slope[i] <- summary(m1)$coefficients[2,1]
  
}

quantile(intercept, c(0.025, 0.975), na.rm=TRUE)
quantile(slope, c(0.025, 0.975), na.rm=TRUE)



#> quantile(intercept, c(0.025, 0.975), na.rm=TRUE)
#     2.5%     97.5%
#0.7569020 0.7975452
#> quantile(slope, c(0.025, 0.975), na.rm=TRUE)
#       2.5%       97.5%
#-0.1261756 -0.1018579


#save(out.new, perm.out2, perm.raw2, m.sel.all, m.neut.all, m.sel.1S, m.sel.2S, m.sel.3S, m.sel.4S, 
#     m.sel.5S, m.neut.1S, m.neut.2S, m.neut.3S, m.neut.4S, m.neut.5S, densities.qtiles, file = "LDx_selSNP_top10.RData")

#save(out.new, perm.out2, perm.raw2, m.sel.all, m.neut.all, m.sel.1S, m.sel.2S, m.sel.3S, m.sel.4S, 
#    m.sel.5S, m.neut.1S, m.neut.2S, m.neut.3S, m.neut.4S, m.neut.5S, densities.qtiles, file = "LDx_selSNP_top0.1.RData")

#save(out.new, perm.out2, perm.raw2, m.sel.all, m.sel.1S, m.sel.2S, m.sel.3S, m.sel.4S, 
#     m.sel.5S, m.neut.1S, m.neut.2S, m.neut.3S, m.neut.4S, m.neut.5S, densities.qtiles, file = "LDx_selSNP_top1.RData")

#save(df_sel_1S, df_sel_2S, df_sel_3S, df_sel_4S, df_sel_5S, df_neut_1S, df_neut_2S, df_neut_3S, df_neut_4S, df_neut_5S, df_sel_all, df_neut_all, file = "LDx_selSNPs_top1_data2.RData")

save(df_sel_all, df_neut_all, m.unsel.neut.all, out.new, perm.out2, perm.raw2, m.sel.all, m.neut.all, m.sel.1S, m.sel.2S, m.sel.3S, m.sel.4S, 
     m.sel.5S, m.neut.1S, m.neut.2S, m.neut.3S, m.neut.4S, m.neut.5S, densities.qtiles, file = "LDx_selSNP_cmhwind.RData")


#load("LDx_selSNP_top10.RData") 
#load("LDx_selSNP_top0.1.RData") 
#load("LDx_selSNP_top1.RData") 

load("LDx_selSNP_cmhwind.RData")

#load("LDx_selSNPs_top1_data2.RData")


#plot
quartz()
plot(0,type='n', xlim=c(1,100), ylim=c(0.5,0.85),
     main="",
     ylab="",
     xlab="",
     cex.lab=0.9, cex.axis=0.7,
     xaxt="n",yaxt="n")

#plot 95% distribution
polygon(x=c(densities.qtiles$bin,rev(densities.qtiles$bin)),
        y=c(densities.qtiles$q05,rev(densities.qtiles$q95)),
        col=alpha("black", alpha=0.3),border=NA)

#lines(densities.qtiles$bin, densities.qtiles$q50, lwd=1.5, lty=1, col='black')

#single models
lines(sort(df_sel_1S$distance, decreasing=FALSE), sort(fitted(m.sel.1S), decreasing=TRUE),
      lwd=1.5, lty=3, col='firebrick3')
lines(sort(df_neut_1S$distance, decreasing=FALSE), sort(fitted(m.neut.1S), decreasing=TRUE),
      lwd=1.5, lty=3, col='royalblue3')
lines(sort(df_sel_2S$distance, decreasing=FALSE), sort(fitted(m.sel.2S), decreasing=TRUE),
      lwd=1.5, lty=3, col='firebrick3')
lines(sort(df_neut_2S$distance, decreasing=FALSE), sort(fitted(m.neut.2S), decreasing=TRUE),
      lwd=1.5, lty=3, col='royalblue3')
lines(sort(df_sel_3S$distance, decreasing=FALSE), sort(fitted(m.sel.3S), decreasing=TRUE),
      lwd=1.5, lty=3, col='firebrick3')
lines(sort(df_neut_3S$distance, decreasing=FALSE), sort(fitted(m.neut.3S), decreasing=TRUE),
      lwd=1.5, lty=3, col='royalblue3')
lines(sort(df_sel_4S$distance, decreasing=FALSE), sort(fitted(m.sel.4S), decreasing=TRUE),
      lwd=1.5, lty=3, col='firebrick3')
lines(sort(df_neut_4S$distance, decreasing=FALSE), sort(fitted(m.neut.4S), decreasing=TRUE),
      lwd=1.5, lty=3, col='royalblue3')
lines(sort(df_sel_5S$distance, decreasing=FALSE), sort(fitted(m.sel.5S), decreasing=TRUE),
      lwd=1.5, lty=3, col='firebrick3')
lines(sort(df_neut_5S$distance, decreasing=FALSE), sort(fitted(m.neut.5S), decreasing=TRUE),
      lwd=1.5, lty=3, col='royalblue3')

#full models
lines(smooth.spline(sort(df_sel_all$distance, decreasing=FALSE), sort(fitted(m.sel.all), decreasing=TRUE)),
      lwd=1.5, lty=1, col='firebrick3')
lines(smooth.spline(sort(df_neut_all$distance, decreasing=FALSE), sort(fitted(m.neut.all), decreasing=TRUE)),
      lwd=1.5, lty=1, col='royalblue3')
lines(smooth.spline(sort(df_unsel_neut_all$distance, decreasing=FALSE), sort(fitted(m.unsel.neut.all), decreasing=TRUE)),
      lwd=1.5, lty=3, col='royalblue3')

#legend
legend("bottomleft",
       legend=c("Selected SNPs Selected Line",
                "Neutral SNPs Selected Line", "Selected SNPs Control Line"),
       col=c("firebrick3","royalblue3", "royalblue3"),
       lty=c(1,1,3), lwd=2.2, cex=0.75, bty = "n")

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7, tcl=-0.2, at=c(1, 20, 40, 60, 80,100)) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)
title(xlab="Distance between SNPs in base pairs", line=1.5, cex.lab=0.9)
title(ylab="Estimated Linkage Disequilibrium", line=1.5, cex.lab=0.9)

dev.off()
