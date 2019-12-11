
library(scales)
library(plyr)
library(dplyr)
library(coin)


setwd("~/Documents/Tigriopus_whole_genome_seq/Seq_divergence")
kaks = read.delim("BR_kaks.txt")
kaks$Feature <- kaks$GeneID
kaks$not_sig <- rep("all",nrow(kaks))

mean(kaks$ka_div_ks_plus1) #genomewide mean is 0.66744

siggenes = read.delim("../cmh_files/cmh_sigwind_geneoverlap.txt") 
siggenes$sig <- rep("overlap",nrow(siggenes))

sig_kaks_overlap = merge(kaks, siggenes, by="Feature")
kaks_sig_table <- data.frame(sig_kaks_overlap$Feature, sig_kaks_overlap$ka_div_ks_plus1)

library("dplyr")
kaks_sig_table2 <- distinct(kaks_sig_table)
mean(kaks_sig_table$sig_kaks_overlap.ka_div_ks_plus1) #mean of significant genes is 0.4787611

write.csv(kaks_sig_table2, file = "kaks_cmhsig")

############### making file for permutation

sig_kaks_all = merge(kaks, siggenes, by="Feature", all.x=TRUE)
kaks_all_table <- data.frame(sig_kaks_all$Feature, sig_kaks_all$ka_div_ks_plus1, sig_kaks_all$not_sig, sig_kaks_all$sig)

write.csv(kaks_all_table, file = "kaks_cmhsig_all")

#couldnt figure out how to replace NAs so I just opened it up in excel and add "all" for places that are NA

kaks_all = read.delim("kaks_cmhsig_all_fixed")



#t test
mod <- t.test(kaks_all$sig_kaks_all.ka_div_ks_plus1 ~ kaks_all$sig_kaks_all.sig)
mod #t = 7.8406, df = 563.31, p-value = 2.256e-14

#permutations
perm_length <- length(which(kaks_all$sig_kaks_all.sig=="overlap")) #for test dataset:sel_kaks$sel
#define or clear list
perm.out2 <- list()



for (i in 1:1000){
  rand.sel <- rep(0, nrow(kaks_all))
  rand.sel[sample(1:nrow(kaks_all), perm_length)]<- 1
  selected.perm <- data.frame( kaks_all$sig_kaks_all.Feature, kaks_all$sig_kaks_all.ka_div_ks_plus1, kaks_all$sig_kaks_all.sig, rand.sel )
  mod <- t.test(kaks_all$sig_kaks_all.ka_div_ks_plus1 ~ rand.sel)
  perm.out2[[(length(perm.out2)+1)]] <- mod$p.value 
}



perm <- do.call(rbind.data.frame, perm.out2)
colnames(perm) <- c("p")



quantile(perm, 0.025, na.rm=TRUE)
#1000 perm: 0.007343216
quantile(perm, 0.975, na.rm=TRUE)
#1000 perm: 0.9815728 

percentile <- ecdf(perm$p)
percentile(2.256e-14)
#my pvalue is in the 0th percentile (very significant!)



