library(scales)
library(plyr)
library(dplyr)

setwd("~/Documents/Tigriopus_whole_genome_seq/LD/LDx_reuslts")

########################### 1S

LD_1S_chr1=read.delim("1S_LDx_chr1")

colnames(LD_1S_chr1) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")

#add chr column
LD_1S_chr1$Chr <- rep(1,nrow(LD_1S_chr1))

LD_1S_chr2=read.delim("1S_LDx_chr2")
colnames(LD_1S_chr2) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr2$Chr <- rep(2,nrow(LD_1S_chr2))

LD_1S_chr3=read.delim("1S_LDx_chr3")
colnames(LD_1S_chr3) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr3$Chr <- rep(3,nrow(LD_1S_chr3))

LD_1S_chr4=read.delim("1S_LDx_chr4")
colnames(LD_1S_chr4) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr4$Chr <- rep(4,nrow(LD_1S_chr4))

LD_1S_chr5=read.delim("1S_LDx_chr5")
colnames(LD_1S_chr5) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr5$Chr <- rep(5,nrow(LD_1S_chr5))

LD_1S_chr6=read.delim("1S_LDx_chr6")
colnames(LD_1S_chr6) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr6$Chr <- rep(6,nrow(LD_1S_chr6))

LD_1S_chr7=read.delim("1S_LDx_chr7")
colnames(LD_1S_chr7) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr7$Chr <- rep(7,nrow(LD_1S_chr7))

LD_1S_chr8=read.delim("1S_LDx_chr8")
colnames(LD_1S_chr8) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr8$Chr <- rep(8,nrow(LD_1S_chr8))

LD_1S_chr9=read.delim("1S_LDx_chr9")
colnames(LD_1S_chr9) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr9$Chr <- rep(9,nrow(LD_1S_chr9))

LD_1S_chr10=read.delim("1S_LDx_chr10")
colnames(LD_1S_chr10) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr10$Chr <- rep(10,nrow(LD_1S_chr10))

LD_1S_chr11=read.delim("1S_LDx_chr11")
colnames(LD_1S_chr11) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr11$Chr <- rep(11,nrow(LD_1S_chr11))

LD_1S_chr12=read.delim("1S_LDx_chr12")
colnames(LD_1S_chr12) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1S_chr12$Chr <- rep(12,nrow(LD_1S_chr12))

#merge dataframs vertically 
LD_1S <- rbind(LD_1S_chr1, LD_1S_chr2, LD_1S_chr4, LD_1S_chr4, LD_1S_chr5, LD_1S_chr6, LD_1S_chr7, LD_1S_chr8, LD_1S_chr9, LD_1S_chr10, LD_1S_chr11)
write.csv(LD_1S, file = "LD_1S")


########################### 1U

LD_1U_chr1=read.delim("1U_LDx_chr1")

colnames(LD_1U_chr1) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")

#add chr column
LD_1U_chr1$Chr <- rep(1,nrow(LD_1U_chr1))

LD_1U_chr2=read.delim("1U_LDx_chr2")
colnames(LD_1U_chr2) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr2$Chr <- rep(2,nrow(LD_1U_chr2))

LD_1U_chr3=read.delim("1U_LDx_chr3")
colnames(LD_1U_chr3) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr3$Chr <- rep(3,nrow(LD_1U_chr3))

LD_1U_chr4=read.delim("1U_LDx_chr4")
colnames(LD_1U_chr4) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr4$Chr <- rep(4,nrow(LD_1U_chr4))

LD_1U_chr5=read.delim("1U_LDx_chr5")
colnames(LD_1U_chr5) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr5$Chr <- rep(5,nrow(LD_1U_chr5))

LD_1U_chr6=read.delim("1U_LDx_chr6")
colnames(LD_1U_chr6) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr6$Chr <- rep(6,nrow(LD_1U_chr6))

LD_1U_chr7=read.delim("1U_LDx_chr7")
colnames(LD_1U_chr7) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr7$Chr <- rep(7,nrow(LD_1U_chr7))

LD_1U_chr8=read.delim("1U_LDx_chr8")
colnames(LD_1U_chr8) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr8$Chr <- rep(8,nrow(LD_1U_chr8))

LD_1U_chr9=read.delim("1U_LDx_chr9")
colnames(LD_1U_chr9) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr9$Chr <- rep(9,nrow(LD_1U_chr9))

LD_1U_chr10=read.delim("1U_LDx_chr10")
colnames(LD_1U_chr10) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr10$Chr <- rep(10,nrow(LD_1U_chr10))

LD_1U_chr11=read.delim("1U_LDx_chr11")
colnames(LD_1U_chr11) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr11$Chr <- rep(11,nrow(LD_1U_chr11))

LD_1U_chr12=read.delim("1U_LDx_chr12")
colnames(LD_1U_chr12) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_1U_chr12$Chr <- rep(12,nrow(LD_1U_chr12))

#merge dataframs vertically 
LD_1U <- rbind(LD_1U_chr1, LD_1U_chr2, LD_1U_chr4, LD_1U_chr4, LD_1U_chr5, LD_1U_chr6, LD_1U_chr7, LD_1U_chr8, LD_1U_chr9, LD_1U_chr10, LD_1U_chr11, LD_1U_chr12)
write.csv(LD_1U, file = "LD_1U")


########################### 2S

LD_2S_chr1=read.delim("2S_LDx_chr1")

colnames(LD_2S_chr1) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")

#add chr column
LD_2S_chr1$Chr <- rep(1,nrow(LD_2S_chr1))

LD_2S_chr2=read.delim("2S_LDx_chr2")
colnames(LD_2S_chr2) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr2$Chr <- rep(2,nrow(LD_2S_chr2))

LD_2S_chr3=read.delim("2S_LDx_chr3")
colnames(LD_2S_chr3) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr3$Chr <- rep(3,nrow(LD_2S_chr3))

LD_2S_chr4=read.delim("2S_LDx_chr4")
colnames(LD_2S_chr4) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr4$Chr <- rep(4,nrow(LD_2S_chr4))

LD_2S_chr5=read.delim("2S_LDx_chr5")
colnames(LD_2S_chr5) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr5$Chr <- rep(5,nrow(LD_2S_chr5))

LD_2S_chr6=read.delim("2S_LDx_chr6")
colnames(LD_2S_chr6) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr6$Chr <- rep(6,nrow(LD_2S_chr6))

LD_2S_chr7=read.delim("2S_LDx_chr7")
colnames(LD_2S_chr7) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr7$Chr <- rep(7,nrow(LD_2S_chr7))

LD_2S_chr8=read.delim("2S_LDx_chr8")
colnames(LD_2S_chr8) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr8$Chr <- rep(8,nrow(LD_2S_chr8))

LD_2S_chr9=read.delim("2S_LDx_chr9")
colnames(LD_2S_chr9) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr9$Chr <- rep(9,nrow(LD_2S_chr9))

LD_2S_chr10=read.delim("2S_LDx_chr10")
colnames(LD_2S_chr10) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr10$Chr <- rep(10,nrow(LD_2S_chr10))

LD_2S_chr11=read.delim("2S_LDx_chr11")
colnames(LD_2S_chr11) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr11$Chr <- rep(11,nrow(LD_2S_chr11))

LD_2S_chr12=read.delim("2S_LDx_chr12")
colnames(LD_2S_chr12) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2S_chr12$Chr <- rep(12,nrow(LD_2S_chr12))

#merge dataframs vertically 
LD_2S <- rbind(LD_2S_chr1, LD_2S_chr2, LD_2S_chr4, LD_2S_chr4, LD_2S_chr5, LD_2S_chr6, LD_2S_chr7, LD_2S_chr8, LD_2S_chr9, LD_2S_chr10, LD_2S_chr11, LD_2S_chr12)
write.csv(LD_2S, file = "LD_2S")


########################### 2U

LD_2U_chr1=read.delim("2U_LDx_chr1")

colnames(LD_2U_chr1) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")

#add chr column
LD_2U_chr1$Chr <- rep(1,nrow(LD_2U_chr1))

LD_2U_chr2=read.delim("2U_LDx_chr2")
colnames(LD_2U_chr2) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr2$Chr <- rep(2,nrow(LD_2U_chr2))

LD_2U_chr3=read.delim("2U_LDx_chr3")
colnames(LD_2U_chr3) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr3$Chr <- rep(3,nrow(LD_2U_chr3))

LD_2U_chr4=read.delim("2U_LDx_chr4")
colnames(LD_2U_chr4) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr4$Chr <- rep(4,nrow(LD_2U_chr4))

LD_2U_chr5=read.delim("2U_LDx_chr5")
colnames(LD_2U_chr5) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr5$Chr <- rep(5,nrow(LD_2U_chr5))

LD_2U_chr6=read.delim("2U_LDx_chr6")
colnames(LD_2U_chr6) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr6$Chr <- rep(6,nrow(LD_2U_chr6))

LD_2U_chr7=read.delim("2U_LDx_chr7")
colnames(LD_2U_chr7) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr7$Chr <- rep(7,nrow(LD_2U_chr7))

LD_2U_chr8=read.delim("2U_LDx_chr8")
colnames(LD_2U_chr8) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr8$Chr <- rep(8,nrow(LD_2U_chr8))

LD_2U_chr9=read.delim("2U_LDx_chr9")
colnames(LD_2U_chr9) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr9$Chr <- rep(9,nrow(LD_2U_chr9))

LD_2U_chr10=read.delim("2U_LDx_chr10")
colnames(LD_2U_chr10) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr10$Chr <- rep(10,nrow(LD_2U_chr10))

LD_2U_chr11=read.delim("2U_LDx_chr11")
colnames(LD_2U_chr11) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr11$Chr <- rep(11,nrow(LD_2U_chr11))

LD_2U_chr12=read.delim("2U_LDx_chr12")
colnames(LD_2U_chr12) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_2U_chr12$Chr <- rep(12,nrow(LD_2U_chr12))

#merge dataframs vertically 
LD_2U <- rbind(LD_2U_chr1, LD_2U_chr2, LD_2U_chr4, LD_2U_chr4, LD_2U_chr5, LD_2U_chr6, LD_2U_chr7, LD_2U_chr8, LD_2U_chr9, LD_2U_chr10, LD_2U_chr11, LD_2U_chr12)
write.csv(LD_2U, file = "LD_2U")


########################### 3S

LD_3S_chr1=read.delim("3S_LDx_chr1")

colnames(LD_3S_chr1) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")

#add chr column
LD_3S_chr1$Chr <- rep(1,nrow(LD_3S_chr1))

LD_3S_chr2=read.delim("3S_LDx_chr2")
colnames(LD_3S_chr2) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr2$Chr <- rep(2,nrow(LD_3S_chr2))

LD_3S_chr3=read.delim("3S_LDx_chr3")
colnames(LD_3S_chr3) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr3$Chr <- rep(3,nrow(LD_3S_chr3))

LD_3S_chr4=read.delim("3S_LDx_chr4")
colnames(LD_3S_chr4) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr4$Chr <- rep(4,nrow(LD_3S_chr4))

LD_3S_chr5=read.delim("3S_LDx_chr5")
colnames(LD_3S_chr5) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr5$Chr <- rep(5,nrow(LD_3S_chr5))

LD_3S_chr6=read.delim("3S_LDx_chr6")
colnames(LD_3S_chr6) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr6$Chr <- rep(6,nrow(LD_3S_chr6))

LD_3S_chr7=read.delim("3S_LDx_chr7")
colnames(LD_3S_chr7) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr7$Chr <- rep(7,nrow(LD_3S_chr7))

LD_3S_chr8=read.delim("3S_LDx_chr8")
colnames(LD_3S_chr8) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr8$Chr <- rep(8,nrow(LD_3S_chr8))

LD_3S_chr9=read.delim("3S_LDx_chr9")
colnames(LD_3S_chr9) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr9$Chr <- rep(9,nrow(LD_3S_chr9))

LD_3S_chr10=read.delim("3S_LDx_chr10")
colnames(LD_3S_chr10) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr10$Chr <- rep(10,nrow(LD_3S_chr10))

LD_3S_chr11=read.delim("3S_LDx_chr11")
colnames(LD_3S_chr11) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr11$Chr <- rep(11,nrow(LD_3S_chr11))

LD_3S_chr12=read.delim("3S_LDx_chr12")
colnames(LD_3S_chr12) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_3S_chr12$Chr <- rep(12,nrow(LD_3S_chr12))

#merge dataframs vertically 
LD_3S <- rbind(LD_3S_chr1, LD_3S_chr2, LD_3S_chr4, LD_3S_chr4, LD_3S_chr5, LD_3S_chr6, LD_3S_chr7, LD_3S_chr8, LD_3S_chr9, LD_3S_chr10, LD_3S_chr11, LD_3S_chr12)
write.csv(LD_3S, file = "LD_3S")


########################### 4S

LD_4S_chr1=read.delim("4S_LDx_chr1")

colnames(LD_4S_chr1) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")

#add chr column
LD_4S_chr1$Chr <- rep(1,nrow(LD_4S_chr1))

LD_4S_chr2=read.delim("4S_LDx_chr2")
colnames(LD_4S_chr2) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr2$Chr <- rep(2,nrow(LD_4S_chr2))

LD_4S_chr3=read.delim("4S_LDx_chr3")
colnames(LD_4S_chr3) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr3$Chr <- rep(3,nrow(LD_4S_chr3))

LD_4S_chr4=read.delim("4S_LDx_chr4")
colnames(LD_4S_chr4) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr4$Chr <- rep(4,nrow(LD_4S_chr4))

LD_4S_chr5=read.delim("4S_LDx_chr5")
colnames(LD_4S_chr5) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr5$Chr <- rep(5,nrow(LD_4S_chr5))

LD_4S_chr6=read.delim("4S_LDx_chr6")
colnames(LD_4S_chr6) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr6$Chr <- rep(6,nrow(LD_4S_chr6))

LD_4S_chr7=read.delim("4S_LDx_chr7")
colnames(LD_4S_chr7) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr7$Chr <- rep(7,nrow(LD_4S_chr7))

LD_4S_chr8=read.delim("4S_LDx_chr8")
colnames(LD_4S_chr8) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr8$Chr <- rep(8,nrow(LD_4S_chr8))

LD_4S_chr9=read.delim("4S_LDx_chr9")
colnames(LD_4S_chr9) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr9$Chr <- rep(9,nrow(LD_4S_chr9))

LD_4S_chr10=read.delim("4S_LDx_chr10")
colnames(LD_4S_chr10) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr10$Chr <- rep(10,nrow(LD_4S_chr10))

LD_4S_chr11=read.delim("4S_LDx_chr11")
colnames(LD_4S_chr11) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr11$Chr <- rep(11,nrow(LD_4S_chr11))

LD_4S_chr12=read.delim("4S_LDx_chr12")
colnames(LD_4S_chr12) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4S_chr12$Chr <- rep(12,nrow(LD_4S_chr12))

#merge dataframs vertically 
LD_4S <- rbind(LD_4S_chr1, LD_4S_chr2, LD_4S_chr4, LD_4S_chr4, LD_4S_chr5, LD_4S_chr6, LD_4S_chr7, LD_4S_chr8, LD_4S_chr9, LD_4S_chr10, LD_4S_chr11, LD_4S_chr12)
write.csv(LD_4S, file = "LD_4S")


########################### 4U

LD_4U_chr1=read.delim("4U_LDx_chr1")

colnames(LD_4U_chr1) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")

#add chr column
LD_4U_chr1$Chr <- rep(1,nrow(LD_4U_chr1))

LD_4U_chr2=read.delim("4U_LDx_chr2")
colnames(LD_4U_chr2) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr2$Chr <- rep(2,nrow(LD_4U_chr2))

LD_4U_chr3=read.delim("4U_LDx_chr3")
colnames(LD_4U_chr3) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr3$Chr <- rep(3,nrow(LD_4U_chr3))

LD_4U_chr4=read.delim("4U_LDx_chr4")
colnames(LD_4U_chr4) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr4$Chr <- rep(4,nrow(LD_4U_chr4))

LD_4U_chr5=read.delim("4U_LDx_chr5")
colnames(LD_4U_chr5) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr5$Chr <- rep(5,nrow(LD_4U_chr5))

LD_4U_chr6=read.delim("4U_LDx_chr6")
colnames(LD_4U_chr6) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr6$Chr <- rep(6,nrow(LD_4U_chr6))

LD_4U_chr7=read.delim("4U_LDx_chr7")
colnames(LD_4U_chr7) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr7$Chr <- rep(7,nrow(LD_4U_chr7))

LD_4U_chr8=read.delim("4U_LDx_chr8")
colnames(LD_4U_chr8) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr8$Chr <- rep(8,nrow(LD_4U_chr8))

LD_4U_chr9=read.delim("4U_LDx_chr9")
colnames(LD_4U_chr9) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr9$Chr <- rep(9,nrow(LD_4U_chr9))

LD_4U_chr10=read.delim("4U_LDx_chr10")
colnames(LD_4U_chr10) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr10$Chr <- rep(10,nrow(LD_4U_chr10))

LD_4U_chr11=read.delim("4U_LDx_chr11")
colnames(LD_4U_chr11) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr11$Chr <- rep(11,nrow(LD_4U_chr11))

LD_4U_chr12=read.delim("4U_LDx_chr12")
colnames(LD_4U_chr12) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_4U_chr12$Chr <- rep(12,nrow(LD_4U_chr12))

#merge dataframs vertically 
LD_4U <- rbind(LD_4U_chr1, LD_4U_chr2, LD_4U_chr4, LD_4U_chr4, LD_4U_chr5, LD_4U_chr6, LD_4U_chr7, LD_4U_chr8, LD_4U_chr9, LD_4U_chr10, LD_4U_chr11, LD_4U_chr12)
write.csv(LD_4U, file = "LD_4U")


########################### 5S

LD_5S_chr1=read.delim("5S_LDx_chr1")

colnames(LD_5S_chr1) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")

#add chr column
LD_5S_chr1$Chr <- rep(1,nrow(LD_5S_chr1))

LD_5S_chr2=read.delim("5S_LDx_chr2")
colnames(LD_5S_chr2) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr2$Chr <- rep(2,nrow(LD_5S_chr2))

LD_5S_chr3=read.delim("5S_LDx_chr3")
colnames(LD_5S_chr3) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr3$Chr <- rep(3,nrow(LD_5S_chr3))

LD_5S_chr4=read.delim("5S_LDx_chr4")
colnames(LD_5S_chr4) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr4$Chr <- rep(4,nrow(LD_5S_chr4))

LD_5S_chr5=read.delim("5S_LDx_chr5")
colnames(LD_5S_chr5) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr5$Chr <- rep(5,nrow(LD_5S_chr5))

LD_5S_chr6=read.delim("5S_LDx_chr6")
colnames(LD_5S_chr6) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr6$Chr <- rep(6,nrow(LD_5S_chr6))

LD_5S_chr7=read.delim("5S_LDx_chr7")
colnames(LD_5S_chr7) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr7$Chr <- rep(7,nrow(LD_5S_chr7))

LD_5S_chr8=read.delim("5S_LDx_chr8")
colnames(LD_5S_chr8) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr8$Chr <- rep(8,nrow(LD_5S_chr8))

LD_5S_chr9=read.delim("5S_LDx_chr9")
colnames(LD_5S_chr9) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr9$Chr <- rep(9,nrow(LD_5S_chr9))

LD_5S_chr10=read.delim("5S_LDx_chr10")
colnames(LD_5S_chr10) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr10$Chr <- rep(10,nrow(LD_5S_chr10))

LD_5S_chr11=read.delim("5S_LDx_chr11")
colnames(LD_5S_chr11) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr11$Chr <- rep(11,nrow(LD_5S_chr11))

LD_5S_chr12=read.delim("5S_LDx_chr12")
colnames(LD_5S_chr12) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5S_chr12$Chr <- rep(12,nrow(LD_5S_chr12))

#merge dataframs vertically 
LD_5S <- rbind(LD_5S_chr1, LD_5S_chr2, LD_5S_chr4, LD_5S_chr4, LD_5S_chr5, LD_5S_chr6, LD_5S_chr7, LD_5S_chr8, LD_5S_chr9, LD_5S_chr10, LD_5S_chr11, LD_5S_chr12)
write.csv(LD_5S, file = "LD_5S")


########################### 5U

LD_5U_chr1=read.delim("5U_LDx_chr1")

colnames(LD_5U_chr1) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")

#add chr column
LD_5U_chr1$Chr <- rep(1,nrow(LD_5U_chr1))

LD_5U_chr2=read.delim("5U_LDx_chr2")
colnames(LD_5U_chr2) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr2$Chr <- rep(2,nrow(LD_5U_chr2))

LD_5U_chr3=read.delim("5U_LDx_chr3")
colnames(LD_5U_chr3) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr3$Chr <- rep(3,nrow(LD_5U_chr3))

LD_5U_chr4=read.delim("5U_LDx_chr4")
colnames(LD_5U_chr4) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr4$Chr <- rep(4,nrow(LD_5U_chr4))

LD_5U_chr5=read.delim("5U_LDx_chr5")
colnames(LD_5U_chr5) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr5$Chr <- rep(5,nrow(LD_5U_chr5))

LD_5U_chr6=read.delim("5U_LDx_chr6")
colnames(LD_5U_chr6) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr6$Chr <- rep(6,nrow(LD_5U_chr6))

LD_5U_chr7=read.delim("5U_LDx_chr7")
colnames(LD_5U_chr7) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr7$Chr <- rep(7,nrow(LD_5U_chr7))

LD_5U_chr8=read.delim("5U_LDx_chr8")
colnames(LD_5U_chr8) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr8$Chr <- rep(8,nrow(LD_5U_chr8))

LD_5U_chr9=read.delim("5U_LDx_chr9")
colnames(LD_5U_chr9) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr9$Chr <- rep(9,nrow(LD_5U_chr9))

LD_5U_chr10=read.delim("5U_LDx_chr10")
colnames(LD_5U_chr10) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr10$Chr <- rep(10,nrow(LD_5U_chr10))

LD_5U_chr11=read.delim("5U_LDx_chr11")
colnames(LD_5U_chr11) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr11$Chr <- rep(11,nrow(LD_5U_chr11))

LD_5U_chr12=read.delim("5U_LDx_chr12")
colnames(LD_5U_chr12) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_5U_chr12$Chr <- rep(12,nrow(LD_5U_chr12))

#merge dataframs vertically 
LD_5U <- rbind(LD_5U_chr1, LD_5U_chr2, LD_5U_chr4, LD_5U_chr4, LD_5U_chr5, LD_5U_chr6, LD_5U_chr7, LD_5U_chr8, LD_5U_chr9, LD_5U_chr10, LD_5U_chr11, LD_5U_chr12)
write.csv(LD_5U, file = "LD_5U")


########################### 6U

LD_6U_chr1=read.delim("6U_LDx_chr1")

colnames(LD_6U_chr1) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")

#add chr column
LD_6U_chr1$Chr <- rep(1,nrow(LD_6U_chr1))

LD_6U_chr2=read.delim("6U_LDx_chr2")
colnames(LD_6U_chr2) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr2$Chr <- rep(2,nrow(LD_6U_chr2))

LD_6U_chr3=read.delim("6U_LDx_chr3")
colnames(LD_6U_chr3) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr3$Chr <- rep(3,nrow(LD_6U_chr3))

LD_6U_chr4=read.delim("6U_LDx_chr4")
colnames(LD_6U_chr4) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr4$Chr <- rep(4,nrow(LD_6U_chr4))

LD_6U_chr5=read.delim("6U_LDx_chr5")
colnames(LD_6U_chr5) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr5$Chr <- rep(5,nrow(LD_6U_chr5))

LD_6U_chr6=read.delim("6U_LDx_chr6")
colnames(LD_6U_chr6) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr6$Chr <- rep(6,nrow(LD_6U_chr6))

LD_6U_chr7=read.delim("6U_LDx_chr7")
colnames(LD_6U_chr7) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr7$Chr <- rep(7,nrow(LD_6U_chr7))

LD_6U_chr8=read.delim("6U_LDx_chr8")
colnames(LD_6U_chr8) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr8$Chr <- rep(8,nrow(LD_6U_chr8))

LD_6U_chr9=read.delim("6U_LDx_chr9")
colnames(LD_6U_chr9) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr9$Chr <- rep(9,nrow(LD_6U_chr9))

LD_6U_chr10=read.delim("6U_LDx_chr10")
colnames(LD_6U_chr10) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr10$Chr <- rep(10,nrow(LD_6U_chr10))

LD_6U_chr11=read.delim("6U_LDx_chr11")
colnames(LD_6U_chr11) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr11$Chr <- rep(11,nrow(LD_6U_chr11))

LD_6U_chr12=read.delim("6U_LDx_chr12")
colnames(LD_6U_chr12) <- c("snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
LD_6U_chr12$Chr <- rep(12,nrow(LD_6U_chr12))

#merge dataframs vertically 
LD_6U <- rbind(LD_6U_chr1, LD_6U_chr2, LD_6U_chr4, LD_6U_chr4, LD_6U_chr5, LD_6U_chr6, LD_6U_chr7, LD_6U_chr8, LD_6U_chr9, LD_6U_chr10, LD_6U_chr11, LD_6U_chr12)
write.csv(LD_6U, file = "LD_6U")



