######This code has been adapted from the paper: Brennan RS, Garrett AD, Huber KE, Hargarten H, Pespeni MH, Brennan RS, Pespeni MH. 2019. Rare genetic variation and balanced polymorphisms are important for survival in global change conditions. Proc. R. Soc. B 286.
#https://github.com/PespeniLab/urchin_single_gen_selection/blob/master/scripts/Fig_03.R

library(scales)
library(plyr)
library(dplyr)
library(tidyr)


setwd("~/Documents/Tigriopus_whole_genome_seq/LD/LDx_results")
###############################read in 10,000bp window cmh results 
cmh <-read.table("../../cmh_files/cmh_10000window_overlap_forLD.txt", header=TRUE)
cmh <- cmh[c(2,3,4,7)]
cmh <- separate(cmh, FeatureID, c("BP_start", "BP_end"))

cmh$id_wind <- paste(cmh$Chromosome, cmh$BP_start, sep=".")
#create column with chromosome and bp merged together
cmh$id <- paste(cmh$Chromosome, cmh$Start, sep=".")

cmh_sig_id_list <- data.frame(c(1.10460000,1.11080000,1.11090000,1.11100000,1.11870000,1.11970000,1.11990000,1.12010000,1.13550000,1.13750000,1.14300000,1.14620000,1.14640000,1.14650000,1.14760000,1.290000,1.3380000,1.3560000,1.370000,1.3800000,1.390000,1.430000,1.4560000,1.4580000,1.4750000,1.620000,1.640000,1.7120000,1.7200000,1.7210000,1.7220000,1.7230000,1.7240000,1.7250000,1.7260000,1.7270000,1.7280000,1.7290000,
                                1.7350000,1.7360000,1.7370000,1.7380000,1.7390000,1.7400000,1.7410000,1.7420000,1.7450000,1.7510000,1.7520000,1.7550000,1.7580000,1.7690000,1.7830000,1.7920000,1.8150000,1.8190000,2.2090000,2.2170000,2.2320000,2.2440000,2.2570000,2.2580000,2.2600000,2.2650000,2.2690000,2.2700000,2.2720000,2.2730000,2.2740000,2.2750000,2.2780000,2.2790000,2.2820000,2.2830000,2.2860000,2.2880000,2.2890000,2.2920000,
                                2.3e+06,2.3030000,2.3070000,2.890000,3.10140000,3.10160000,3.10380000,3.10390000,3.10570000,3.14360000,3.14370000,3.14440000,3.14490000,3.14590000,3.14700000,3.2550000,3.6360000,4.1e+07,4.10010000,4.10050000,4.10060000,4.10070000,4.10080000,4.10090000,4.10110000,4.10120000,4.10130000,4.10140000,4.10170000,4.10220000,4.10350000,4.10360000,4.10370000,4.10390000,4.10410000,4.10430000,4.10440000,
                                4.10450000,4.10470000,4.10490000,4.10500000,4.10510000,4.10520000,4.10530000,4.10540000,4.10550000,4.10560000,4.10570000,4.10580000,4.10590000,4.10600000,4.10610000,4.10620000,4.10630000,4.10640000,4.10690000,4.10710000,4.10720000,4.10740000,4.10750000,4.10760000,4.10770000,4.10790000,4.10800000,4.11890000,4.11900000,4.11910000,4.11920000,4.11930000,4.11940000,4.11950000,4.11960000,4.12070000,
                                4.12080000,4.12120000,4.12130000,4.12160000,4.12170000,4.12180000,4.12210000,4.12230000,4.12260000,4.12270000,4.12300000,4.12310000,4.12320000,4.12330000,4.12340000,4.12350000, 4.12370000,4.12380000,4.12390000,4.12410000,4.12420000,4.12430000,4.12480000,4.12490000,4.12580000,4.12600000,4.12650000,4.13320000,4.9070000,4.9120000,4.9130000,4.9140000,4.9160000,4.9180000,4.9200000,4.9220000,
                                4.9290000,4.9310000,4.9330000,4.9360000,4.9380000,4.9400000,4.9460000,4.9490000,4.9560000,4.9590000,4.9630000,4.9640000,4.9650000,4.9660000,4.9670000,4.9680000,4.9690000,4.9710000,4.9720000,4.9730000,4.9750000,4.9800000,4.9810000,4.9830000,4.9850000,4.9860000,4.9870000,4.9880000,4.9910000,4.9920000,4.9940000,4.9970000,5.12220000,5.12440000,5.320000,5.390000,5.4290000,5.4680000,5.5e+06,5.5250000,
                                5.5280000,5.5610000,5.5780000,5.5840000,5.5900000,5.6150000,5.6250000,5.7060000,6.10640000,6.10820000,6.10890000,6.11130000,6.12170000,6.12310000,6.12420000,6.12480000,6.12560000,6.12590000,6.13820000,6.2010000,6.8250000,6.8260000,6.8270000,6.8280000,6.8290000,6.8360000,6.8380000,6.8440000,6.8460000,6.8510000,6.8520000,6.8700000,6.9050000,6.9630000,7.10080000,7.10140000,7.10490000,7.10500000,7.10780000,
                                7.10790000,7.10850000,7.10880000,7.10920000,7.10930000,7.10940000,7.10960000,7.10970000,7.10980000,7.10990000,7.11020000,7.11030000,7.11040000,7.11060000,7.11090000,7.11160000,7.11240000,7.11270000,7.11320000,7.11330000,7.11400000,7.11440000,7.11460000,7.11530000,7.11690000,7.11710000,7.11730000,7.11780000,7.11860000,7.11940000,7.11950000,7.12050000,7.12080000,7.12120000,7.12130000,7.12230000,
                                7.12590000,7.12660000,7.12980000,7.13830000,7.13860000,7.1580000,7.1710000,7.1810000,7.1880000,7.2090000,7.4890000,7.4990000,7.5180000,7.5220000,7.5490000,7.5500000,7.5510000,7.5530000,7.5540000,7.5550000,7.5560000,7.5580000,7.5590000,7.5600000,7.5630000,7.5670000,7.5750000,7.5840000,7.9330000,7.9770000,7.9810000,7.9840000,8.1e+07,8.10010000,8.10020000,8.10030000,8.10040000,8.10050000,8.10060000,
                                8.10080000,8.10090000,8.1010000,8.10100000,8.10110000,8.10130000,8.10140000,8.1040000,8.1160000,8.1190000,8.1200000,8.1250000,8.1290000,8.1300000,8.14420000,8.1480000,8.1580000,8.1810000,8.2470000,8.2910000,8.2920000,8.2930000,8.2940000,8.2960000,8.2970000,8.3090000,8.3160000,8.3230000,8.3280000,8.3340000,8.3350000,8.3360000, 8.3430000,8.3610000,8.3630000,8.3690000,8.3740000,8.3750000,8.3760000,
                                8.3780000,8.3800000,8.3810000,8.3820000,8.3840000,8.3860000,8.3880000,8.3900000,8.3910000,8.3920000,8.3940000,8.3950000,8.3980000,8.4e+06,8.4010000,8.4040000,8.4080000,8.4090000,8.4100000,8.4120000,8.4130000,8.4140000,8.4150000,8.4160000,8.4170000,8.4190000,8.4260000,8.4270000,8.4280000,8.4300000,8.4330000,8.4340000,8.4350000,8.4360000,8.4370000,8.4380000,8.4390000,8.4400000,8.4410000,8.4430000,
                                8.4440000,8.4450000,8.4460000,8.4470000,8.4480000,8.4490000,8.4500000,8.4510000,8.4530000,8.4550000,8.4590000,8.4640000,8.4710000,8.4750000,8.4760000,8.4780000,8.4790000,8.4800000,8.4810000,8.4820000,8.4840000,8.4850000,8.4860000,8.4900000,8.4910000,8.4940000,8.5e+06,8.6350000,8.6460000,8.6550000,8.6560000,8.6590000,8.6630000,8.6650000,8.6660000,8.6710000,8.6760000,8.6770000,8.6780000,
                                8.6800000,8.6830000,8.6840000,8.6850000,8.6860000,8.6870000,8.6920000,8.6960000,8.6980000,8.6990000,8.7e+05,8.7e+06,8.7020000,8.7030000,8.7060000,8.7070000,8.7080000,8.7090000,8.7450000,8.7480000,8.7570000,8.7670000,8.7790000,8.780000,8.8420000,8.8980000,8.920000,8.9330000,8.9340000,8.9360000,8.9380000,8.9390000,8.9410000,8.9420000,8.9440000,8.9490000,8.9500000,8.9510000,8.9550000,8.9560000,
                                8.9600000,8.9610000,8.9630000,8.9650000,8.9660000,8.9670000,8.9680000,8.9720000,8.9740000,8.9750000,8.9760000,8.9770000,8.9780000,8.9790000,8.980000,8.9800000,8.9810000,8.9830000,8.9850000,8.9860000,8.9890000,8.9900000,8.9910000,8.9940000,8.9950000,8.9960000,8.9970000,8.9980000,9.10040000,9.13180000,9.13730000,9.13760000,9.13780000,9.13920000,9.14060000,9.14700000,9.260000,9.4600000,9.4650000,
                                9.4670000,9.4680000,9.4690000,9.4700000,9.4710000,9.4730000,9.4750000,9.4810000,9.4830000,9.4850000,9.4890000,9.4900000,9.5090000,9.5100000,9.5150000,9.5170000,9.5180000,9.5200000,9.5230000,9.5330000,9.5340000,9.5350000,9.5380000,9.9770000,10.11340000,10.12050000,10.12210000,10.14800000,10.16110000,10.7490000,10.7540000,10.9530000,11.11620000,11.13330000,11.1730000,11.4230000,11.4690000,11.4700000,11.7420000,11.7830000,11.7960000,12.110000,12.1340000,12.15030000,12.16470000,12.1710000,12.17280000,12.3070000,12.3680000,12.8910000))

cmh_sig_id_list2 <- c(1.10460000,1.11080000,1.11090000,1.11100000,1.11870000,1.11970000,1.11990000,1.12010000,1.13550000,1.13750000,1.14300000,1.14620000,1.14640000,1.14650000,1.14760000,1.290000,1.3380000,1.3560000,1.370000,1.3800000,1.390000,1.430000,1.4560000,1.4580000,1.4750000,1.620000,1.640000,1.7120000,1.7200000,1.7210000,1.7220000,1.7230000,1.7240000,1.7250000,1.7260000,1.7270000,1.7280000,1.7290000,
                      1.7350000,1.7360000,1.7370000,1.7380000,1.7390000,1.7400000,1.7410000,1.7420000,1.7450000,1.7510000,1.7520000,1.7550000,1.7580000,1.7690000,1.7830000,1.7920000,1.8150000,1.8190000,2.2090000,2.2170000,2.2320000,2.2440000,2.2570000,2.2580000,2.2600000,2.2650000,2.2690000,2.2700000,2.2720000,2.2730000,2.2740000,2.2750000,2.2780000,2.2790000,2.2820000,2.2830000,2.2860000,2.2880000,2.2890000,2.2920000,
                      2.3e+06,2.3030000,2.3070000,2.890000,3.10140000,3.10160000,3.10380000,3.10390000,3.10570000,3.14360000,3.14370000,3.14440000,3.14490000,3.14590000,3.14700000,3.2550000,3.6360000,4.1e+07,4.10010000,4.10050000,4.10060000,4.10070000,4.10080000,4.10090000,4.10110000,4.10120000,4.10130000,4.10140000,4.10170000,4.10220000,4.10350000,4.10360000,4.10370000,4.10390000,4.10410000,4.10430000,4.10440000,
                      4.10450000,4.10470000,4.10490000,4.10500000,4.10510000,4.10520000,4.10530000,4.10540000,4.10550000,4.10560000,4.10570000,4.10580000,4.10590000,4.10600000,4.10610000,4.10620000,4.10630000,4.10640000,4.10690000,4.10710000,4.10720000,4.10740000,4.10750000,4.10760000,4.10770000,4.10790000,4.10800000,4.11890000,4.11900000,4.11910000,4.11920000,4.11930000,4.11940000,4.11950000,4.11960000,4.12070000,
                      4.12080000,4.12120000,4.12130000,4.12160000,4.12170000,4.12180000,4.12210000,4.12230000,4.12260000,4.12270000,4.12300000,4.12310000,4.12320000,4.12330000,4.12340000,4.12350000, 4.12370000,4.12380000,4.12390000,4.12410000,4.12420000,4.12430000,4.12480000,4.12490000,4.12580000,4.12600000,4.12650000,4.13320000,4.9070000,4.9120000,4.9130000,4.9140000,4.9160000,4.9180000,4.9200000,4.9220000,
                      4.9290000,4.9310000,4.9330000,4.9360000,4.9380000,4.9400000,4.9460000,4.9490000,4.9560000,4.9590000,4.9630000,4.9640000,4.9650000,4.9660000,4.9670000,4.9680000,4.9690000,4.9710000,4.9720000,4.9730000,4.9750000,4.9800000,4.9810000,4.9830000,4.9850000,4.9860000,4.9870000,4.9880000,4.9910000,4.9920000,4.9940000,4.9970000,5.12220000,5.12440000,5.320000,5.390000,5.4290000,5.4680000,5.5e+06,5.5250000,
                      5.5280000,5.5610000,5.5780000,5.5840000,5.5900000,5.6150000,5.6250000,5.7060000,6.10640000,6.10820000,6.10890000,6.11130000,6.12170000,6.12310000,6.12420000,6.12480000,6.12560000,6.12590000,6.13820000,6.2010000,6.8250000,6.8260000,6.8270000,6.8280000,6.8290000,6.8360000,6.8380000,6.8440000,6.8460000,6.8510000,6.8520000,6.8700000,6.9050000,6.9630000,7.10080000,7.10140000,7.10490000,7.10500000,7.10780000,
                      7.10790000,7.10850000,7.10880000,7.10920000,7.10930000,7.10940000,7.10960000,7.10970000,7.10980000,7.10990000,7.11020000,7.11030000,7.11040000,7.11060000,7.11090000,7.11160000,7.11240000,7.11270000,7.11320000,7.11330000,7.11400000,7.11440000,7.11460000,7.11530000,7.11690000,7.11710000,7.11730000,7.11780000,7.11860000,7.11940000,7.11950000,7.12050000,7.12080000,7.12120000,7.12130000,7.12230000,
                      7.12590000,7.12660000,7.12980000,7.13830000,7.13860000,7.1580000,7.1710000,7.1810000,7.1880000,7.2090000,7.4890000,7.4990000,7.5180000,7.5220000,7.5490000,7.5500000,7.5510000,7.5530000,7.5540000,7.5550000,7.5560000,7.5580000,7.5590000,7.5600000,7.5630000,7.5670000,7.5750000,7.5840000,7.9330000,7.9770000,7.9810000,7.9840000,8.1e+07,8.10010000,8.10020000,8.10030000,8.10040000,8.10050000,8.10060000,
                      8.10080000,8.10090000,8.1010000,8.10100000,8.10110000,8.10130000,8.10140000,8.1040000,8.1160000,8.1190000,8.1200000,8.1250000,8.1290000,8.1300000,8.14420000,8.1480000,8.1580000,8.1810000,8.2470000,8.2910000,8.2920000,8.2930000,8.2940000,8.2960000,8.2970000,8.3090000,8.3160000,8.3230000,8.3280000,8.3340000,8.3350000,8.3360000, 8.3430000,8.3610000,8.3630000,8.3690000,8.3740000,8.3750000,8.3760000,
                      8.3780000,8.3800000,8.3810000,8.3820000,8.3840000,8.3860000,8.3880000,8.3900000,8.3910000,8.3920000,8.3940000,8.3950000,8.3980000,8.4e+06,8.4010000,8.4040000,8.4080000,8.4090000,8.4100000,8.4120000,8.4130000,8.4140000,8.4150000,8.4160000,8.4170000,8.4190000,8.4260000,8.4270000,8.4280000,8.4300000,8.4330000,8.4340000,8.4350000,8.4360000,8.4370000,8.4380000,8.4390000,8.4400000,8.4410000,8.4430000,
                      8.4440000,8.4450000,8.4460000,8.4470000,8.4480000,8.4490000,8.4500000,8.4510000,8.4530000,8.4550000,8.4590000,8.4640000,8.4710000,8.4750000,8.4760000,8.4780000,8.4790000,8.4800000,8.4810000,8.4820000,8.4840000,8.4850000,8.4860000,8.4900000,8.4910000,8.4940000,8.5e+06,8.6350000,8.6460000,8.6550000,8.6560000,8.6590000,8.6630000,8.6650000,8.6660000,8.6710000,8.6760000,8.6770000,8.6780000,
                      8.6800000,8.6830000,8.6840000,8.6850000,8.6860000,8.6870000,8.6920000,8.6960000,8.6980000,8.6990000,8.7e+05,8.7e+06,8.7020000,8.7030000,8.7060000,8.7070000,8.7080000,8.7090000,8.7450000,8.7480000,8.7570000,8.7670000,8.7790000,8.780000,8.8420000,8.8980000,8.920000,8.9330000,8.9340000,8.9360000,8.9380000,8.9390000,8.9410000,8.9420000,8.9440000,8.9490000,8.9500000,8.9510000,8.9550000,8.9560000,
                      8.9600000,8.9610000,8.9630000,8.9650000,8.9660000,8.9670000,8.9680000,8.9720000,8.9740000,8.9750000,8.9760000,8.9770000,8.9780000,8.9790000,8.980000,8.9800000,8.9810000,8.9830000,8.9850000,8.9860000,8.9890000,8.9900000,8.9910000,8.9940000,8.9950000,8.9960000,8.9970000,8.9980000,9.10040000,9.13180000,9.13730000,9.13760000,9.13780000,9.13920000,9.14060000,9.14700000,9.260000,9.4600000,9.4650000,
                      9.4670000,9.4680000,9.4690000,9.4700000,9.4710000,9.4730000,9.4750000,9.4810000,9.4830000,9.4850000,9.4890000,9.4900000,9.5090000,9.5100000,9.5150000,9.5170000,9.5180000,9.5200000,9.5230000,9.5330000,9.5340000,9.5350000,9.5380000,9.9770000,10.11340000,10.12050000,10.12210000,10.14800000,10.16110000,10.7490000,10.7540000,10.9530000,11.11620000,11.13330000,11.1730000,11.4230000,11.4690000,
                      11.4700000,11.7420000,11.7830000,11.7960000,12.110000,12.1340000,12.15030000,12.16470000,12.1710000,12.17280000,12.3070000,12.3680000,12.8910000)


colnames(cmh_sig_id_list) <- c("id_wind")

cmh$id_wind <- as.numeric(cmh$id_wind)

cmh_top <- merge(cmh, cmh_sig_id_list, by="id_wind")

cmh_neut <- cmh[!cmh$id_wind %in% cmh_sig_id_list2,]



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
