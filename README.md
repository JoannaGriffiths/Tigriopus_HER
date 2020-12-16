# Scripts and data from Joanna S Griffiths, Yasmeen Kawji, Morgan W Kelly. 2020. An Experimental Test of Adaptive Introgression in Locally Adapted Populations of Splash Pool Copepods, Molecular Biology and Evolution, https://doi.org/10.1093/molbev/msaa289

## Gene Expression Analysis

*Script*: DESeq2.R

*Input files*: Tig_RSEM_merged_matrix, column_Tig_CT.txt, column_Tig_HS.txt, column_Tig_all.txt, column_Tig_selected.txt, column_Tig_unselected.txt

*Output files*: Tig_output_DEG_CT_sel_vs_unsel.txt, Tig_output_DEG_HS_sel_vs_unsel.txt, Tig_output_DEG_all_HSvsCT.txt, Tig_output_DEG_all_Sel_HS.txt, Tig_output_DEG_all_sel_vs_unsel.txt, Tig_output_DEG_sel_HS_vs_CT.txt, Tig_output_DEG_unsel_HS_vs_CT.txt

*Description*: Script contains code for performing differential gene expression analysis among treatment conditions. The merged matrix contains all the gene counts for the genes. The different column text files contain treatment information associated with each sample name.



## KaKs Analysis

*Script*: kaks.py

*Input files*: Bodega annotated vcf file

*Output files*: BR_kaks.txt

*Description*: kaks.py is a python script that converts annotated vcf file (for the Bodega population) into a tab delimited file with KaKs values for each gene. The output from this script is the BR_kaks.txt file contains kaks values for the Bodega population compared to the San Diego reference genome.


*Script*: KaKs_permutation.R

*Input files*: BR_kaks.txt, cmh_sigwind_geneoverlap.txt

*Output files*: kaks_cmhsig_all_fixed

*Description*: Script runs a permutation and t-test to determine if the kaks values for genes under heat tolerance selection have a significantly different mean kaks value than the genome-wide average.



## LD Analysis
Note: LD analysis not included in final manuscript

*Script*: mergeChr_LD_results.R

*Input files*: 1S_LDx_chr1, 1S_LDx_chr10, 1S_LDx_chr11, 1S_LDx_chr2, 1S_LDx_chr3, 1S_LDx_chr4, 1S_LDx_chr5, 1S_LDx_chr6, 1S_LDx_chr7, 1S_LDx_chr8, 1S_LDx_chr9

*Output files*: LD_1S

*Description*: Script contains takes output from LDx.pl and merges all chromosome level results into a single file for each sample. Files above are an example for one of the lines under selection (1S). The following was repeated for all selected and control lines.


*Script*: LD_vs_Distance.R

*Input files*: cmh_10000window_overlap_forLD.tx, LD_1S, LD_2S, LD_3S, LD_4S, LD_5S, LD_1U, LD_2U, LD_4U, LD_5U, LD_6U

*Output*: figure displaying LD vs bp distance between SNPs

*Description*: Script models LD as a function of distance between SNPs and displays results as a graph. Hypothesis was that SNPs under heat tolerance selection had a higher LD than SNPs not under heat tolerance selection in selected lines. We also hypothesized SNPs under heat tolerance selection in selected lines had higher LD than these same SNPs in the control line.



## LT50 Analysis

*Script*: LT50_ANOVA.R

*Input files*: LT50_data

*Description*: Script runs ANOVA to compare LT50 values among treatment groups.


