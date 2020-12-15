# Scripts and data from Griffiths, J.S. et al. (2020) Differential responses to ocean acidification between populations of *Balanophyllia elegans* corals from high and low upwelling environments. Mol. Ecol. DOI: 10.1111/mec.15050

## Gene Expression Analysis

*Script*: DESeq2.R

*Input files*: Tig_RSEM_merged_matrix, column_Tig_CT.txt, column_Tig_HS.txt, column_Tig_all.txt, column_Tig_selected.txt, column_Tig_unselected.txt

*Output files*: Tig_output_DEG_CT_sel_vs_unsel.txt, Tig_output_DEG_HS_sel_vs_unsel.txt, Tig_output_DEG_all_HSvsCT.txt, Tig_output_DEG_all_Sel_HS.txt, Tig_output_DEG_all_sel_vs_unsel.txt, Tig_output_DEG_sel_HS_vs_CT.txt, Tig_output_DEG_unsel_HS_vs_CT.txt

*Description*: Script contains code for performing differential gene expression analysis among treatment conditions. The merged matrix contains all the gene counts for the genes. The different column text files contain treatment information associated with each sample name.



## KaKs Analysis

*Script*: kaks.py

*Input files*: 

*Output files*: BR_kaks.txt

*Description*: kaks.py is a python script that converts annotated vcf file (for the Bodega population) into a tab delimited file with KaKs values for each gene. The output from this script is the BR_kaks.txt file contains kaks values for the Bodega population compared to the San Diego reference genome.


*Script*: KaKs_permutation.R

*Input files*: BR_kaks.txt, cmh_sigwind_geneoverlap.txt

*Output files*: kaks_cmhsig_all_fixed

*Description*: Script runs a permutation and t-test to determine if the kaks values for genes under heat tolerance selection have lower kaks values than the genome-wide average.



## LD Analysis

*Script*: DeSeq2.R

*Input files*: Orthoblast_RSEM_merged_matrix, column_2transcriptomes.txt

*Output files*: output_DEG_GOLtp1.txt, output_DEG_GOLtp2.txt, output_DEG_PACtp1.txt, output_DEG_PACtp2.txt

*Description*: Script contains analysis using DeSeq2. The merged matrix file contains the total counts of contains for each “Orthogroup”. The column file gives an explanation for the headers in the matrix file (the treatment conditions and population names for each read file. Output files contain the pvalues for log fold changes in Orthogroups in response to pH for each population timepoint (day 9 and 29).



## LT50 Analysis

*Script*: Please see scripts and explanations located at: https://github.com/z0on/GO_MWU\

*Input files*: Orthoblast_interproresults_nonredun.csv, pvalue_GOL_tp1, pvalue_GOL_tp2, pvalue_PAC_tp1, pvalue_PAC_tp2

*Output files*: BP_pvalue_GOL_tp1, BP_pvalue_GOL_tp2, BP_pvalue_PAC_tp1, BP_pvalue_PAC_tp2, MF_pvalue_GOL_tp1, MF_pvalue_GOL_tp2, MF_pvalue_PAC_tp1, MF_pvalue_PAC_tp2

*Description*: Script contains analysis for functional enrichment of logfold contig changes. The Orthoblast_interproresults_nonredun.csv input file contains the GO terms associated with each Orthogroup. The pvalue input files contain the signed logfold pvalues for expression changes derived from the DEG analysis using DeSeq2. Output files correspind to each input file name and whether the Biological Processes (BP) or Molecular Functions (MF) were analyzed.


