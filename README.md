## Handy R functions for genetics research 

Originally hosted at http://evachan.org/rscripts.html, these R functions were primarily written for my own research. However, (I'd like to think) I've written them to be generic enough to be useful to others. 

Please drop me an email if you spot a bug or have ideas for improvement.

[Statistical Functions] (../blob/master/README.md#stats)
[Plotting Functions] (../blob/master/README.md#plots)  
[Example Data] (https://github.com/ekfchan/evachan.org-Rscripts/blob/master/README.md#examples)  

###Statistical Functions

[geno_to_allelecnt.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/geno_to_allelecnt.R)  
A function to convert biallelic unphased SNP genotypes, such as {AA,CC,GG,TT,AC,AG,AT,CG,CT,GT}, to number of copies/counts {0,1,2} of the reference (or arbitrary) allele.  
[See example II and simgeno.R for example and usage.]  
[simgeno.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/simgeno.R)  
Very simple function to generate a biallelic unphased SNP genotype matrix in the format {AA,CC,GG,TT,AC,AG,AT,CG,CT,GT}. Used predominantly to test geno_to_allelecnt.R.  
[See example II for usage and purpose.]  
[calc_EHHS.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_EHHS.R)  
A function to calculate the normalised homozygosity between the i-th and j-th loci, EHHS(geno)i,j, for a given chromosome / linkage group (Tang, Thornton, Stoneking 2007)  
[calc_iES.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_iES.R)  
A function to calculate the integrated EHHS statistic, iES, as described in Tang, Thornton and Stoneking (2007).  You'd probably want to calculate the EHHS first!   
[calc_LD.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_LD.R)  
Given a biallelic genotype matrix, calculates one or more measures of linkage disequilibrium between all locus-pairs. The available LD measures include: D, D', r2, X2 (chi-square), X2' (chi-square-prime).  
[See example I for example data and usage.]  
[calc_snp_stats.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_snp_stats.R)  
A function to calculate basic SNP stats, including: allele frequency (p), MAF (minor allele frequency), MGF (minor genotype frequency), and tests for deviation from HWE (X2 test and Fisher's Exact test).  
[See example I for example data and usage.]  
[gwas_lm.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/gwas_lm.R)  
Performs single-locus (SNP) genome-wide association tests for one or more traits simultaneously under one or more of five inheritance models (additive, co-dominance, dominance, recessive, over-dominance) using linear regression.  
[calc_hwe_fisher.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_hwe_fisher.R)  
A script to test for deviation from HWE using Fisher's Exact test. This test is also incorporated into calc_snp_stats.R.  
[See example I for example data and usage.]  
[calc_hwe_chisq.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_hwe_chisq.R)  
A script to test for deviation from HWE using Pearson's Chi-Squared test. This test is also incorporated into calc_snp_stats.R.  
[See example I for example data and usage.]  
[calc_neiFis_multispop.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_neiFis_multispop.R)  
A script to calculate inbreeding coefficients, Fis, for each sub-population using a given set of SNP markers.  
[See example I for example data and usage.]  
[calc_neiFis_onepop.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_neiFis_onepop.R)  
A script to calculate inbreeding coefficients, Fis,  for a given population using a given set of SNP markers.  
[See example I for example data and usage.]  
[calc_wcFstats.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_wcFstats.R)  
A script to estimate the variance components and fixation indices as described in  Weir & Cockerham 1984 Evolution 38(6) : 1358-1370.  
[See example I for example data and usage. ]  
[calc_wcFst_spop_pairs.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_wcFst_spop_pairs.R)  
A script to estimate Fst (theta) values for each pair of sub-populations using the method of Weir & Cockerham 1984 Evolution 38(6): 1358-1370.  
[See example I for example data and usage.]  
[calc_allele_sharing.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/calc_allele_sharing.R)  
Calculates allele sharing distances between pairs of individuals (c.f. Gao & Stramer 2007 BMC Genetics 8:34).  
[See example I for example data and usage.]  


###Plotting Functions###

[plclust_in_colour.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/plclust_in_colour.R)  
A modification of (wrapper to) plclust for plotting hclust (hierarchical cluster) objects with coloured leaf labels.    
[plot_marker_lox.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/plot_marker_lox.R)  
Generates a visual representation of the genetic positions of a set of markers.  
[plot_markers_by_set.R](https://github.com/ekfchan/evachan.org-Rscripts/blob/master/rscripts/plot_markers_by_set.R)  
A function to plot sets of markers on a map where the markers are coloured based on a defined variable.   


###Example Data and Usage###

Example I
exampleI.R
Download and read exampleI.R first. This script contains several very simple lines of codes for creating a geno and a subpop object, and their usages in the following scripts: 
    calc_wcFstats(geno, subpop)
    calc_wcFst_spop_pairs(geno, subpop)
    calc_neiFis_onepop(geno)
    calc_snp_stats(geno)
    calc_neiFis_multispop(geno,subpop)
    calc_LD(geno)
    calc_allele_sharing(geno)  [See NJ tree of using resulting allele-sharing distance matrix.]
    calc_hwe_chisq(geno)
    calc_hwe_fisher(geno)
exampleI_data.RData  
A R workspace containing an instance of a geno and subpop objects used in exampleI.R; i.e. the actual datasets corresponding to the outputs in exampleI.R.
exampleI_functions.RData  
A R workspace containing all functions used in exampleI.R.

Example II
geno <- simgeno()
alleleCount <- geno_to_allelecnt(geno)
