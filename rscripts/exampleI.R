
### ----------------------------------- ###
###            Create Data              ###
### ----------------------------------- ###

## Create data of 20 individuals belonging to three subpopulations (A,B,and C)
subpop <- rep(LETTERS[1:3],c(8,5,7) )

## Simulate genotype matrix of 100 markers by 20 individuals, based on subpop. 
## Simulation assumes the 100 markers are unlinked.
geno <- matrix( NA, ncol=20, nrow=100,
dimnames=list(paste('M',sapply(mapply(rep,3-nchar(1:100),MoreArgs=list(x=0)),paste,collapse=''),1:100,sep=''),paste('S',sapply(mapply(rep,2-nchar(1:20),MoreArgs=list(x=0)),paste,collapse=''),1:20,sep='')))

for( i in 1:100 ) {
	for( s in LETTERS[1:3] ) {
		sinds <- which(subpop==s)
		tmp.prob = runif(1)
		geno[i,sinds] <-
			sample(0:1,length(sinds),replace=T,p=c(tmp.prob,1-tmp.prob)) +
			sample(0:1,length(sinds),replace=T,p=c(tmp.prob,1-tmp.prob))
			rm(sinds, tmp.prob)
	}; rm(s)
}; rm(i)


### ----------------------------------- ###
###            Example Usage            ###
### ----------------------------------- ###

### calc_wcFstats.R
### ===============
wcFstats <- calc_wcFstats(geno,subpop)
sapply(wcFstats,head)	#peek at the output
	#$perloc
	#           a_hat        b_hat c_hat      F_hat   theta_hat        f_hat
	#[1,] 0.134044839 -0.025315126 0.175 0.38321581 0.472438497 -0.169122807
	#[2,] 0.066264674  0.001680672 0.200 0.25357912 0.247306679  0.008333333
	#[3,] 0.252778401  0.034926471 0.050 0.85194173 0.748518669  0.411255411
	#[4,] 0.134213227 -0.045168067 0.200 0.30806660 0.464333073 -0.291723202
	#[5,] 0.001800148  0.014548319 0.200 0.07556544 0.008320593  0.067809058
	#[6,] 0.089654564 -0.038130252 0.125 0.29188224 0.507887911 -0.438935913
	#
	#$global
	#      F_hat   theta_hat       f_hat 
	#0.343454963 0.341996174 0.002216991 


### calc_wcFst_spop_pairs.R
### =======================
wcFstats.spop.pairs <- calc_wcFst_spop_pairs( geno, subpop, TRUE )
wcFstats.spop.pairs	#output:  FST between the three pairs of subpopulations: A vs B, A vs C, and B vs C.
	#   A         B         C
	#A NA 0.3850011 0.3405093
	#B NA        NA 0.2961459
	#C NA        NA        NA


### calc_neiFis_onepop.R
### ====================
neiFis.popA <- calc_neiFis_onepop(geno[,which(subpop=='A')])
sapply(neiFis.popA,head)	#peek at the output:
	#$aveloc
	#[1] 0.05077399
	#
	#$perloc
	#       M001        M002        M003        M004        M005        M006 
	#        NaN -0.27272727  0.00000000 -0.07692308  0.30000000         NaN 
geno[1:6,which(subpop=='A')]	#Note:  markers at which FIS cannot be calculated have NA (missing) values. In this example, markers M001 & M006 are monomorphic in subpopulation A.
	#     S01 S02 S03 S04 S05 S06 S07 S08
	#M001   0   0   0   0   0   0   0   0
	#M002   1   0   1   0   0   1   1   0
	#M003   0   1   0   0   0   0   0   0
	#M004   1   2   1   2   2   2   2   2
	#M005   0   1   1   0   2   1   0   2
	#M006   2   2   2   2   2   2   2   2


### calc_snp_stats.R
### ================
snp.stats <- calc_snp_stats(geno)
head(snp.stats)	#peek at the output
	#      n n0 n1 n2     p   maf  mgf  mono   loh   hwe.chisq  hwe.chisq.p hwe.fisher hwe.fisher.p
	#M001 20  9  7  4 0.625 0.375 0.20 FALSE FALSE  1.28355556 0.2572389852   2.828469 0.3563467492
	#M002 20  8  8  4 0.600 0.400 0.20 FALSE FALSE  0.55555556 0.4560565403   1.929814 0.6479161705
	#M003 20  8  2 10 0.450 0.450 0.10 FALSE FALSE 12.73543516 0.0003587923  50.725102 0.0009228388
	#M004 20  4  8  8 0.400 0.400 0.20 FALSE FALSE  0.55555556 0.4560565403   1.929814 0.6479161705
	#M005 20 10  8  2 0.700 0.300 0.10 FALSE FALSE  0.04535147 0.8313590555   1.235838 1.0000000000
	#M006 20  1  5 14 0.175 0.175 0.05 FALSE FALSE  0.36018815 0.5484017591   2.218801 0.5087719298

# Some of the statistics are actually not appropriate in the presence of known population substrucutre. 
# In our case, since our *geno* corresponded to individuals from three different subpopulations, we may wish to re-estimate the marker statistics on a by-population basis.
snp.stats.popA <- calc_snp_stats(geno[,which(subpop=='A')])	#statistics for subpopulation A
head(snp.stats.popA)	#Note that markers M001, M006 (& others) were monomorphic and so chi-square could not be calculated. 
	#     n n0 n1 n2      p    maf  mgf  mono   loh  hwe.chisq hwe.chisq.p hwe.fisher hwe.fisher.p
	#M001 8  8  0  0 1.0000 0.0000 0.00  TRUE  TRUE        NaN         NaN   0.000000            1
	#M002 8  4  4  0 0.7500 0.2500 0.00 FALSE FALSE 0.88888889   0.3457786   0.000000            1
	#M003 8  7  1  0 0.9375 0.0625 0.00 FALSE FALSE 0.03555556   0.8504363   0.000000            1
	#M004 8  0  2  6 0.1250 0.1250 0.00 FALSE FALSE 0.16326531   0.6861678   0.000000            1
	#M005 8  3  3  2 0.5625 0.4375 0.25 FALSE FALSE 0.45351474   0.5006706   2.601683            1
	#M006 8  0  0  8 0.0000 0.0000 0.00  TRUE  TRUE        NaN         NaN   0.000000            1


### calc_neiFis_multispop.R
### =======================
neiFis.allpops <- calc_neiFis_multispop(geno,subpop)
sapply(neiFis.allpops,head)	#peek at output
	#$aveloc	## average FIS for each subpopulation, of combined pop, and averaged across the 3 subpops
	#           A            B            C        total      average 
	# 0.050773994 -0.096671949  0.017681729  0.265446224 -0.008713363 
	#
	#$perloc	## FIS values for EACH marker
	#               A             B             C      total     average
	#[1,]         NaN -6.000000e-01  1.428571e-01 0.27717391 -0.22953762
	#[2,] -0.27272727 -1.428571e-01  3.684211e-01 0.19148936  0.00980153
	#[3,]  0.00000000  6.000000e-01           NaN 0.80710660  0.46406559
	#[4,] -0.07692308 -2.220446e-16 -5.000000e-01 0.19148936 -0.28099627
	#[5,]  0.30000000 -1.428571e-01 -2.000000e-01 0.07317073  0.03715450
	#[6,]         NaN -6.000000e-01  2.220446e-16 0.15929204 -0.47610485


### calc_LD.R
## ==========
LDrsq.allpops.SNPs1to50 <- calc_LD( geno, inds=1:50, get.D=F, get.Dprime=F, get.rsq=T, get.chisq=F, get.chisq_prime=F )
sapply(LDrsq.allpops.SNPs1to50,'[',1:6,1:6)	#peek at the output: notice that only $rsq is not NULL as we had asked calc_LD to return only r-square values. 
	#$D
	#NULL
	#
	#$Dprime
	#NULL
	#
	#$rsq	##symmetircal matrix of r-squares estimated for all 100x100 pairs of markers
	#            [,1]        [,2]        [,3]      [,4]       [,5]       [,6]
	#[1,] 1.000000000 0.006734007 0.486531987 0.2045455 0.01010101 0.04306220
	#[2,] 0.006734007 1.000000000 0.001736111 0.1666667 0.16666667 0.03508772
	#[3,] 0.486531987 0.001736111 1.000000000 0.0234375 0.00000000 0.03508772
	#[4,] 0.204545455 0.166666667 0.023437500 1.0000000 0.00000000 0.21052632
	#[5,] 0.010101010 0.166666667 0.000000000 0.0000000 1.00000000 0.05263158
	#[6,] 0.043062201 0.035087719 0.035087719 0.2105263 0.05263158 1.00000000
	#
	#$chisq
	#NULL
	#
	#$chisq_prime
	#NULL
	#
	#$chisq_df
	#NULL

# LD estimated using samples with known population structure may not be appropriate. 
# In our case here, we may wish to estimate LD using only samples from subpopulation A. 
LDrsq.popA.SNPs1to50 <- calc_LD( geno[,which(subpop=='A')], inds=1:50, get.D=F, get.Dprime=F, get.rsq=T, get.chisq=F, get.chisq_prime=F )
LDrsq.popA.SNPs1to50$rsq[1:6,1:6]	#peek at output:  Note the large numbers of missing data due to (1) small sample size and (2) monomophic loci
	#     [,1]       [,2]       [,3] [,4]       [,5] [,6]
	#[1,]  NaN        NaN        NaN  NaN        NaN  NaN
	#[2,]  NaN 1.00000000 0.14285714  NaN 0.06666667  NaN
	#[3,]  NaN 0.14285714 1.00000000  NaN 0.08571429  NaN
	#[4,]  NaN        NaN        NaN  NaN        NaN  NaN
	#[5,]  NaN 0.06666667 0.08571429  NaN 1.00000000  NaN
	#[6,]  NaN        NaN        NaN  NaN        NaN  NaN


### calc_allele_sharing.R
### =====================
ASdist.allpop <- calc_allele_sharing(geno)
rownames(ASdist.allpop) <- colnames(ASdist.allpop) <- colnames(geno)	#Note that calc_allele_sharing does not preserve marker labels. 

#To use allele-sharing distance matrix as input for estiamting Neighbour-Joining Tree: 
library(ape)	#the library with nj()
ASdist.allpop.nj <- nj(ASdist.allpop)
plot( ASdist.allpop.nj, tip.color=match(subpop,LETTERS[1:3]), type='unroot' )


### calc_hwe_chisq.R
### ================
HWEchisq.allpop <- calc_hwe_chisq(geno)
head(HWEchisq.allpop)	#peek at output
	#           chisq      chisq.p
	#M001  1.28355556 0.2572389852
	#M002  0.55555556 0.4560565403
	#M003 12.73543516 0.0003587923
	#M004  0.55555556 0.4560565403
	#M005  0.04535147 0.8313590555
	#M006  0.36018815 0.5484017591

# Testing of deviation from HWE in data with known population substructure may not be appropriate. 
# Here we re-estimate with only data from subpopulation A.
HWEchisq.popA <- calc_hwe_chisq(geno[,which(subpop=='A')])
head(HWEchisq.popA)	#Note: chi-square cannot be calculated on monomorphic markers (e.g. M001 & M006)
	#          chisq   chisq.p
	#M001        NaN       NaN
	#M002 0.88888889 0.3457786
	#M003 0.03555556 0.8504363
	#M004 0.16326531 0.6861678
	#M005 0.45351474 0.5006706
	#M006        NaN       NaN


### calc_hwe_fisher.R
### =================
HWEfisher.allpop <- calc_hwe_fisher(geno)
head(HWEfisher.allpop)	#peek at output
	#     odds.ratio     p.values
	#[1,]   2.828469 0.3563467492
	#[2,]   1.929814 0.6479161705
	#[3,]  50.725102 0.0009228388
	#[4,]   1.929814 0.6479161705
	#[5,]   1.235838 1.0000000000
	#[6,]   2.218801 0.5087719298

# Testing of deviation from HWE in data with known population substructure may not be appropriate. 
# Here we re-estimate with only data from subpopulation A.
HWEfisher.popA <- calc_hwe_fisher(geno[,which(subpop=='A')])
head(HWEfisher.popA)
	#     odds.ratio p.values
	#[1,]   0.000000        1
	#[2,]   0.000000        1
	#[3,]   0.000000        1
	#[4,]   0.000000        1
	#[5,]   2.601683        1
	#[6,]   0.000000        1

