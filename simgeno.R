simgeno <- function(M=100, N=30, propNA=0) {

	## http://evachan.org
	## Eva KF Chan
	## Created:  7 July 2014
	##
	## Very simple function to simulate a matrix of biallelic unphased SNP genotypes in the format: {AA,CC,GG,TT,AC,AG,AT,CG,CT,GT}.  This is written predominantly for the purpose of demonstrating the geno_toallelecnt.R function. 
	## Inputs:
	##	M:	The number of SNP markers to simulate.
	##	N:	The number of samples to simulate. 
	##	propNA:	The proportion of missing data to simulate. 
	## Output:
	##	A matrix of genotypes. 

	alleles = c('A','C','G','T')
	a1 <- sample(alleles,M,replace=T)
	a2 <- sample(alleles,M,replace=T)
	g0=paste(a1,a1,sep='')
	g1=paste(a1,a2,sep='')
	g2=paste(a2,a2,sep='')

	geno <- matrix( NA, ncol=N, nrow=M, dimnames=list(paste("marker",1:M,sep=""),paste("sample",1:N,sep="")) )
	for( i in 1:M ) { geno[i,] <- sample(c(g0[i],g1[i],g2[i]),N,replace=T) }
	if( propNA>0 ) { geno[sample(1:(M*N),ceiling(M*N*propNA))] <- NA }

	geno
}
