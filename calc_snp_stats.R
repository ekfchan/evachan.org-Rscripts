calc_snp_stats <- function(geno)
{
     ## Eva KF Chan
     ## http://evachan.org
     ##
     ## Created: 21/08/07
     ## Last Modified: 21/10/12
     ##
     ## Function to calculate basic stats on SNPs, including: allele frequency, MAF, and exact estimate of HWE
     ##
     ## geno: snp-by-individual matrix of genotypes, {0,1,2}.
     ##       NOTE:: any other values are ignored
     ##
     ## OUTPUT: data.frame of 
     ##         n, n0, n1, n2: number of samples with total non-missing genotype, and geno=0,1,or 2
     ##         p: allele frequency
     ##         maf & mgf: minor allele & genotype frequencies
     ##         mono: {T,F} indicating if marker is monomorphic (MAF<0%)
     ##         loh: {T,F} indicating if marker has loss of heterozygote
     ##         hwe.chisq & hwe.chisq.p: chi-square test statistic for deviation from HWE and correp p-value
     ##         hwe.fisher & hwe.fisher.p: Fisher's Exact test statistic for deviation from HWE and correp p-value
     ##

     m <- nrow(geno)     ## number of snps
     n <- ncol(geno)     ## number of individuals

     ## assign all non {0,1,2} to NA
     geno[(geno!=0) & (geno!=1) & (geno!=2)] <- NA
     geno <- as.matrix(geno)

     ## calc_n
     n0 <- apply(geno==0,1,sum,na.rm=T)
     n1 <- apply(geno==1,1,sum,na.rm=T)
     n2 <- apply(geno==2,1,sum,na.rm=T)

     n <- n0 + n1 + n2

     ## calculate allele frequencies
     p <- ((2*n0)+n1)/(2*n)
     q <- 1 - p
     maf <- pmin(p, q)
     mgf <- apply(cbind(n0,n1,n2),1,min) / n

     ## HWE: Chi-Square test
     obs <- cbind(n0=n0,n1=n1,n2=n2)
     exp <- cbind(p*p, 2*p*q, q*q)
     exp <- exp*n
     chisq <- (obs-exp)
     chisq <- (chisq*chisq) /exp
     hwe.chisq <- apply(chisq,1,sum)
     hwe.chisq.p <- 1-pchisq(hwe.chisq,df=1)

     ## HWE: Fisher's Exact test
     z <- cbind(n0, ceiling(n1/2), floor(n1/2), n2)
     z <- lapply( split( z, 1:nrow(z) ), matrix, ncol=2 )
     z <- lapply( z, fisher.test )
     hwe.fisher <- as.numeric(unlist(lapply(z, "[[", "estimate")))
     hwe.fisher.p <- as.numeric(unlist(lapply(z, "[[", "p.value")))

	# MODIFIED 21 Oct 2012:  prior to this version, we had "mono=(mgf<0)" instead of "mono<(maf<0)"
     res <- data.frame( n=n, n0=n0, n1=n1, n2=n2, p=p, maf=maf, mgf=mgf,
                        mono=(maf<=0), loh=(n1<=0), 
                        hwe.chisq=hwe.chisq, hwe.chisq.p=hwe.chisq.p,
                        hwe.fisher=hwe.fisher, hwe.fisher.p=hwe.fisher.p, 
                        stringsAsFactors=F )
     row.names(res) <- row.names(geno)
     res
}