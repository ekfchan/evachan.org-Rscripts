calc_hwe_chisq <- function(geno) {

    ## Copyright Eva Chan 2008
    ## http://evachan.org
    ## 
    ## This is a function for testing the significance of deviation from HWE using Pearson's Chi-Squared test.
    ## chisq = [((Obs(AA)-Exp(AA))^2)/Exp(AA)] + [((Obs(Aa)-Exp(Aa))^2)/Exp(Aa)] + [((Obs(aa)-Exp(aa))^2)/Exp(aa)]
    ## df = 1 (# phenotypes - # alleles; i.e. 3 genotypes - 2 alleles)
    ##
    ## Input:
    ## geno: SNP-by-sample matrix of genotypes {0,1,2}; any other values are ignored. 
    ##
    ## Output: two column matrix of Chi-square values and corresponding P-values for each SNP in geno.

    ## assign all non {0,1,2} to NA
    geno[(geno!=0) & (geno!=1) & (geno!=2)] <- NA
    geno <- as.matrix(geno)

    n0 <- apply(geno==0,1,sum,na.rm=T)
    n1 <- apply(geno==1,1,sum,na.rm=T)
    n2 <- apply(geno==2,1,sum,na.rm=T)
    n <- n0+n1+n2
    obs <- cbind(n0, n1, n2)
    p <- ((2*n0)+n1)/(2*n)
    q <- (1-p)
    expected <- cbind(p*p, 2*p*q, q*q)
    expected <- expected*n
    chisq <- (obs-expected)
    chisq <- (chisq*chisq) /expected
    chisq <- apply(chisq,1,sum)
    chisq.p <- 1-pchisq(chisq,df=1)
    
    res <- cbind(chisq=chisq, chisq.p=chisq.p)
    rownames(res) <- rownames(geno)
    res

}
