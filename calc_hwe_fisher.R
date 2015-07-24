calc_hwe_fisher <- function(geno) {
    
    ## Copyright Eva Chan 2008
    ## htp://evachan.org
    ##
    ## This is a function for testing the significance of deviation from HWE
    ## using Fisher's Exact test.
    ## Note that the observed number of Aa and aA genotypes are identical if
    ## their sum is even, else Aa is always one more than aA.
    ##
    ## Input:
    ## geno: SNP-by-sample matrix of genotypes {0,1,2}; any other values are ignored. 
    ##
    ## Output: two column matrix of Odds Ratio and corresponding P-values for each SNP in geno.

    ## assign all non {0,1,2} to NA
    geno[(geno!=0) & (geno!=1) & (geno!=2)] <- NA
    geno <- as.matrix(geno)

    n0 <- apply(geno==0, 1, sum, na.rm=T)
    n1 <- apply(geno==1, 1, sum, na.rm=T)
    n2 <- apply(geno==2, 1, sum, na.rm=T)
    
    z <- cbind(n0, ceiling(n1/2), floor(n1/2), n2)
    z <- lapply( split( z, 1:nrow(z) ), matrix, ncol=2 )
    z <- lapply( z, fisher.test )
    
    res <- cbind( odds.ratio = as.numeric(unlist(lapply(z, "[[", "estimate"))),
           p.values = as.numeric(unlist(lapply(z, "[[", "p.value"))) )
    rownames(res) <- rownames(geno)
    res

}
