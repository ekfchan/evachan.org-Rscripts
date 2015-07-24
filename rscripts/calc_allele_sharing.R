calc_allele_sharing <- function(geno)
{
    ## Eva KF Chan
    ## http://evachan.org
    ##
    ## A script to calculate allele sharing distance between pairs of individuals
    ## (c.f. Gao & Stramer 2007 BMC Genetics 8:34)
    ## D_ij = (1/L) * sum(d_ij(l))
    ## where d_ij(l) = { 0 if individuals i & j have 2 alleles in common at l-th locus
    ##                 { 1 if individuals i & j have only 1 allele in common at l-th locus
    ##                 { 2 if individuals i & j have no allele in common at l-th locus
    ## and L = number of SNP loci.
    ##
    ## Input:
    ## geno: SNP-by-sample matrix of genotypes {0,1,2}; any other values are ignored.
    ##
    ## Output: 
    ## symmetrical matix of allele-sharing distance between each pair of individuals (columns of geno)
    ##
    ## NOTE:: if one wants to use this distance matrix to obtain Ward's Minimum Variance 
    ##        Hierarchical Clustering as in Gao & Stramer 2007, simply use the following 
    ##        command:
    ##        plot(allele.sharing.hclust <- hclust(as.dist(allele.sharing), method="ward"))

    n <- ncol(geno)                      ## number of individuals
    d <- matrix(NA, ncol=n, nrow=n, dimnames=list(colnames(geno),colnames(geno)))      ## distance
    
    for(i in 1:n)
    {
        cat(i,"\n")
        z <- abs(geno - geno[,i])
        d[,i] <- apply(z, 2, mean, na.rm=T)
    }

    d

}
