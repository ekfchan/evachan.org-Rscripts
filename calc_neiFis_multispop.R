calc_neiFis_multispop <- function (geno, spop) {
    
    ## Copyright Eva Chan 2008
    ## eva@evachan.org
    ##
    ## A script to calculate inbreeding coefficients, Fis (Nei 1977 Ann Hum Genet 41:225-233), 
    ## for each sub-population from a given set of SNP markers.
    ##
    ## Input:
    ## geno: SNP-by-sample matrix of genotypes {0,1,2}; any other values are ignored.
    ## spop: a factor indicating the sub-population to which the corresponding samples 
    ##       (columns) in geno belong.
    ##
    ## Output: list of 
    ##         1) aveloc: numeric vector of Fis averaged over all loci for each sub-population, 
    ##                    the total population (2nd last value),
    ##                    and average of total population (last value)
    ##         2) perloc: matrix of Fis per SNP (row) for each sub-population, 
    ##                    the total population (2nd last column),
    ##                    and average of total population (last column)

    ## assign all non {0,1,2} to NA
    geno[(geno!=0) & (geno!=1) & (geno!=2)] <- NA
    geno <- as.matrix(geno)
    
    m = nrow(geno)   ## number of markers
    N = ncol(geno)   ## number of samples
    
    if( length(spop) != N ) { stop( "Number of samples with genotypes does not match provided number of spop.\n" ) }
    spop <- as.factor(as.character(spop))
    unique.spop <- levels(spop)
    nspop <- length(unique.spop)

    ## determine numbers of each genotypes for each spop at each locus
    nNA <- nAA <- nAa <- naa <- matrix(NA, ncol=nspop, nrow=m, dimnames=list(NULL,unique.spop))
    for(i in 1:nspop) {
         inds <- which(spop == unique.spop[i])
         nAA[,i] <- apply(geno[,inds]==0,1,sum,na.rm=T)
         nAa[,i] <- apply(geno[,inds]==1,1,sum,na.rm=T)
         naa[,i] <- apply(geno[,inds]==2,1,sum,na.rm=T)
    }
    n <- nAA + nAa + naa
    nAA <- cbind(nAA, total=apply(nAA[,unique.spop],1,sum))
    nAa <- cbind(nAa, total=apply(nAa[,unique.spop],1,sum))
    naa <- cbind(naa, total=apply(naa[,unique.spop],1,sum))
    n <- cbind(n, total=apply(n[,unique.spop],1,sum))
    
    Ho <- (nAa/n)                                   ## observed het
    p <- ((2*nAA)+nAa)/(2*n)                        ## allele freq
    He <- (n/(n-1)) * ((2*p*(1-p)) - (Ho/(2*n)))    ## Nei's expected het
    
    s <- apply(!is.na(n[,unique.spop]),1,sum)        ## number of spop per marker
    n_tilda <- s/apply((1/n[,unique.spop]),1,sum)    ## harmonic mean of sample sizes
    Ho <- cbind(Ho, average=(apply(Ho[,unique.spop],1,sum,na.rm=T)/s))   ## Ho averged over samples
    He <- cbind(He, average=( (n_tilda/(n_tilda-1)) * ((apply(2*p[,unique.spop]*(1-p[,unique.spop]),1,sum,na.rm=T)/s) - (Ho[,"average"]/(2*n_tilda))) ))                                          ## Nei's averaged He
    
    ncFis <- 1 - (apply(Ho,2,mean,na.rm=T) / apply(He,2,mean,na.rm=T))
    ncFis.perloc <- 1 - (Ho/He)
    
    list(aveloc = ncFis, perloc = ncFis.perloc)

}
