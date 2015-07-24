calc_neiFis_onepop <- function (geno) {
    
    ## Copyright Eva Chan 2008
    ## eva@evachan.org
    ##
    ## A script to calculate inbreeding coefficients, Fis (Nei 1977 Ann Hum Genet 41:225-233), 
    ## for total population from a given set of SNP markers.
    ##
    ## Input:
    ## geno: SNP-by-sample matrix of genotypes {0,1,2}; any other values are ignored.
    ##
    ## Output: list of 
    ##         1) aveloc: single Fis value averaged over all loci for the given population
    ##         2) perloc: numeric vector of Fis per SNP (row) for the given population

    ## assign all non {0,1,2} to NA
    geno[(geno!=0) & (geno!=1) & (geno!=2)] <- NA
    geno <- as.matrix(geno)
    
    m = nrow(geno)   ## number of markers
    N = ncol(geno)   ## number of samples
    
    ## determine numbers of each genotypes for the given pop at each locus
    nNA <- apply(is.na(geno),1,sum)
    nAA <- apply(geno==0,1,sum,na.rm=T)
    nAa <- apply(geno==1,1,sum,na.rm=T)
    naa <- apply(geno==2,1,sum,na.rm=T)
    n <- nAA + nAa + naa

    Ho <- (nAa/n)                                   ## observed het
    p <- ((2*nAA)+nAa)/(2*n)                        ## allele freq
    He <- (n/(n-1)) * ((2*p*(1-p)) - (Ho/(2*n)))    ## Nei's expected het
    
    ncFis <- 1 - (mean(Ho,na.rm=T) / mean(He,na.rm=T))
    ncFis.perloc <- 1 - (Ho/He)
    
    list(aveloc = ncFis, perloc = ncFis.perloc)

}
