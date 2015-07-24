calc_wcFst_spop_pairs <- function(geno, spop, plot.nj=F) {

    ## Copyright Eva Chan 2008
    ## eva@evachan.org
    ##
    ## A script to estimate Fst (theta) values for each pair of sub-populations
    ## using the method of Weir & Cockerham 1984 Evolution 38(6): 1358-1370.
    ##
    ## Arguments
    ## =========
    ## geno:   matrix of genotypes with rows corresp. to markers and columns to individuals;
    ##         notation for genotyeps are {0,1,2} indicating the number of one of the two alleles
    ## subpop: vector indicting the sub-popln to which the individuals belong to
    ## plot.nj: logical indicating whether a Neighbouring-Joining Tree of the results should be plotted.
    ##          (note that the R/ape library is required for this); defaults to FALSE
    ##
    ## Side effects
    ## ============
    ## output: Symmetical matix (in which only upper triangle is filled) of theta (Fst) values for each
    ## pair of unique sub-populaitons.
    ## plot: neighbouring-joining tree of results. 

    ## assign all non {0,1,2} to NA
    geno[(geno!=0) & (geno!=1) & (geno!=2)] <- NA
    geno <- as.matrix(geno)

    N = ncol(geno)   ## sample size
    if( length(spop) != N ) { stop( "Number of samples with genotypes does not match provided number of spop.\n" ) }
    spop <- as.factor(as.character(spop))
    unique.spop <- levels(spop)
    nspop = length(unique.spop)

    n0 <- n1 <- n <- matrix(NA, ncol=nspop, nrow=nrow(geno), dimnames=list(NULL,unique.spop))
    for(i in 1:nspop) {
        inds <- which(spop == unique.spop[i])
        n0[,i] <- apply(geno[,inds]==0,1,sum,na.rm=T)
        n1[,i] <- apply(geno[,inds]==1,1,sum,na.rm=T)
        n[,i] <- apply(!is.na(geno[,inds]),1,sum,na.rm=T)
    }
    p <- ((2*n0)+n1)/(2*n)   ## allele freq
    Ho <- (n1/n)             ## observed het

    pairwise.wcFst <- matrix(NA, ncol=nspop, nrow=nspop)
    r=2   ## now, only two spops are examined at a time
    for( i in 1:(nspop-1) ) {
        for( j in (i+1):nspop ) {

            n_bar <- apply(n[,unique.spop[c(i,j)]],1,sum,na.rm=T)/r
            nc <- ((r*n_bar) - (apply((n[,unique.spop[c(i,j)]]*n[,unique.spop[c(i,j)]])/(r*n_bar),1,sum,na.rm=T))) / (r-1)
            p_bar <- apply( (n[,unique.spop[c(i,j)]]*p[,unique.spop[c(i,j)]])/(r*n_bar), 1, sum, na.rm=T )
            s_square <- apply( (n[,unique.spop[c(i,j)]]*((p[,unique.spop[c(i,j)]]-p_bar)^2)) / ((r-1)*n_bar), 1, sum, na.rm=T )
            h_bar <- apply((n[,unique.spop[c(i,j)]]*Ho[,unique.spop[c(i,j)]])/(r*n_bar), 1, sum, na.rm=T)
            
            a_hat <- (n_bar/nc) * ( s_square - ((1/(n_bar-1))*((p_bar*(1-p_bar)) - (((r-1)/r)*s_square) - ((1/4)*h_bar))) )
            b_hat <- (n_bar/(n_bar-1)) * ((p_bar*(1-p_bar)) - (((r-1)/r)*s_square) - ((((2*n_bar)-1)/(4*n_bar))*h_bar))
            c_hat <- h_bar/2
            
            inds <- which(is.finite(a_hat) & is.finite(b_hat) & is.finite(c_hat))
            pairwise.wcFst[i,j] <- sum(a_hat[inds],na.rm=T) / sum(apply(cbind(a_hat,b_hat,c_hat)[inds,],1,sum,na.rm=T),na.rm=T)
            
            rm(n_bar, nc, p_bar, s_square, h_bar, a_hat, b_hat, c_hat,inds)
        }
    }
    colnames(pairwise.wcFst) <- rownames(pairwise.wcFst) <- unique.spop

    ## plot Neighbouring-Joining Tree
    if(plot.nj) {
        library(ape)
        pairwise.wcFst.nj <- nj(as.dist(t(pairwise.wcFst)))
        plot(pairwise.wcFst.nj, main="Weir & Cockerham's Fst",sub="neighbor joining",type="unrooted")
    }

    pairwise.wcFst

}
