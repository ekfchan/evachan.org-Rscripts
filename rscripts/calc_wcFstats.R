calc_wcFstats <- function(geno, subpop) {

## Copyright Eva Chan 2008
## eva@evachan.org
##
## A script to estimate the variance components and fixation indices as described in
## Weir & Cockerham 1984 Evolution 38(6): 1358-1370.
##
## Arguments
## =========
## geno:   matrix of genotypes with rows corresp. to markers and columns to individuals; 
##         notation for genotyeps are {0,1,2} indicating the number of one of the two alleles
## subpop: vector indicting the sub-popln to which the individuals belong to
##
## Output
## ======
## list of two objects: perloc and global
## perloc: matrix of 6 columns and as many rows as markers in geno
##         the 6 columns contain the estiamted variance components and fixation indices per locus
##         a = component of variance between subpops
##         b = component of variance between individuals within subpops
##         c = component of variance between gametes within individuals
## global: numeric vector of three values corresponding to the esimated F (Fit), theta (Fst), & 
##         f (Fis) across all loci
##
## Note
## ====
## R/HIERFSTAT also estimate F-statistics using variance component estimation.
## Results from that package is not too different to those from this function; I suspect
## there are two sources of differences:
## 1) all estimates of variance component from HIERFSTAT are doubled in magnitude to those
##    from this function (i.e. scaled by factor of 2);
## 2) rounding off variations may also be present.
## The scaled difference in the estimates of variance components poses no problem when
## calcualting fixation indicies as the values scaling factor is cancelled out in the
## calculation of the ratios.

     spop <- unique(as.character(subpop))    ## unique spops
     r <- length(spop)

     n11 <- n12 <- n22 <- matrix(NA, ncol=r, nrow=nrow(geno))
     for(i in 1:r) {
          inds <- which(subpop == spop[i])
          n11[,i] <- apply(geno[,inds]==0,1,sum,na.rm=T)
          n12[,i] <- apply(geno[,inds]==1,1,sum,na.rm=T)
          n22[,i] <- apply(geno[,inds]==2,1,sum,na.rm=T)
     }
     ni <- n11 + n12 + n22
     pi_tilda <- ((2 * n11) + n12) / (2 * ni)
     hi_tilda <- n12 / ni
     n_bar <- apply(ni,1,sum,na.rm=T)/r
#      C_square <- ( apply(ni*ni,1,sum,na.rm=T) - (n_bar*n_bar*r) ) / ( (n_bar*n_bar) * (r-1) )    ## mod 2/3/2008
#      nc <- n_bar * (1 - (C_square/r))
     nc <- ((r*n_bar) - apply(((ni*ni)/(r*n_bar)),1,sum,na.rm=T)) / (r - 1)
     p_bar <- apply( (ni*pi_tilda)/(r*n_bar), 1, sum, na.rm=T )
     s_square <- apply( (ni*((pi_tilda-p_bar)^2)) / ((r-1)*n_bar), 1, sum, na.rm=T )
     h_bar <- apply((ni*hi_tilda)/(r*n_bar), 1, sum, na.rm=T)

#      F_hat = 1 - ( (h_bar*(1-(C_square/r))) / ( (2*p_bar*(1-p_bar)*(1-((n_bar*C_square)/(r*(n_bar-1))))) + (2*(s_square/r)*(1+(((r-1)*(n_bar*C_square))/(r*(n_bar-1))))) + ((h_bar/2)*(C_square/(r*(n_bar-1)))) ))
#      theta_hat <- (s_square - ((1/(n_bar-1))*((p_bar*(1-p_bar)) - (((r-1)/r)*s_square) - (h_bar/4)))) / (((1-((n_bar*C_square)/(r*(n_bar-1))))*p_bar*(1-p_bar)) + ((1+(((r-1)*n_bar*C_square)/(r*(n_bar-1))))*(s_square/r)) + ((C_square/(r*(n_bar-1)))*(h_bar/4)))
#      f_hat <- 1 - (h_hat / ((((2*n_bar)/(n_bar-1))*p_bar*(1-p_bar)) - (((2*n_bar*(r-1))/(r*(n_bar-1)))*s_square) - ((1/(n_bar-1))*(h_bar/2))))

     a_hat <- (n_bar/nc) * ( s_square - ((1/(n_bar-1))*((p_bar*(1-p_bar)) - (((r-1)/r)*s_square) - ((1/4)*h_bar))) )
     b_hat <- (n_bar/(n_bar-1)) * ((p_bar*(1-p_bar)) - (((r-1)/r)*s_square) - ((((2*n_bar)-1)/(4*n_bar))*h_bar))
     c_hat <- h_bar/2

     F_hat <- 1 - (c_hat / (a_hat + b_hat + c_hat))
     theta_hat <- a_hat / (a_hat + b_hat + c_hat)
     f_hat <- 1 - (c_hat / (b_hat + c_hat))

     F_hat_w <- 1 - (sum(c_hat,na.rm=T) / sum((a_hat + b_hat + c_hat),na.rm=T))
     theta_hat_w <- sum(a_hat,na.rm=T) / sum((a_hat + b_hat + c_hat),na.rm=T)
     f_hat_w <- 1 - (sum(c_hat,na.rm=T) / sum((b_hat + c_hat),na.rm=T))
     
     list( perloc=cbind(a_hat=a_hat, b_hat=b_hat, c_hat=c_hat, F_hat=F_hat, theta_hat=theta_hat, f_hat=f_hat),
           global=c(F_hat=F_hat_w, theta_hat=theta_hat_w, f_hat=f_hat_w) )

}

