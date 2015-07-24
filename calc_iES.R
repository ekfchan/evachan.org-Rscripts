calc_iES <- function(EHHS, lox) {

     ## Eva KF Chan 
     ## Nov 2008
     ## http://evachan.org
     ##
     ## A function to calculate the iES statistics as per Tang K, Thornton KR, Stoneking M (2007) A New Approach for Using Genome Scans to Detect Recent Positive Selection in the Human Genome . PLoS Biol 5(7): e171 doi:10.1371/journal.pbio.0050171
     ## iES integrates the area under the curve of EHHS against distance: 
     ## iES_i = sum_for_j_from_a+1_to_b = { (EHHS_i,j-1 + EHHS_i,j) * (Pos_j - Pos_j-1) } / 2
     ## where:
     ## 	a & b are the two ending positions where EHHS < X 
     ## 	Pos_j is the physical position of site j
     ##
     ## Parameters:
     ##		EHHS: matrix of size MxM (M=number of markers=nrow(geno)) of EHHS values calcualted for all i-th marker (row) to each j-th marker (colum) until EHH < thresh
     ##		      <<MARKERS SHOULD ALREADY BE IN GENOMIC ORDER!!>>
     ##		lox: genomic location of the markers (row) in EHHS; THIS SHOULD BE IN SAME ORDER AS EHHS
     
     if( nrow(EHHS) != length(lox) ) { stop("Number of positions given does not agree with number of markers.\n") }
     
     M <- length(lox) 
     iES <- rep(NA, M)
     
     x = lox[2:M] - lox[1:(M-1)]
     for(i in 1:M) {
          y = EHHS[i,1:(M-1)] + EHHS[i,2:M]
          if( !all(is.na(y)) ) {
               iES[i] = sum(y*x, na.rm=T) / 2
          }; rm(y)
     }
     
     iES

}
