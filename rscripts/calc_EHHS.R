calc_EHHS <- function(geno, thresh=0.1) {

     ## March 2010
     ## Eva KF Chan 
     ## http://evachan.org 
     ##
     ## Function to calculate the EHHS(geno)i,j values for a given chromosome as described in:
     ## Tang K, Thornton KR, Stoneking M (2007) A New Approach for Using Genome Scans to Detect Recent Positive Selection in the Human Genome . PLoS Biol 5(7): e171 doi:10.1371/journal.pbio.0050171
     ##
     ## EHHS is the haplotype homozygosity between sites i and j, normalised by the homozygosity at site i:
     ## EHHS(geno)i,j = sum_k=1..n{ Ik,(hap ij) [1 if hap1 = hap2] } / suml=1..n{ I_i,(alle i) [1 if alle1 = alle2]}
     ## where:
     ##		Ik,(hap ij) = identity of the two haplotypes between site i & j in one individual
     ##		Il,(alle i) = identity of the alleles at site i
     ##
     ##                  EHHS_(geno)i,j = ( sum(k=1..n) {I_k,(hap ij)  [1 if hap1 = hap2] )
     ##                                   -------------------------------------------------
     ##					  ( sum(l=1..n) {I_l,(alle i)  [1 if alle1 = alle2] )
     ##					  
     ##					=      number of individuals where hap1 = hap2 
     ##					  ---------------------------------------------------
     ##					  number of individuals where alle1 = alle2 at site i
     ## Parameters:
     ##		geno: matrix of genotypes (0,1,2,NA) of size marker (row) by sample (column)
     ##		      <<MAKE SURE THE MARKERS ARE IN GENOMIC ORDER!!!>>
     ##		thresh: the threshold [0,1] to which EHHS is calculated for all j moving away 
     ##                 from site i until EHHS < thresh (0.1 by default)
     ## 
     ## Output: 
     ## 	Returns a matrix of size MxM (M=number of markers=nrow(geno)) of EHHS values 
     ##         calcualted for all i-th marker (row) to each j-th marker (colum) until EHH < thresh

     geno[geno==2] <- 0		## 0=homozygous; 1=heterozygous
     M = nrow(geno)
     EHH <- matrix(NA, ncol=M, nrow=M, dimnames=list(rownames(geno),rownames(geno)))

     for(i in 1:M) {
	  initial_list = which(geno[i,]==0)
	  Ii = length(initial_list)
	  EHH[i,i] = 1

	  ## left-flank
	  cur_list = initial_list
	  j = i-1
	  while( j >= 1 ) {
	       tmp_list = which(geno[j,]==0)
	       cur_list = intersect( cur_list, tmp_list )
	       Ij = length(cur_list)
	       cur.EHH = Ij/Ii
	       if (is.na(cur.EHH) | cur.EHH < thresh) { break } else {
		    EHH[i,j] = cur.EHH
	       }
	       j = j-1
	  }

	  ## right-flank
	  cur_list = initial_list
	  j = i + 1
	  while( j <= M ) {
	       tmp_list = which(geno[j,]==0)
	       cur_list = intersect( cur_list, tmp_list )
	       Ij = length(cur_list)
	       cur.EHH = Ij/Ii
	       if (is.na(cur.EHH) | cur.EHH < thresh) { break } else {
		    EHH[i,j] = cur.EHH
	       }
	       j = j+1
	  }
     }

     return(EHH)

}
