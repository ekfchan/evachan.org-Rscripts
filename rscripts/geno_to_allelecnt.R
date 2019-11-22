geno_to_allelecnt <- function(geno, ref=NULL, info=FALSE) {

	## http://evachan.org
	## Eva KF Chan
	## Created:  7 July 2014
	##
	## Converts a matrix of genotypes into a matrix of allele counts. 
	## It essentially converts the bi-allele SNP data format of {AA,AG,GG,CC,...} 
	## to the number of copies of the ref (or alphabetically "smaller") allele {0,1,2}
	## Inputs:
	##	geno:	Matrix of genotypes with rows corresponding to markers and columns 
	##			to samples. NA is allowed.
	##	ref:	Character vector of same size the number of rows in geno, representing
	##			the reference "allele".  If absent, then conversion will be based on 
	##			the alphabetically smaller allele. 
	## Output:
	##	If info is FALSE (default), the function returns a single matrix of the same size as geno, containing the counts of the reference/common allele at each marker (rows). 
	##	If info is TRUE, a list will be return containing the matrix of allele counts as well as a data.frame of marker information. This is useful for checking which alelles are counted. 

	if(!is.matrix(geno) | !mode(geno)=="character") { stop("geno must be of 'matrix' class and 'character' mode.\n") }
	if( !all(nchar(as.character(geno[!is.na(geno)]))==2) ) { stop("geno should contain bi-allelic genotypes, e.g. {AA,CC,GG,TT,AC,AG,AT,CG,CT,GT}\n") }

	markers <- data.frame( N=rowSums(!is.na(geno)) ) 
	
	alleles <- apply(cbind(substr(geno,1,1),substr(geno,2,2)),1,unique) 
	if( is.matrix(alleles) ) { alleles <- lapply(apply(alleles,2,as.list),as.character) }	#2017-03-15: corrected apply direction
	alleles <- lapply(alleles,sort)
	markers$numAlleles = sapply(alleles,length)
	if( any(markers$numAlleles>2) ) { stop("markers {",paste(which(markers$numAlleles>2),collapse=","),"} contains more than two alleles.\n") }

	markers$A1 = NA
	inds <- which(markers$numAlleles>0)
	markers$A1[inds] <- sapply(alleles[inds],'[[',1)
	markers$A2 = NA
	inds <- which(markers$numAlleles>1)
	markers$A2[inds] <- sapply(alleles[inds],'[[',2)

	if(is.null(ref)) { ref <- markers$A1; markers$input_ref=NA } else { markers$input_ref=ref }
	# If ref allele was not known, the alphabetically smaller allele is used
	if(length(inds<-which(is.na(ref)))>0) { ref[inds] = markers$A1[inds] }
	alt <- rep(NA,length(ref))
	inds <- which(ref==markers$A1); alt[inds] <- markers$A2[inds]
	inds <- which(ref==markers$A2); alt[inds] <- markers$A1[inds]
	inds <- which(is.na(alt)); alt[inds] = markers$A1[inds]	#if neither alleles is the reference, arbitrarily assign the alphabetically smaller allele as the alt
	markers$ref = ref 
	markers$alt = alt

	#if( any(ref!=markers$A1 & ref!=markers$A2) ) { warning("ref allele not present in geno for some markers. Conversions for these markers cannot be performed and will be coerced to NA.\n") }
	
	markers$G2 = paste(ref,ref,sep="")	#2 copies of ref
	markers$G1.1 = paste(ref,alt,sep="")	#1 copy of ref, ref allele coded first
	markers$G1.2 = paste(alt,ref,sep="")	#1 copy of ref, reversed coding
	markers$G0 = paste(alt,alt,sep="")	#0 copy of ref
	markers$G2[is.na(ref)] <- NA
	markers$G1.1[is.na(alt)] <- NA
	markers$G1.2[is.na(alt)] <- NA
	markers$G0[is.na(alt)] <- NA

	geno.as.num <- matrix( 0, ncol=ncol(geno), nrow=nrow(geno), dimnames=dimnames(geno) )
	geno.as.num[geno==markers$G2] <- 2 
	geno.as.num[geno==markers$G1.1 | geno==markers$G1.2] <- 1
	geno.as.num[geno==markers$G0] <- 0
	geno.as.num[which(is.na(markers$ref)),] = NA
	geno.as.num[is.na(geno)] = NA

	if( info ) { 
		return( list(allelecnt=geno.as.num, markers=markers[,c("N","numAlleles","input_ref","ref","alt")]) ) 
	} else {
		return(geno.as.num)
	}
}