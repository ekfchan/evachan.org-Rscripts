plot_marker_lox <- function(chr, lox) {

     ## Copyright 2006-2008 Eva Chan
     ## eva@evachan.org
     ##
     ## This funciton generates a visual representation of a set of markers
     ## onto the genome.
     ##
     ## Inputs
     ## chr: vector of chromosomes
     ## lox: numeric vector of markers' positions on the corresponsing chrs
     ## NOTE:: chr and lox are assumed to be in same marker order!!
     ##

     ## remove markers with missing chr or pos
     inds <- which( is.na(chr) | is.na(lox) )
     if(length(inds)>0) {
          warning(length(inds), " SNPs are missing map information; they are ignored.\n")
          chr <- chr[-inds]
          lox <- lox[-inds]
     }

     ## set non-integer chromosomes as integers
     unique.chrs <- unique(chr)
     suppressWarnings(unique.chrs.as.num <- as.integer(unique(unique.chrs)))
     non.int.chrs.ind <- which(is.na(unique.chrs.as.num))
     chr.as.num <- chr
     if(length(non.int.chrs.ind)>0) {
          num.int.chrs <- length(unique.chrs) - length(non.int.chrs.ind)
          for(i in 1:length(non.int.chrs.ind)) {
               unique.chrs.as.num[non.int.chrs.ind[i]] <- num.int.chrs + i
               chr.as.num[which(chr==unique.chrs[non.int.chrs.ind[i]])] <- num.int.chrs + i
          }
     }
     chr.as.num <- as.integer(chr.as.num)
     unique.chrs <- unique.chrs[order(unique.chrs.as.num)]
     unique.chrs.as.num <- unique.chrs.as.num[order(unique.chrs.as.num)]

     ## set lox to Mb if in bases
     new.lox <- lox / 1000000
     if( max(new.lox) > 1 ) {
         lox <- new.lox
         yunit <- "(Mb)"
         rm(new.lox)
     } else { yunit = "(bases)" }

     ## calculate chromosome range
     chr.len <- rep(NA, length(unique.chrs))
     for(i in 1:length(chr.len)) {
         chr.len[i] <- max(lox[which(chr==unique.chrs[i])],na.rm=T)
     }

     ## plot frame
     plot( unique.chrs.as.num, chr.len, ylim=c(0, max(chr.len)), pch="_", xlab="Chromosome", ylab=paste("position",yunit), las=1, axes=F, cex=1.2, col="dark grey" )
     points( unique.chrs.as.num, rep(0, length(unique.chrs.as.num)), pch="_", cex=1.2, col="dark grey" )
     axis(1, at=1:length(unique.chrs.as.num), labels=unique.chrs, las=1)
     axis(2, las=1)
     for(i in 1:length(chr.len)) {
         points( rep(unique.chrs.as.num[i],2), c(0,chr.len[i]), type="l" )
     }

     ## plot markers
     points( chr.as.num, lox, pch="_" )

}
