plot_markers_by_set <- function(set, chrom, pos, horiz=F, usr.colours=NULL, pt.cex=1) {

     ## Copyright 2008 Eva Chan 
     ##	eva@evachan.org
     ##
     ## Generates a map of markers with different sets of markers marked with different colours.
     ## Parameters:
     ## set: factor indicating the set to which the corresponding markers belong.
     ## chrom: vector of chromosome to which the markers are located; sex chromosomes (X,Y) allowed.
     ## pos: numeric vector of chromosome position of markers.
     ## horiz: logical indicating if chromosomes should be represented horizontally.
     ## usr.colours: colour vector for each marker set; if NULL, default rainbow colours are used;
     ##              if length is less than number of unique marker sets, colours will be recycled.
     ## pt.cex: expansion factor the markers and legend.
     ## NOTE:: Ordering of markers between the three parameters are assumed to be the same.

     ## check to ensure parameters are of the same length
     if( (length(chrom) != length(set)) || (length(pos) != length(set)) ) {
          stop("Parameters are of differnet length.\n")
     }

     ## exclude markers without known mapping info
     excl.inds <- unique(which(is.na(chrom) | is.na(pos)))
     if(length(excl.inds)>0) {
          warning(length(excl.inds)," markers have missing position and are excluded in plot.\n")
          set <- set[-excl.inds]
          chrom <- chrom[-excl.inds]
          pos <- pos[-excl.inds]
     }

     ## set "set" as factor variable
     set <- as.factor(as.character(set))
     
     ## check to ensure positions are given in numeric variables
     if(!is.numeric(pos)) {
          if(sum(is.na(suppressWarnings(as.numeric(pos))))>0) {
               stop("Inappropriate chromosome positions.\n")
          }
          pos <- as.numeric(pos)
     }
     
     ## check for and recode sex chromosomes
     chrom.ori <- chrom
     chrom.as.num <- suppressWarnings(as.numeric(chrom))
     if(sum(is.na(chrom.as.num)) > 0) {
          ## check for chrom X
          inds <- grep( "X", chrom, ignore.case=T )
          if(length(inds)>0) {
               chrom[inds] <- sum(!is.na(unique(chrom.as.num))) + 1
               chrom.as.num <- suppressWarnings(as.numeric(chrom))
          }
     }
     if(sum(is.na(chrom.as.num)) > 0) {
          ## check for chrom Y
          inds <- which( (chrom=="Y") | (chrom=="y") )
          if(length(inds)>0) {
               chrom[inds] <- sum(!is.na(unique(chrom.as.num))) + 1
               chrom.as.num <- suppressWarnings(as.numeric(chrom))
          }
     }
     if(sum(is.na(chrom.as.num)) > 0) {
          stop("Inappropriate chromosome:\n", paste(unique(chrom[which(is.na(chrom.as.num))]),collapse=", "), "\n")
     }

     ## set chromosome length
     pos.lab=""
     if( median(z <- (pos/1000000)) > 1 ) {
          pos <- z
          pos.lab <- "(Mb)"
     } else {
          if( median(z <- (pos/1000)) > 1 ) {
               pos <- z
               pos.lab <- "(kb)"
          }
     }

     unique.chroms <- sort(unique(chrom.as.num))
     if(is.null(usr.colours)) {
          usr.colours <- rainbow(nlevels(set))
     }
     if(length(usr.colours) < nlevels(set)) {
          warning("Colours provided is fewer than marker sets-- colours will be recycled\n")
          usr.colours <- rep(usr.colours,len=nlevels(set))
     }
     if(horiz) {
          plot( c(0,max(pos[which(chrom.as.num==unique.chroms[1])])), c(unique.chroms[1],unique.chroms[1]), type="l", ylim=c(0,length(unique.chroms)), xlim=c(0,ceiling(1.1*max(pos))), ylab="Chromosomes", xlab=paste("Position",pos.lab), las=1, axes=F )
          axis.lab <- unique(chrom.ori)[match(unique.chroms,unique(chrom.ori))]
          if(sum(is.na(axis.lab))==1) {
               axis.lab[(length(axis.lab))] <- "X"
          }
          if(sum(is.na(axis.lab))==2) {
               axis.lab[(length(axis.lab)-1):length(axis.lab)] <- c("X","Y")
          }
          axis(2, at=unique.chroms, labels=axis.lab, las=1 )
          axis(1, las=1)
          for(i in 2:length(unique.chroms)) {
               points( c(0,max(pos[which(chrom.as.num==unique.chroms[i])])), c(unique.chroms[i],unique.chroms[i]), type="l" )
          }
          for(i in 1:nlevels(set)) {
               text( pos[which(set==levels(set)[i])], chrom.as.num[which(set==levels(set)[i])], labels="|", col=usr.colours[i], cex=pt.cex )
          }
          legend(max(pos),length(unique.chroms),legend=levels(set), col=usr.colours, pch=45, bty="o", horiz=F, pt.cex=pt.cex, ncol=1)
     } else {
          plot( c(unique.chroms[1],unique.chroms[1]), c(0,max(pos[which(chrom.as.num==unique.chroms[1])])), type="l", xlim=c(0,length(unique.chroms)), ylim=c(0,ceiling(1.1*max(pos))), xlab="Chromosomes", ylab=paste("Position",pos.lab), las=1, axes=F )
          axis.lab <- unique(chrom.ori)[match(unique.chroms,unique(chrom.ori))]
          if(sum(is.na(axis.lab))==1) {
               axis.lab[(length(axis.lab))] <- "X"
          } 
          if(sum(is.na(axis.lab))==2) {
               axis.lab[(length(axis.lab)-1):length(axis.lab)] <- c("X","Y")
          }
          axis(1, at=unique.chroms, labels=axis.lab, las=1 )
          axis(2, las=1)
          for(i in 2:length(unique.chroms)) {
               points( c(unique.chroms[i],unique.chroms[i]), c(0,max(pos[which(chrom.as.num==unique.chroms[i])])), type="l" )
          }
          for(i in 1:nlevels(set)) {
               text( chrom.as.num[which(set==levels(set)[i])], pos[which(set==levels(set)[i])], labels="--", col=usr.colours[i], cex=pt.cex )
          }
          legend(2,max(pos),legend=levels(set), col=usr.colours, pch=45, bty="o", horiz=T, pt.cex=pt.cex)
     }

}
