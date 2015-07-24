gwas_lm <- function(pheno, geno, model = NULL) {
### Copyright 2006 - 2008 Eva Chan
### http://evachan.org
### Created: August 2006
### Last modified: July 2014
###
### For each of the traits in pheno, perform linear regression on one or more 
### allelic models, using the genotype data provided in geno
###
### *** inputs ***
### pheno: matrix, or data.frame, of trait values: one trait per column with 
###        rownames(pheno) being sample IDs
### geno: matrix of genotypes: {0,1,2,NA}: one marker per row with rownames(geno) being 
###       marker IDs and one individual per column with colnames(geno) being sample IDs
### model: character vector listing the inheritence models to be tested;
###        avaiable options are: logadditive, codominance; dominance, overdominance, recessive,
###        or all (which is all of these five models): first three letter matches
###
### *** output ***
### res$
###     trait$
###           model
### Outputs a list of traits of lists of models of m-by-6 matrix,
### where rows of matrix corrrespond to each marke and
### columns of matrix to correspond to "f.stat", "df1", "df2", "p.val", "r.sq", "adj.r.sq"

     stopifnot( nrow(pheno) == ncol(geno) )

     ## Get phenotype names
     if(is.null(colnames(pheno))) { traits <- paste("trait",1:ncol(pheno),sep="") } else { traits <- colnames(pheno) }
     num.traits <- length(traits)

     ## Get marker names
     if(is.null(rownames(geno))) { markers <- paste("marker",1:nrow(geno),sep="") } else { markers <- rownames(geno) }
     num.markers <- length(markers)

     ## get sample names
     pheno.samples <- rownames(pheno)
     geno.samples <- colnames(geno)
     if(is.null(pheno.samples) | is.null(pheno.samples)) {
          if(nrow(pheno) == ncol(geno)) { 
               samples <- paste("sample",1:nrow(pheno),sep="") 
          } else {
               stop("Different sample numbers in pheno and geno!\n") 
          }
     } else {
          if(!all(is.element(pheno.samples, geno.samples)) ) { 
               stop("Different samples in pheno and geno!\n") 
          } else {
               if(!all(pheno.samples == geno.samples) ) {   #samples in different order between geno & pheno
                    pheno <- pheno[match(geno.samples, pheno.samples),]
               }
               samples <- geno.samples
          }
     }
     rm(geno.samples, pheno.samples)

     ## Determine genetic models to test
     do.codom = do.logadd = do.dom = do.rec = do.overdom = F
     if(is.null(model) | is.na(model) | model=="") stop("No inheritance model selected\n")
     if(is.element("all", model)) { do.codom = T; do.logadd = T; do.dom = T; do.rec = T; do.overdom = T }
     submodel <- substr(model,1,3)
     if(is.element('add',submodel) | is.element('log',submodel)) { do.logadd=T }
     if(is.element('cod',submodel)) { do.codom=T }
     if(is.element('dom',submodel)) { do.dom=T }
     if(is.element('rec',submodel)) { do.rec=T }
     if(is.element('ove',submodel)) { do.overdom=T }

     res <- list()

     ## Perform GWAS on each trait
     for(i in 1:num.traits)
     {
          cat(traits[i],"\n")
          cur.phval <- as.vector(pheno[,i])
          inds <- which(!is.na(cur.phval))
          if (length(inds)<3) { next }  #skip trait if fewer than 2 datapoints
          cur.phval <- cur.phval[inds]

          if(do.codom) codom.mat <- matrix(NA, ncol=6, nrow=num.markers, dimnames=list(markers,c("f.stat", "df1", "df2", "p.val", "r.sq", "adj.r.sq")))
          if(do.logadd) logadd.mat <- matrix(NA, ncol=6, nrow=num.markers, dimnames=list(markers,c("f.stat", "df1", "df2", "p.val", "r.sq", "adj.r.sq")))
          if(do.dom) dom.mat <- matrix(NA, ncol=6, nrow=num.markers, dimnames=list(markers,c("f.stat", "df1", "df2", "p.val", "r.sq", "adj.r.sq")))
          if(do.overdom) overdom.mat <- matrix(NA, ncol=6, nrow=num.markers, dimnames=list(markers,c("f.stat", "df1", "df2", "p.val", "r.sq", "adj.r.sq")))
          if(do.rec) rec.mat <- matrix(NA, ncol=6, nrow=num.markers, dimnames=list(markers,c("f.stat", "df1", "df2", "p.val", "r.sq", "adj.r.sq")))

          for(j in 1:num.markers)
          {
               cur.geval <- geno[j,inds]
               if(nlevels(as.factor(cur.geval))<=1) {next}

               if(do.codom)
               {
                    z <- summary( lm( cur.phval ~ as.factor(cur.geval) ) )
                    codom.mat[j,] <- c( z$fstatistic[1], z$fstatistic[2], z$fstatistic[3],
                                    1-pf(z$fstatistic[1],z$fstatistic[2],z$fstatistic[3]),
                                    z$r.squared, z$adj.r.squared )
               }

               if(do.logadd)
               {
                    z <- summary( lm( cur.phval ~ as.numeric(cur.geval) ) )
                    logadd.mat[j,] <- c( z$fstatistic[1], z$fstatistic[2], z$fstatistic[3],
                                     1-pf(z$fstatistic[1],z$fstatistic[2],z$fstatistic[3]),
                                     z$r.squared, z$adj.r.squared )
               }

               if(do.dom)
               {
                    g <- as.factor(cur.geval == 0)
                    if(nlevels(g)>1)
                    {
                         z <- summary( lm( cur.phval ~ g ) )
                         dom.mat[j,] <- c( z$fstatistic[1], z$fstatistic[2], z$fstatistic[3],
                                         1-pf(z$fstatistic[1],z$fstatistic[2],z$fstatistic[3]),
                                         z$r.squared, z$adj.r.squared )
                    }
               }

               if(do.rec)
               {
                    g <- as.factor(cur.geval == 2)
                    if(nlevels(g)>1)
                    {
                         z <- summary( lm( cur.phval ~ g ) )
                         rec.mat[j,] <- c( z$fstatistic[1], z$fstatistic[2], z$fstatistic[3],
                                         1-pf(z$fstatistic[1],z$fstatistic[2],z$fstatistic[3]),
                                         z$r.squared, z$adj.r.squared )
                    }
               }

               if(do.overdom)
               {
                    g <- as.factor(cur.geval == 1)
                    if(nlevels(g)>1)
                    {
                         z <- summary( lm( cur.phval ~ g ) )
                         overdom.mat[j,] <- c( z$fstatistic[1], z$fstatistic[2], z$fstatistic[3],
                                         1-pf(z$fstatistic[1],z$fstatistic[2],z$fstatistic[3]),
                                         z$r.squared, z$adj.r.squared )
                    }
               }
          }

          res[[traits[i]]] <- list()
          if(do.codom) res[[traits[i]]][["codom"]] <- codom.mat
          if(do.logadd) res[[traits[i]]][["logadd"]] <- logadd.mat
          if(do.dom) res[[traits[i]]][["dom"]] <- dom.mat
          if(do.rec) res[[traits[i]]][["rec"]] <- rec.mat
          if(do.overdom) res[[traits[i]]][["overdom"]] <- overdom.mat
     }

     res

}