#!/usr/bin/Rscript

# THIS PROGRAM DRAWS VSE PLOTS WITH MATRICES (intersection heatmaps)

args <- commandArgs(TRUE)
ID <- args[1]
PDFOUTPUT <- args[2]

library(car)

# Set normalization parameters:
rald    <- 1     # 1 raSNPs / 2 ldSNPs
norm    <- 1     # 1 yes / 0 no
toorder <- 1     # 1 yes / 0 no
fixp    <- 1     # Find best p. 1 fixed; else not fixed, i.e. search.
norm_p  <- 1  # If fixed == 1, the value of p. 0.75
Nall    <- 100   # Number of tests for Bonferroni.
null_size <- 100 # WAS 1000 ORIGINALLY


# Set graphing parameters:
load_names  <- TRUE
plot_matrix <- FALSE
plot_scale  <- FALSE

x <- read.table( ID, as.is = TRUE )
colnames(x) <- c( "AVS", sprintf( "%04d", 0:99), "BED") # WAS 0:999 ORIGINALLY

N <- dim(x)[1]
gray <- "gray30"
red  <- "red"

rav_all <- matrix(data = NA, nrow = 1, ncol = N)
null_all <- matrix(data = NA, nrow = null_size, ncol = N)
maxVal <- matrix(data = NA, nrow = 1, ncol=N)
for ( n in 1:N ){
  rav_all[,n] <-  x[ n, 1 ]
  null_all[,n] <-  t(x[ n, (1:null_size) + 1])
  maxVal[,n] <- max(null_all[,n])
  if (rav_all[,n] > maxVal[,n]){
    maxVal[,n] <- rav_all[,n]
  }
}
names <- x[["BED"]]
names <- gsub( "\\/.*\\/", "", names, perl = TRUE)
names <- sprintf( "%-15s", names)

pdf( paste( ID, "_density.pdf", sep = ""), height=3, width=3.5, pointsize=10)
for (n in 1:N){
  par(mar=c(4,4,2,1))
    #hist(null_all[,n], main=names[n], xlab="Overlaps", ylab="Frequency", xlim=c(0,maxVal[,n]+3))
    plot(density(null_all[,n]), main=names[n],xlab="No. of Overlaps", xlim=c(0,maxVal[,n]+3))
    points(rav_all[,n], 0, pch = 23, bg = "red", col = "red", lwd = 1)
 #   text(rav_all[,n], 5, labels = names[n],col = "red" )
}
dev.off()

# SD values of p-values

seq <- seq(1,5,by = 0.01)

pthresh <- -log10( 0.05 )                                # 1.30103
ptmp <- cbind( -log10(1 - pnorm(seq)) - pthresh, seq )
pthresh_sd <- ptmp[order( ptmp[,1] ^ 2),][1,2]           # 1.64

xthresh <- -log10( 1 / null_size )                       
xtmp <- cbind( -log10(1 - pnorm(seq)) - xthresh, seq )   
xthresh_sd <- xtmp[order( xtmp[,1] ^ 2),][1,2]           
                                                         
# The Bonferroni threshold                               
bthresh <- -log10( 0.05 / Nall )                         
btmp <- cbind( -log10(1 - pnorm(seq)) - bthresh, seq )   
bthresh_sd <- btmp[order( btmp[,1] ^ 2),][1,2]           

print( paste( "p-value threshold: ", bthresh, "-log10;", bthresh_sd, "SD. Based on", Nall, "comparisons." ))

# NORMALIZATION

pvals <- matrix(data = NA, nrow =         1, ncol = N)
 ravs <- matrix(data = NA, nrow =         1, ncol = N)
nulls <- matrix(data = NA, nrow = null_size, ncol = N)
nonor <- matrix(data = NA, nrow =         1, ncol = N)

print(paste("AVS","Enrichment","p-value","normalized_p-value","KST_p-value","-log10_normalize_p-value","sample"), sep="\t", quote = FALSE )
for (n in 1:N ){
    if( norm == 1 ){
        norm_p <- 1
        xran <- null_all[,n] + runif( null_size, 0, 1) # adding a uniform probability number to each null value
        x <- bcPower( xran,norm_p) #box-cox power transformation, lambda=1 (usually between -2 and 2)
        kst <- ks.test(x, "pnorm", mean = mean(x), sd = sd(x) , exact = TRUE) #one-sampled ks test
        kst_pvalue <- kst$p.value
        if( kst_pvalue < 0.05 ){
            scan <- matrix( NA, nrow = 0, ncol = 2 )
            xran <- null_all[,n] + runif( null_size, 0, 1)
            for (norm_p in seq( 0.5,1.5, by = 0.1)){
                x <- bcPower( xran,norm_p)
                kst <- ks.test(x, "pnorm", mean = mean(x), sd = sd(x) , exact = TRUE)
                scan <- rbind(scan, c( norm_p, kst$p.value ))
                print( paste( norm_p, kst$p.value, sep = " "))
            }
            scan <- scan[order(scan[,n], decreasing = TRUE),]
            norm_p <- scan[1,1]
            kst_pvalue <- scan[1,2]
        }
        avs  <- rav_all[,n]
        null <- null_all[,n]
        enrich <- avs / mean( null )
        perm_larger <- length( null[null >= avs] )
        perm_pvalue <- length( null[null >= avs] ) / null_size
        ravs[,n]  <- bcPower(  rav_all[,n] + 1,1)
        nulls[,n] <- bcPower( null_all[,n] + 1 ,1)

        # normalization
        nm <- median( nulls[,n] ) # null median
        ravs[,n]  <- ravs[,n] - nm
        nulls[,n] <- nulls[,n] - nm

        # scale <- quantile(nulls[,n])[5] 
        # null quantile. 4 is upper quartile, 5 is max value
        scale     <- sd( nulls[,n] )
        ravs[,n]  <- ravs[,n] / scale
        nulls[,n] <- nulls[,n] / scale

        norm_pvalue <- 1 - pnorm( ravs[,n] )
        pvals[,n] <- norm_pvalue
        nonor[,n] <- FALSE

        print( paste( 
            sprintf("%2d",   avs), 
            sprintf("%5.2f", enrich),
            sprintf("%6.4f", perm_pvalue), 
            sprintf("%4.2f", norm_p), 
            sprintf("%4.2f", kst_pvalue), 
            sprintf("%5.3f", -log10(norm_pvalue)), 
            names[n]),sep="\t",quote = FALSE)
    } else {
         ravs[,n] <- rav_all[,n]
        nulls[,n] <- null_all[,n]
        pvals[,n] <- length( nulls[,n][ nulls[,n] >= ravs[,n] ]) / null_size
           # the distribution is centered on either the mean or the median
          nm <- median( nulls[,n] )
         ravs[,n]  <- ravs[,n] - nm
        nulls[,n] <- nulls[,n] - nm

        # scale <- quantile(nulls[,n])[5] # null quantile. 4 is upper quartile, 5 is max value
          scale <- sd( nulls[,n] )
        if( scale >= 0.1){
             ravs[,n]  <- ravs[,n] / scale
            nulls[,n] <- nulls[,n] / scale
            nonor[,n] <- FALSE
        } else {
            max <- max( nulls[,n] )
            ravs[,n]  <- ravs[,n] / max * xthresh_sd
            nulls[,n] <- nulls[,n] / max * xthresh_sd
            nonor[,n] <- TRUE
        }
    }
}

# MRSVs DUNE PLOTS

rav_dune <- matrix(data = NA, nrow =         1, ncol = 0)
null_dune <- matrix(data = NA, nrow = null_size, ncol = 0)
nonor_dune <- matrix(data = NA, nrow =         1, ncol = 0)

# ordering
if( toorder == 1 ){
  on <- order( ravs, na.last = FALSE ) # order of null
} else {
  on  <- 1:N
}

names_tmp  <- paste( names, sprintf("%5.2f", -log10( pvals )))
names_dune <- names_tmp[on]
rav_dune <-     ravs[,on]
null_dune <-    nulls[,on]
nonor_dune <-    nonor[,on]

min <- min( c( nulls[!nulls == "NaN"] ) )
max <- max( c( ravs[!ravs == "NaN"] ) )

pdf(PDFOUTPUT)
mar.orig <- par()$mar # save the original values
par(mar = c(20,4,4,4)) # set your new values
boxplot(as.data.frame(null_dune), 
    horizontal = FALSE, las = 2, 
     ylim = c( min, max ), 
    family = "mono", las = 1,
    range = 0, names = names_dune, border = gray)
par(mar = mar.orig) # put the original values back

abline( h = pthresh_sd, col="grey" )
# abline( h = xthresh_sd )
abline( h = bthresh_sd, col="grey" )

color <- "green" # this should never appear
for (i in 1:N ){
    if(nonor_dune[i]){
        bg = "white"
        color = gray
    } else {
        if( is.finite(rav_dune[i]) ){
            bg = gray
            color = gray
        } else {
            bg = "green"
            color = "green"
        }
    }
    points(i,rav_dune[i], pch = 23, bg = bg, col = color, lwd = 1)
}
dev.off()
q( status = 0 )


