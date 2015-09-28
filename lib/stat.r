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

  rav_all <- matrix(data = NA, nrow =         1, ncol = N)
 null_all <- matrix(data = NA, nrow = null_size, ncol = N)

for ( n in 1:N ){
         rav_all[,n] <-  x[ n, 1 ]
        null_all[,n] <-  t(x[ n, (1:null_size) + 1])
}                            





   names <- x[["BED"]]
   names <- gsub( "\\/.*\\/", "", names, perl = TRUE)
   #names <- gsub( "\\..*", "", names, perl = TRUE)
   #names <- gsub( "\\.bed", "", names, perl = TRUE)
   #names <- gsub( "Human_", "", names, perl = TRUE)
   #names <- gsub( "_", " ", names, perl = TRUE)
   names <- sprintf( "%-15s", names)


if( plot_matrix ){
    matrix_file <- paste( ID, "_matrix.txt", sep = "")
}







#if(FALSE){
    
    # this chunk of code generates a single histogram with a null distribution
    n <- 1
    pdf( paste( ID, "_single.pdf", sep = ""), height = 8, width = 6)
    hist(null_all[,n], xlim = c(0,20))
    points(rav_all[,n], 0, pch = 23, bg = "red", col = "red", lwd = 1)
    text(rav_all[,n], 25, labels = names[n],col = "red" )
    dev.off()

    pdf( paste( ID, "_boxplot.pdf", sep = ""), height = 6, width = 2)
    boxplot(null_all[,n], ylim = c(0,20))
    points( 1 ,rav_all[,n], pch = 23, bg = "red", col = "red", lwd = 1)
    dev.off()
#}




# Fix N for all nine panels

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

# print( paste( "p-value threshold: ", bthresh, "-log10;", bthresh_sd, "SD. Based on", Nall, "comparisons." ))




# NORMALIZATION


# the larger the exponent the larger the deformation
# for FoxA1, going from 0.4 to 1.2 changes dramatically
# there's some randoness here... probably in the p-values
# do I get consistency if I use the function?
# I'm using a random normal function for the ks.test...


pvals <- matrix(data = NA, nrow =         1, ncol = N)
 ravs <- matrix(data = NA, nrow =         1, ncol = N)
nulls <- matrix(data = NA, nrow = null_size, ncol = N)
nonor <- matrix(data = NA, nrow =         1, ncol = N)


# print( " A     E    >     PP    p  KSP    NP BED", quote = FALSE )
for (n in 1:N ){

    if( norm == 1 ){
        norm_p <- 1
        xran <- null_all[,n] + runif( null_size, 0, 1)
        x <- bcPower( xran,1)
        kst <- ks.test(x, "pnorm", mean = mean(x), sd = sd(x) , exact = TRUE)
        kst_pvalue <- kst$p.value

        if( kst_pvalue < 0.05 ){

            scan <- matrix( NA, nrow = 0, ncol = 2 )
            xran <- null_all[,n] + runif( null_size, 0, 1)
            for (norm_p in seq( 0.5,1.5, by = 0.1)){
                x <- bcPower( xran,1)
                kst <- ks.test(x, "pnorm", mean = mean(x), sd = sd(x) , exact = TRUE)
                scan <- rbind(scan, c( norm_p, kst$p.value ))
                # print( paste( p, kst$p.value, sep = " "))
            }
            scan <- scan[order(scan[,2], decreasing = TRUE),]
            norm_p <- scan[1,1]
            kst_pvalue <- scan[1,2]
        }

        avs  <- rav_all[,n]
        null <- null_all[,n]
        enrich <- avs / mean( null )
        perm_larger <- length( null[null >= avs] )
        perm_pvalue <- length( null[null >= avs] ) / null_size

        # if( kst_pvalue >= 0.05 ){
        # I'm no longer using this as I also have the perm_pvalue.
    
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
            sprintf("%4d",   perm_larger), 
            sprintf("%6.4f", perm_pvalue), 
            sprintf("%4.2f", norm_p), 
            sprintf("%4.2f", kst_pvalue), 
            sprintf("%5.3f", -log10(norm_pvalue)), 
            names[n]),quote = FALSE)

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

if( FALSE){
    # For Batool
    names2 <- c("G9", "G8", "G7", "G6", "G5", "G4", "G3", "G2", "G1", "CO", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9")
    colnames(ravs) <- names2
    colnames(nulls) <- names2
    rownames(ravs) <- "CRC"
    rownames(nulls) <- 1:1000
    write.table( rbind(ravs,nulls), file = "VEL_VSE.txt",
        quote = FALSE, sep = "\t",
        row.names = TRUE, col.names = TRUE )
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


out_pdf <- paste( ID, "_", null_size,"_rald", rald, "_norm", norm, "_fixp", fixp, "_toorder", toorder, ".VSE.pdf", sep = "" )

pdf( out_pdf, width= N/4 + 1, height= 10, family = "mono")
# par(fig = c( 0.1,1.0,0.0,1.0 ), mfrow = c(3,1))
par( mfrow = c(2,1))



min <- min( c( nulls[!nulls == "NaN"] ) )
max <- max( c( ravs[!ravs == "NaN"] ) )

pdf(PDFOUTPUT)
mar.orig <- par()$mar # save the original values
par(mar = c(20,4,4,4)) # set your new values
boxplot(as.data.frame(null_dune), 
    horizontal = FALSE, las = 2, 
#    ylim = c( -2, 10 ), 
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





if( FALSE ){
    options(width = 275)
    ID <- "Breast_core"
    matrix_file <- paste( ID, "_matrix.txt", sep = "")
    order_file <- "../NHGRI_Breast_cancer_order.txt"
}



if( plot_matrix ){
    jpeg("heatmapmatrix.jpg")
    print("IN PLOT MATRIX")    
    # Reorder the raSNPs:
    # The following is intended to allow an arbitrary ordering.
    # matrix_file <- "lost_matrix.txt"
    M <- read.table( matrix_file )
    M[M > 1] <- 1

    order_genome <- read.table( order_file )
    colnames(order_genome) <- c("genome", "rsid")
    
    # Store the current order (wrong or arbitrary).
    order_wrong <- data.frame( wrong =  1:dim(M)[2], rsid = colnames(M))

    orders <- merge( order_wrong, order_genome) # Merge both orders.
    orders <- orders[ order( orders$wrong ),]   # Make sure order is the original (wrong).
    M <- M[ ,order( orders$genome )]            # Order M with the genomic order.

    # Set the colors:
    ncols <- 100
    value <- 1.0
    # boy_colors  <- hsv(h = 6/6, s = seq.int( 0, 0.9, length.out = ncols)%%1, v = value)
    gray_colors  <- gray( seq.int( 1, 0.5 , length.out = ncols) )

    # Plot the heatmap:
    # heatmap( log2(as.matrix( t(M[on,] + 1) )), 

    
    heatmap( as.matrix( t(M[on,] + 1) ), 
       col = gray_colors,
       Rowv = NA, 
       Colv = NA , 
       scale = "none"

    )
    dev.off()
}




if( plot_scale ){
    print("IN PLOT MATRIX SCALE")
    scale <- seq.int(min(M), max(M), length.out = 10)
    scaleM <- cbind( scale, scale, scale)

    jpeg("heatmaplogmatrix.jpg")
    heatmap( log2(as.matrix( t( scaleM + 1) )), 
       col = gray_colors,
       Rowv = NA, 
       Colv = NA , 
       scale = "none",
       labCol = round(scale, digits = 0)
    )
    dev.off()
}


q( status = 0 )


