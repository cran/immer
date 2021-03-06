## File Name: immer_opcat.R
## File Version: 0.15

#######################################################
# computing integer item discriminations
# given harmonic mean and minimum and maximum discrimination
immer_opcat <- function(a, hmean, min=1, max=10, maxiter=200 )
{
    # compute init factor
    fac <- hmean / harm_mean(a)
    g0 <- harm_mean(round_squeeze( a * fac, digits=0, min=min, max=max )    )
    fac1 <- fac/4
    g1 <- harm_mean( round_squeeze( a * fac1, digits=0, min=min, max=max ))
    fac2 <- 4*fac
    g2 <- harm_mean( round_squeeze( a * fac2, digits=0, min=min, max=max ))
    iter <- 0
    aint <- rep( 1E5, length(a))
    conv <- FALSE
    #********
    # iterations
    while( ! conv ){
        aint_old <- aint
        if ( g0 > hmean ){
            fac2 <- fac
            g2 <- g0
        } else {
            fac1 <- fac
            g1 <- g0
        }

        fac <- ( fac1 + fac2 ) / 2
        aint <- round_squeeze( a * fac, digits=0, min=min, max=max )
        g0 <- harm_mean( aint    )
        change <- max( abs( aint - aint_old ))
        if ( change==0 ){ conv <- TRUE }
        if ( iter==maxiter ){ conv <- TRUE }
        iter <- iter + 1
    }
    return(aint)
}
#**************************************************

