## File Name: prior_hrm.R
## File Version: 0.212

##########################################################
# prior distributions for HRM
prior_hrm <- function( prior, b, a, phi,est_settings, sd_init=1 )
{
    prior0 <- prior
    est.a <- est_settings$est.a
    est.sigma <- est_settings$est.sigma
    est.mu <- est_settings$est.mu
    est.phi <- est_settings$est.phi
    est.psi <- est_settings$est.psi
    prior <- list()
    psi <- phi
    # b and a parameters
    prior$b$M <- 0*b
    prior$b$SD <- 10 + 0*b
    # prior$a$M <- 0*a
    # prior$a$SD <- 10 + 0*a
    prior$a$M <- 0*a
    if ( ! est.a ){
        prior$a$SD <- 1E-6 + 0*a
    } else {
        prior$a$SD <- 5 + 0*a
    }

    # phi and psi parameters
    prior$phi$M <- 0 + 0*phi
    prior$psi$M <- .4 + 0*psi
    prior$phi$SD <- 10 + 0*phi
    prior$psi$SD <- 10 + 0*psi
    if ( est.psi=='n' ){
        prior$psi$M <- 1E-10+ 0*psi
        prior$psi$SD <- 1E-30 +0*psi
    }

    # distribution parameters
    prior$mu$M <- 0
    if ( ! est.mu ){
        prior$mu$SD <- 1E-10
    } else {
        prior$mu$SD <- 1E5
    }
    if ( est.sigma ){
        prior$sigma2$w0 <- .001
    } else {
        prior$sigma2$w0 <- 1E10
    }
    prior$sigma2$sig02 <- sd_init

    #****************
    # defaults from prior
    if ( ! is.null( prior0 )){
        L1 <- length(prior0)
        for (vv in 1L:L1){
            # vv <- 1
            vv_label <- names(prior0)[vv]
            prior0.vv <- prior0[[vv]]
            L2 <- length( prior0[[vv]])
            for (zz in 1L:L2){
                zz_label <- names(prior0.vv)[[zz]]
                prior[[ vv_label ]][[ zz_label ]] <- prior0.vv[[zz]]
            }
        }
    }
    return(prior)
}
##########################################################
