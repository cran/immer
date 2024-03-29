## File Name: IRT.likelihood.immer_hrm.R
## File Version: 0.092



###########################################################
# likelihood
IRT.likelihood.immer_hrm <- function( object, ... )
{
    ll <- object$f.yi.qk
    attr(ll,'theta') <- object$theta_like
    attr(ll,'prob.theta') <- object$pi.k
    attr(ll,'G') <- 1
    return(ll)
}
#############################################################


###########################################################
# posterior
IRT.posterior.immer_hrm <- function( object, ... )
{
    ll <- object$f.qk.yi
    attr(ll,'theta') <- object$theta_like
    attr(ll,'prob.theta') <- object$pi.k
    attr(ll,'G') <- 1
    return(ll)
}
#############################################################
