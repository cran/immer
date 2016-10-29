


###########################################################
# likelihood
IRT.likelihood.immer_HRM <- function( object , ... ){
    ll <- object$f.yi.qk
    base::attr(ll,"theta") <- object$theta_like
	base::attr(ll,"prob.theta") <- object$pi.k
	base::attr(ll,"G") <- 1
    base::return(ll)
}
#############################################################		


###########################################################
# posterior
IRT.posterior.immer_HRM <- function( object , ... ){
    ll <- object$f.qk.yi
    base::attr(ll,"theta") <- object$theta_like
	base::attr(ll,"prob.theta") <- object$pi.k
	base::attr(ll,"G") <- 1
    base::return(ll)
}
#############################################################	