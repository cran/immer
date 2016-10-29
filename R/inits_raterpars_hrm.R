
###########################################
# inits rater parameters
inits_raterpars_hrm <- function( rater , I , est_settings ){
	R <- base::length( base::unique(rater) )
	phi <- base::matrix( 0 , nrow=I , ncol=R )
	psi <- .3 + phi
	if ( est_settings$est.psi == "n" ){
		psi <- 1E-10 + 0*phi
	}
	res <- base::list(R=R , phi=phi , psi=psi)
	base::return(res)
}
#############################################			