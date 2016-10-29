
#####################################################################
inits_theta_1dim <- function( dat.resp , pid , eps=.05 ){
	N <- base::length( base::unique( pid ) )
	# initial values
	maxK <- base::apply( dat.resp ,2 , base::max , na.rm=TRUE )
	I <- base::ncol(dat.resp)
	K <- base::max(maxK)
	# ability inits
	theta <- stats::aggregate( dat.resp , base::list(pid) , base::mean , na.rm=TRUE )[,-1] 
	theta <- theta / base::matrix( maxK , nrow=N , ncol=I , byrow=TRUE )
	theta <- base::rowMeans(theta , na.rm=TRUE)
	theta <- ( theta + eps ) / ( 1 + 2*eps )
	theta <- stats::qlogis( theta )
	base::return(theta)
}
#####################################################################		