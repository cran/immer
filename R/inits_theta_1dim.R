
#####################################################################
inits_theta_1dim <- function( dat.resp , pid , eps=.05 ){
	N <- length( unique( pid ) )
	# initial values
	maxK <- apply( dat.resp ,2 , max , na.rm=TRUE )
	I <- ncol(dat.resp)
	K <- max(maxK)
	# ability inits
	theta <- aggregate( dat.resp , list(pid) , mean , na.rm=TRUE )[,-1] 
	theta <- theta / matrix( maxK , nrow=N , ncol=I , byrow=TRUE )
	theta <- rowMeans(theta , na.rm=TRUE)
	theta <- ( theta + eps ) / ( 1 + 2*eps )
	theta <- qlogis( theta )
	return(theta)
		}
#####################################################################		