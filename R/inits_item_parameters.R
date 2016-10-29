
################################################
# inits item parameters
inits_itempars <- function( dat , prior ){
	maxK <- base::apply( dat , 2 , base::max , na.rm=TRUE )
	K <- base::max(maxK)
    I <- base::ncol(dat)
    b <- base::matrix( NA , nrow=I , ncol=K)
	b1 <- base::colMeans( dat , na.rm=TRUE ) / maxK
	for (ii in 1:I){
		# ii <- 1
		b[ii , base::seq(1,maxK[ii],1) ] <- b1[ii] + base::seq( -2 , 2 , length= maxK[ii]  )
	}
					
	if ( ! is.null( prior$b$M ) ){					
		b <- prior$b$M
	}
	if ( ! is.null( prior$a$M ) ){					
		a <- prior$a$M
	}					
	a <- base::rep(1,I)
	res <- base::list( b=b , maxK = maxK , K = K , a = a , I=I)
	base::return(res)
}
#################################################				