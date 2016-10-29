
#############################################################
# simulating hierarchical rater model
simulate_HRM <- function( theta , a , b , phi , psi ){

	RR <- base::ncol(phi)
    I <- base::nrow(b)
	#*** simulate partial credit items
	K <- base::ncol(b)
	N <- base::length(theta)
	oneN <- base::rep(1,N)
	KM <- base::outer( oneN , 0:K )
	xiM <- base::matrix(0,nrow=N, ncol=I)
	for (ii in 1:I){
		# ii <- 1
#		probs1 <- a[ii] * theta * ( 1:K ) - b[ ii , ]
	    KM <- base::matrix( 0:K , nrow = N , ncol= K+1 , byrow=TRUE)
        b0 <- c( 0 , b[ii , 1:K] )
	    bM <- base::matrix( b0 , nrow = N , ncol= K+1 , byrow=TRUE)
        probs <- base::exp( a * KM *  theta - bM )					
		probs <- probs / base::rowSums( probs )	
		probs <- sirt::rowCumsums.sirt(probs)
		vals <- sirt::rowIntervalIndex.sirt( matr = probs, rn = stats::runif(N) )
		xiM[,ii] <- vals - 1
	}
	#*** simulate items for all raters
	items <- base::paste0("I" , 1:I)
	dat <- NULL
	for (rr in 1:RR){
		# rr <- 1
		dat.rr <- base::matrix( NA , nrow=N , ncol=I)
		base::colnames(dat.rr) <- items
		for (ii in 1:I){
			# ii <- 1
			probs <- base::exp( - ( KM - ( xiM[,ii] + phi[ii,rr] ) )^2 / psi[ii,rr] / 2 )
			probs <- probs / base::rowSums(probs )
			probs <- sirt::rowCumsums.sirt(probs)
			vals <- sirt::rowIntervalIndex.sirt(matr=probs, rn= stats::runif(N))
			dat.rr[,ii] <- vals - 1
		}
		dat <- base::rbind( dat , dat.rr )
	}
	dat1 <- base::data.frame( "pid" = base::rep(1:N, RR) , 
								"rater" = base::rep(1:RR , each=N) , dat )
	dat1 <- dat1[ base::order( dat1$pid ) , ]
	base::rownames(dat1) <- NULL
	base::return(dat1)
}