

###****************
# probs GPCM
# INPUT:
#  x ... vector of categories
#  theta ... vector of abilities
#  b ... vector of item parameters
#  a ... item discrimination
#  K ... maximum category (1,2,...)
probs_gpcm <- function( x , theta , b , a , K , x_ind = NULL , useRcpp=FALSE){
	if ( ! useRcpp ){
		probs <- probs_gpcm_R( x , theta , b , a , K , x_ind  )
	}
	if ( useRcpp ){
	    if ( base::is.null( x_ind ) ){ 
			x_ind <- base::rep(1, base::length(theta)) 
		}
		probs <- base::.Call( "probs_gpcm_rcpp" , x , theta , b , a , K , x_ind  ,
							PACKAGE="immer")
	}															
    base::return(probs)
}
##################################################################			
			
#*****************
# R version
probs_gpcm_R <- function( x , theta , b , a , K , x_ind = NULL ){
    N <- base::length(theta)
	KM <- base::matrix( 0:K , nrow = N , ncol= K+1 , byrow=TRUE)
    b0 <- base::c( 0 , b[1:K] )
	bM <- base::matrix( b0 , nrow = N , ncol= K+1 , byrow=TRUE)
    probs <- base::exp( a * KM *  theta - bM )
    probs <- probs / base::rowSums(probs , na.rm=TRUE)
	if ( ! base::is.null(x) ){
		ind <- base::cbind( 1:N , x+1)
		probs <- probs[ ind ]
	}
	if ( ! base::is.null( x_ind) ){
		probs <- base::ifelse( x_ind == 0 , 1 , probs )
	}								
    base::return(probs)
}