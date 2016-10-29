###########################################
# converts probabilities to logits
probs2logits <- function(probs){
	include_zero <- TRUE
	#*** initial values
	NP <- base::length(probs)
	q1 <- stats::qlogis(probs)
	x0 <- q1-q1[1]
	x0 <- x0[-1]
	#*** define optimization function
	fr <- function(x) {  
		NX <- base::length(x)  
		NP <- NX+1
		y <- base::rep(0,NP)      
		y[ base::seq(2,NX+1) ] <- x  
		p1 <- base::exp(y) / base::sum( base::exp(y) )
		p1 <- base::sum( (probs - p1)^2 )
		base::return(p1)
	}
	#*** optimization
	a1 <- stats::optim( x0 , fr)
	y0 <- a1$par
	if (include_zero){
		y0 <- base::c(0,y0)
	}
	base::return(y0)
}
##################################################			
logits2probs <- function(y){
	x <- base::exp(y)
	x <- x / base::sum(x)
	base::return(x)
}
##################################################
		