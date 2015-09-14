
#####################################
# truncated normal distribution
rtrnorm <- function( N , mean , sd , lower = rep(-Inf,N) , upper=rep(Inf,N) ){
	t1 <- pnorm( lower , mean = mean , sd = sd )
	t2 <- pnorm( upper , mean = mean , sd = sd )
	rn <- runif( N , t1 , t2 )
	qnorm( rn , mean = mean , sd = sd )
				}
#########################################				