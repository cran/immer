
################################################
# simulation theta values
sampling_hrm_theta_1dim <- function( theta , N , I , maxK , a , b , xi , xi_ind , 
		dat , dat_ind , pid , MHprop , mu_theta , SD_theta , useRcpp , eps=1E-20 ){

	# refresh count
	MHprop$refresh_count$theta <- MHprop$refresh_count$theta + 1

	theta_new <- theta_old <- theta
	theta_new <- stats::rnorm(N , mean= theta , sd = MHprop$SD$theta )

	p_new <- base::log( stats::dnorm( theta_new , mean = mu_theta , sd = SD_theta ) )
	p_old <- base::log( stats::dnorm( theta_old , mean = mu_theta , sd = SD_theta ) )

	for (ii in 1:I){
		# ii <- 1    
		ll_new <- probs_gpcm( x = xi[,ii] , theta = theta_new , b = b[ii,] , a = a[ii] , K = maxK[ii] ,
					x_ind = xi_ind[,ii] , useRcpp )
		ll_new <- base::log( ll_new  + eps )               
		ll_old <- probs_gpcm( x = xi[,ii] , theta = theta_old , b = b[ii,] , a = a[ii] , K = maxK[ii] ,
					x_ind = xi_ind[,ii] , useRcpp )
		ll_old <- base::log( ll_old  + eps )
		p_new <- p_new + ll_new
		p_old <- p_old + ll_old
	}
					
	ratio <- base::exp( p_new - p_old )    
	zuf <- stats::runif(N)

	theta <- base::ifelse( ratio > zuf , theta_new , theta_old )
	MHprop$accept$theta <- MHprop$accept$theta + 1*( ratio > zuf )
	res <- base::list( theta = theta , MHprop = MHprop )
	base::return(res)
}
###################################################################	

		
		