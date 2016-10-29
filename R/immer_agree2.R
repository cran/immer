
###########################################
# agreement statistics for two raters
immer_agree2 <- function( y , w = rep(1,nrow(y) ) ,
		symmetrize = FALSE , tol= c(0,1) ){		
	CALL <- base::match.call()
	res <- immer_unique_patterns( dat = y , w=w )
	y <- res$y
	w <- res$w
	#***
	# symmetrize frequency table
	if ( symmetrize ){
		y <- base::rbind( y , y[,c(2,1) ] )
		w <- base::c( w/2 , w/2 )
		res <- immer_unique_patterns( dat = y , w=w )
		y <- res$y
		w <- res$w    
	}
	w0 <- w / base::sum(w)

	#*** create frequency table
	categs <- base::unique( base::c( y[,1] , y[,2] ) )
	CC <- base::length(categs)

	agree_table <- base::matrix( 0 , nrow=CC , ncol=CC)
	for (ii in 1:CC){
		for (jj in 1:CC){
			agree_table[ii,jj] <- base::sum( ( y[,1] == categs[ii] ) * 
					( y[,2] == categs[jj] ) * w0 )
		}
	}
	
	base::rownames(agree_table) <- base::paste0( base::colnames(y)[1] , "_Cat" , categs)
	base::colnames(agree_table) <- base::paste0( base::colnames(y)[2] , "_Cat" , categs)

	#****
	# compute absolute agreement (with tolerance specified in vector)
	TT <- base::length(tol)
	agree_raw <- base::rep(NA,TT)
	base::names(agree_raw) <- base::paste0("tol" , tol)
	for (tt in 1:TT){
		# tt <- 1
		agree_raw[tt] <- base::sum( w0[ base::abs( y[,1] - y[,2] ) <= tol[tt] ]  )
					}

	#*** marginal probabilities
	marg <- base::matrix( 0 , nrow=3 , ncol=CC)
	base::rownames(marg) <- base::c( base::colnames(y) , "Mean" )
	base::colnames(marg) <- base::paste0( "Cat" , categs )
	marg[1,] <- base::rowSums( agree_table )
	marg[2,] <- base::colSums( agree_table )
	marg[3,] <- ( marg[1,] + marg[2,] ) / 2

	#*** overall agreement
	Pa <- base::sum( base::diag( agree_table ) )

	#**** chance agreement
	Pe <- base::c()
	# Scott's Pi
	Pe["pi"] <- base::sum( marg[3,]^2 )
	# Cohen's kappa
	Pe["kappa"] <- base::sum( marg[1,] * marg[2,] )
	# AC1 Gwet
	Pe["AC1"] <- base::sum( marg[3,] * ( 1 - marg[3,] ) ) / ( CC - 1 )

	#*** algorithm Aicken's alpha
	PAk <- marg[1,]
	PBk <- marg[2,]
	res0 <- agree_aicken( PAk=PAk , PBk=PBk , Pa=Pa )

	Pe["Aicken"] <- res0$Pe
	agree_stats <- ( Pa - Pe ) / ( 1 - Pe )
	agree_stats["Aicken"] <- res0$alpha
	PAH <- res0$PAH
	PBH <- res0$PBH
	# hard to classify probs
	PH <- base::matrix( 0 , nrow=2 , ncol=CC )
	base::colnames(PH) <- base::colnames(marg)
	base::rownames(PH) <- base::rownames(marg)[1:2]
	PH[1,] <- PAH
	PH[2,] <- PBH

	#*****
	agree_stats <- ( Pa - Pe ) / ( 1 - Pe )
		
	#-----
	# output
	res <- base::list( "agree_raw" = agree_raw ,
			"agree_stats" = agree_stats ,
			"agree_table" = agree_table ,
			"marg" = marg , 
			"Pe" = Pe , "Pa" = Pa , 
			"alpha" = res0$alpha	, "PH" = PH,
			"nobs" = base::sum(w) , "ncat" = CC , "tol" = tol , 
			"y" = y , "w" = w , 
			"CALL" = CALL
			)
	base::class(res) <- "immer_agree2"
	base::return(res)
}
#############################################################