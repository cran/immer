

immer_reshape_wideformat <- function( y , pid , rater , Nmin_ratings = 1 ){	
	y_dfr <- FALSE
	if ( base::is.data.frame(y) ){
		y_dfr <- TRUE			
	}
	if ( ! y_dfr ){
		dfr1 <- immer_reshape_wideformat_vector(y, pid,rater, Nmin_ratings )
	}
	if ( y_dfr ){
		NV <- base::ncol(y)		
		for (vv in 1:NV){		
			y1 <- base::as.vector(y[,vv])
			dfr2 <- immer_reshape_wideformat_vector(y=y1,  pid=pid,
							rater=rater, Nmin_ratings=Nmin_ratings )		
			base::colnames(dfr2)[-1] <- base::paste0( base::colnames(y)[vv] , 
											"_" , base::colnames(dfr2)[-1] ) 
			if ( vv == 1 ){ dfr1 <- dfr2 }
			if ( vv > 1 ){
				dfr1 <- base::merge( x = dfr1 , y = dfr2 , by = "pid" , all=TRUE )
			}
		}
	}	
	base::return(dfr1)
}
		

immer_reshape_wideformat_vector <- function(y, pid,rater, Nmin_ratings ){		
	rater <- base::paste(rater)
	Nobs <- base::rowsum( 1 - is.na(y) ,  pid  )
	Nobs <- Nobs[ Nobs >= Nmin_ratings , ]
	persons <- base::names(Nobs)
	NP <- base::length(persons)
	data <- base::data.frame( "pid" = pid , "rater" = rater , "y" = y )
	data <- data[ data$pid %in% persons , ]
	raters <- base::sort( base::unique( base::paste(data$rater )))
	RR <- base::length(raters)
	y <- base::matrix( NA , nrow=NP , ncol=RR+1 )
	y <- base::as.data.frame(y)
	colnames(y) <- c("pid" , raters )
	y$pid <- persons
	indM <- base::cbind( base::match( base::paste(data$pid) , persons) , 
				base::match( base::paste(data$rater) , raters)+1 ) 
	y[ indM ] <- data$y
	base::return(y)		
}