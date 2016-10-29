
####################################################################
immer_identifiers_relabel <- function( dat , pid , rater ){

			pid <- base::paste( pid )
			rater <- base::paste( rater )

			dat0 <- base::data.frame( pid , rater , dat )
			pid_unique <- base::unique( base::paste(pid) )
			rater_unique <-  base::unique( base::paste( rater ) )			
			dat0$pid <- base::match( pid , pid_unique )
			dat0$rater <- base::match( rater , rater_unique )
			
			# item parameters
			item0 <- base::data.frame( "item" = base::colnames(dat) , 
						"N_Rat" = base::colSums( 1 - base::is.na( dat ) ) ,
						"M"= base::colMeans( dat , na.rm=TRUE )
								)
			# item and rater combinations
			rater_pars0 <- NULL
			I <- base::ncol(dat)
			R <- base::length(rater_unique)
			for (ii in 1:I){
				dfr <- base::data.frame( "item" = base::colnames(dat)[ii] , 
				    "rater" = rater_unique , 
					"rid" = 1:R )
				N_Rat <- stats::aggregate( 1 - base::is.na(dat[ , ii  ]) , 
										base::list(rater ) , base::sum )
				M <- stats::aggregate( dat[,ii] , base::list(rater) ,
											base::mean , na.rm=TRUE ) 
				ind <- match( dfr$rater , M[,1] )							
				dfr$N_Rat <- N_Rat[ ind , 2]
				dfr$M <- M[ ind , 2]				
				rater_pars0 <- base::rbind( rater_pars0 , dfr )
			}
		
			res <- base::list( pid = dat0$pid , rater = dat0$rater ,
							dat = dat0[ , - base::c(1,2) , drop=FALSE ] ,
							pid_unique = pid_unique , 
							rater_unique = rater_unique , item0 = item0 ,
							rater_pars0 = rater_pars0						
									)
			base::return(res)		
}
####################################################################				