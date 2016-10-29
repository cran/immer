
################################################
# design matrix LPCM
lpcm_generate_design <- function( pars_info , irtmodel , W , 
			b_const , normalization , I , maxK , nullcats ){
	
		W0 <- W
		b0 <- b_const
		n1 <- base::nrow(pars_info)

		if ( base::is.null(irtmodel) ){
			irtmodel <- "PCM"
		}
		
		if ( base::is.null( b_const ) ){
				b_const <- base::rep( 0 , n1 )
		}		
		
		
		pars_info$estpar <- 1	
		if ( nullcats == "zeroprob" ){
			pars_info$estpar <- 1*(pars_info$Freq > 0	)
			b_const[ pars_info$estpar == 0 ] <- 99
								}
								
		#*********************
		if ( nullcats != "zeroprob" ){
		    n2 <- n1
			W <- base::matrix( 0 , nrow=n1 , ncol=n2-1)
			base::rownames(W) <- base::paste0( pars_info$item , "_Cat" , pars_info$cat )
		}

		if ( nullcats == "zeroprob" ){
		    n2 <- base::sum(pars_info$estpar)
			W <- base::matrix( 0 , nrow=n1 , ncol=n2-1)
			base::rownames(W) <- base::paste0( pars_info$item , "_Cat" , pars_info$cat )
			if ( base::is.null( b_const ) ){
				b_const[ pars_info$estpar == 0 ] <- 99
			}
			irtmodel <- "PCM"
		}
							
		#------------------------
		# irtmodel == "PCM"
		if (irtmodel == "PCM"){
			pinfo2 <- pars_info[ pars_info$estpar == 1 , ]
			n1 <- base::nrow(pinfo2)
			n2 <- n1	
			index <- pinfo2$index
			# PCM: normalization = "first"
			if ( normalization == "first" ){
				W[ base::cbind( index[2:n1] , 1:(n2-1) ) ] <- 1
				base::colnames(W) <- base::rownames(W)[index[-1]	]
			}
											
			# PCM: normalization = "sum"
			if ( normalization == "sum" ){
				W[ base::cbind( index[1:(n1-1)] , 1:(n1-1) ) ] <- 1
				W[ index[n1] , ] <- -1 
				base::colnames(W) <- base::rownames(W)[index[-n1]]	
			}
		}
		
		#--------------------------
		# irtmodel == "PCM2"
		if (irtmodel == "PCM2" ){
		    items <- base::unique( base::paste(pars_info$item))
			I <- base::max( pars_info$itemid )
			base::colnames(W) <- paste0("w",1:(n1-1))
			#--- normalization == "first"
			if ( normalization == "first"){
				# items
				p1 <- pars_info[ pars_info$itemid > 1	, ]
				W[ base::cbind( p1$index , p1$itemid - 1 ) ] <- p1$cat
				base::colnames(W)[ base::seq(1 , I-1 ) ] <- items[-1]
			}
			#--- normalization == "sum"
			if ( normalization == "sum"){
				# items
				p1 <- pars_info[ pars_info$itemid < I	, ]
				W[ base::cbind( p1$index , p1$itemid  ) ] <- p1$cat
 		        base::colnames(W)[ base::seq(1 , I-1 ) ] <- items[-I]				
				p1b <- pars_info[ pars_info$itemid == I , , drop=FALSE ]
				for ( kk in base::seq(1,nrow(p1b) ) ){
					W[ p1b$index[kk] , 1:(I-1) ] <- - p1b$cat[kk]
				}
			}				
			
			vv <- I 
			p2 <- pars_info
			p2$param <- 0
			p2$param[ p2$cat < p2$maxK ] <- 1
			p2$param <- ( p2$param > 0 ) * ( base::cumsum( p2$param ) + ( vv - 1 ) )
			W[ base::cbind( p2$index , p2$param ) ] <- 1
			p2a <- p2[ p2$param > 0 , ]
			base::colnames(W)[ p2a$param ] <- base::paste0( p2a$item , "_Step" , p2a$cat )					
		}

		if ( ! base::is.null(W0) ){
			W <- W0 
		}
		if ( ! base::is.null(b0) ){
			b_const <- b0
		}			
														
		if ( base::is.null( base::colnames(W) ) ){
			base::colnames(W) <- base::paste0("par" , base::seq(1,base::ncol(W)) )
		}
				
		#*********************
		# output
		res <- base::list(W=W , b_const=b_const, irtmodel=irtmodel)
		base::return(res)		
}
###############################################