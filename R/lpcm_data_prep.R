
###################################################
# data preparation function
# linear logistic partial credit model
lpcm_data_prep <- function( dat , weights , a ){
    used_persons <- base::as.vector( base::which( base::rowSums( 1 - base::is.na(dat) ) > 1 ))
	dat <- dat[ used_persons , ]
	if ( ! base::is.null(weights) ){
		weights <- weights[ used_persons ]
	}
	N <- base::nrow(dat)
	I <- base::ncol(dat)
	if ( is.null(a) ){
		a <- base::rep(1,I)
	}
	aM <- base::matrix( a , nrow=N , ncol=I , byrow=TRUE )
	dat <- dat * aM
	
	dat_ind <- 1 * base::is.na(dat)
	res <- sirt::md.pattern.sirt(dat_ind)
	resp_patt <- res$resp_patt
	
	
    if ( base::is.null(weights) ){
		weights <- base::rep(1 , N )
							}
	patt_unique <- base::unique( resp_patt )
	patt <- base::match( resp_patt , patt_unique )
	NP <- base::length(patt_unique)
	maxK <- base::apply( dat , 2 , base::max , na.rm=TRUE ) 
	item_index <- base::sapply( 1:NP , FUN = function(pp){
			# pp <- 1
			dat0 <- dat_ind[  patt == pp , , drop=FALSE ][1,]
			base::which( dat0 == 0 )
					} , simplify = FALSE )
	maxscore <- base::sapply( 1:NP , FUN = function(pp){
		base::sum(maxK[ item_index[[pp]] ] )
			} , simplify=FALSE )
	pars <- base::rep( 1:I , maxK )		
	pars_info <- base::data.frame( "item" = base::rep( base::colnames(dat),maxK) , 
							"itemid" = pars )
	pars_info$cat <- base::unlist( base::sapply( 1:I, FUN = function(ii){
						base::seq( 1 , maxK[ii] ) 
							}  , simplify=FALSE) )
	pars_info$index <- base::seq(1 , base::nrow(pars_info) )
	pars_info$maxK <- maxK[ pars_info$itemid ]
	pars_info$Freq <- 0
	
	K <- base::max( maxK )
	for (kk in 1:K){
    	f1 <- base::colSums( ( dat == kk ) * weights , na.rm=TRUE)
		ind <- pars_info[ pars_info$cat == kk , "itemid" ]
		pars_info[ pars_info$cat == kk , "Freq" ] <- f1[ind]
					}
	parm_index <- base::sapply( 1:NP , FUN = function(pp){
		 base::which( pars_info$itemid %in% item_index[[pp]] ) 
			} , simplify=FALSE )
			
			
	# calculate raw score
	rs <- base::rowSums( dat , na.rm=TRUE)
	score_freq <- base::sapply( 1:NP , FUN = function( pp ){
			# pp <- 2
			base::sapply( 0:maxscore[[pp]] , FUN = function(ss){
								base::sum( ( rs == ss ) * ( patt == pp ) * weights , na.rm=TRUE) } ,
								simplify = TRUE )
						} , simplify=FALSE )

	suffstat <- base::sapply( 1:NP , FUN = function(pp){
			base::unlist( base::sapply( item_index[[pp]] , FUN = function(ii){
						base::sapply( 1:maxK[ii] , FUN = function(kk){
								base::sum( ( dat[,ii] == kk ) * weights * ( patt == pp ) , na.rm=TRUE) 
												} , simplify=FALSE )
											}
							) ) } , simplify = FALSE )

	splitvec <- base::sapply( 1:NP , FUN = function(pp){
					mv <- item_index[[pp]] 
					oj_max <- maxK[ mv ]
					base::rep.int( mv , oj_max)  
						} , simplify=FALSE )

	# generate pararameter names
	parnames <- base::paste0( pars_info$item , "_Cat" , pars_info$cat )	
	base::rownames(pars_info) <- parnames					
						
	res <- base::list( N=N , I=I , NP=NP , dat=dat , 
			      patt=patt , weights=weights ,
				  suffstat=suffstat, splitvec=splitvec,
				  item_index = item_index , parm_index=parm_index ,
				  pars_info=pars_info , maxscore=maxscore ,
				  maxK=maxK , score=rs , used_persons = used_persons ,
				  score_freq = score_freq , a = a , parnames = parnames )	
	base::return(res)	    
}
###################################################