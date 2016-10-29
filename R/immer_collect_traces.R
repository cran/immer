
################################################################
# collect all traces
immer_collect_traces <- function( traces , est_settings ){

	eps <- 1E-5	
	# exclude non-estimated parameters from traces
	
	#*** a
	aM <- traces[[ "a" ]]
	sd1 <- base::max( base::apply( aM , 1 , stats::sd ) )
	if ( sd1 < eps ){
		traces$a <- NULL
					}
	if ( ! est_settings$est.a ){
		traces$a <- NULL
					}									
	#*** mu			
	if ( stats::sd(traces$mu) < eps ){
		traces$mu <- NULL
					}
	if ( ! est_settings$est.mu ){
		traces$mu <- NULL
					}					
	#*** sigma			
	if ( stats::sd(traces$sigma) < eps ){
		traces$sigma <- NULL
					}
	if ( ! est_settings$est.sigma ){
		traces$sigma <- NULL
					}
	#*** phi
	if ( est_settings$est.phi == "r"){
	    traces$phi <- traces$phi[1,,]
					}
	if ( est_settings$est.phi == "n"){
	    traces$phi <- NULL
					}					
	#*** psi
	if ( est_settings$est.psi == "r"){
	    traces$psi <- traces$psi[1,,]
					}
	if ( est_settings$est.psi == "n"){
	    traces$psi <- NULL
					}					
	if ( est_settings$est.psi == "i"){
	    traces$psi <- traces$psi[,1,]
					}
	if ( est_settings$est.psi == "e"){
	    traces$psi <- traces$psi[1,1,]
					}	
	
	
	#----------------------------				
	NV <- base::length( traces )
	BB <- attr( traces , "NSamples" )
	burnin <- attr(traces , "burnin")
	iter <- attr( traces , "iter" )
	
	dfr0 <- NULL
	for (vv in 1:NV){
		# vv <- 2
		trace_vv <- traces[[vv]]
		var_vv <- base::names(traces)[[vv]]    
		dfr1 <- base::data.frame( "index" = vv , "par" = var_vv )
		dim_vv <- base::dim( trace_vv )
		dfr1$dim <- base::paste0( dim_vv , collapse="-")
		dfr1$facets <- base::length(dim_vv ) 
		dfr1$Npars <- base::prod( dim_vv[ - ( dfr1$facets   ) ] )
		dfr0 <- base::rbind( dfr0 , dfr1 )
				}

	dfr0$start <- base::cumsum( base::c(0 , dfr0$Npars[-NV] ) )+1
	dfr0$end <- base::cumsum( dfr0$Npars ) 

	NP <- base::sum( dfr0$Npars)
	dfr <- base::matrix( NA , nrow=BB , ncol=NP)
	
	parnames <- NULL
	sds <- NULL

	for (vv in 1:NV){
		# vv <- 5
		trace_vv <- traces[[vv]]
		var_vv <- base::names(traces)[[vv]]    
		ldvv <- base::length( base::dim(trace_vv) )
		if ( ldvv == 3 ){	
			res_vv <- trace_array3( trace_vv , dfr0 , vv , var_vv , BB )		
								}
		if ( ldvv == 2 ){	
			res_vv <- trace_matrix( trace_vv , dfr0 , vv , var_vv , BB )		
								}							
		if ( ldvv == 0 ){	
			res_vv <- trace_scalar( trace_vv , dfr0 , vv , var_vv , BB )		
								}							
		ind <-  base::seq( dfr0$start[vv] , dfr0$end[vv] ) 	
		dfr[ ,ind ] <- res_vv$av1
		parnames <- base::c( parnames , res_vv$pars_vv )
		sds <- base::c( sds , res_vv$sd_vv )
			}
	base::colnames(dfr) <- parnames		
	dfr <- base::as.matrix(dfr)
	selec <- base::which( sds > eps )
	dfr <- dfr[ , selec ]	
	base::class(dfr) <- "mcmc"
	mcmcobj <- dfr
	# colnames(mcmcobj) <- paste0("b",1:3)
	# mcmcobj <- as.matrix(mcmcobj)
	# class(mcmcobj) <- "mcmc"
	base::attr(mcmcobj, "mcpar") <- base::c( burnin+1 , iter , ( iter  - burnin ) / BB )
	mcmcobj <- coda::as.mcmc.list( mcmcobj )		
	res <- base::list( mcmcobj = mcmcobj ,  info = dfr0 , BB = BB  )
	base::return(res)
}
	
	
#########################################
# helper functions

trace_array3 <- function( trace_vv , dfr0 , vv , var_vv , BB ){
	av <- base::aperm( trace_vv  , base::c(3,1,2) )
	av1 <- base::matrix( av , nrow=BB , ncol=dfr0$Npars[vv] )
	dim1 <- base::dim(trace_vv)[-3]
	pars_vv <- base::apply( base::expand.grid( 1:dim1[1]  , 1:dim1[2]  ) , 1 , 
						FUN = function(hh){
							base::paste0( var_vv  , "[" , hh[1] , "," , hh[2] , "]" )
						}   )
	base::colnames(av1) <- pars_vv	
	sd_vv <- base::apply( av1 , 2 , stats::sd )			
	res0 <- base::list( av1 = av1 , pars_vv = pars_vv , sd_vv=sd_vv )
	base::return(res0)
}
	
	
trace_matrix <- function( trace_vv , dfr0 , vv , var_vv , BB ){
	av <- base::t( trace_vv  )
	av1 <- base::matrix( av , nrow=BB , ncol=dfr0$Npars[vv] )
	pars_vv <- 	paste0( var_vv  , "[" , 1:dfr0$Npars[vv] ,  "]" )
	colnames(av1) <- pars_vv	
	sd_vv <- base::apply( av1 , 2 , stats::sd )	
	res0 <- base::list( av1 = av1 , pars_vv = pars_vv , sd_vv = sd_vv )
	base::return(res0)
}	


	
trace_scalar <- function( trace_vv , dfr0 , vv , var_vv , BB ){
	av1 <- trace_vv	
	pars_vv <- 	var_vv
	res0 <- base::list( av1 = av1 , pars_vv = pars_vv , sd_vv = stats::sd(av1) )
	base::return(res0)
}				