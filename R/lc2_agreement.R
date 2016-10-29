
######################################################
# latent class agreement
# Schuster & Smith (2006)
lc2_agreement <- function( y , w = rep(1,nrow(y)) , type="homo" , 
   method = "BFGS" , ...  ){

	CALL <- base::match.call()
	s0 <- base::Sys.time()	
	#*** preprocessing for dataset		
	if ( base::is.null( base::colnames(y) ) ){
		base::colnames(y) <- base::c("Var1","Var2")
	}	
	
	res <- immer_unique_patterns( dat = y , w = w )
	y <- res$y
	w <- res$w
		
	y0 <- y	
	values <- base::sort( base::unique( base::c(y[,1] , y[,2] ) ) )
	I <- base::ncol(y)
	for (ii in 1:I){
		y[,ii] <- base::match( y[,ii] , values )
	}
	
	#*** preprocessing
	# The input is a weighted frequency table.
	W <- base::sum(w)
	eps <- 1E-20
	# define starting values
	p0 <- base::sum( w[ y[,1] == y[,2] ] )  / W
	gamma0 <- p0 * 2/3
	# compute tables
	vals <- base::unique( c(y[,1] ,y[,2] )  )
	V <- base::length(vals)
	# compute frequencies
	relFreq <- base::rep(0,V)
	base::names(relFreq) <- vals
	for (vv in 1:V){
		# vv <- vals[1]
		relFreq[vv] <- ( base::sum( w[  y[,1] == vv ] ) / W + 
		                       base::sum( w[  y[,2] == vv ] ) / W ) / 2
	}
	#--------------------------------------				
	# define parameter table
	#**** tau parameters
	tau_i <- relFreq
	tau_i_logits <- probs2logits( tau_i )
	dfr1 <- base::data.frame( "pargroup" = "tau" , "parnum" = NA , 
		        "parname" = paste0( "tau_Cat" , values ) ,
			    "est" = TRUE , "est_parindex" = 0 , 
			    "rater" = 0  , "value" = NA , "value_logits" = NA , "start_logits" = NA)
	dfr1[1,"est"] <- FALSE
	dfr1$est_parindex <- base::cumsum(dfr1$est)
	dfr1$start_logits <- tau_i_logits
	dfr0 <- dfr1

	#**** phi parameters	
	dfr1 <- base::data.frame( "pargroup" = "phi" , "parnum" = NA , 
				"parname" = paste0( "phi_Cat" , values ) ,
				"est" = TRUE , "est_parindex" = 0 , 
				"rater" = 0  , "value" = NA , "value_logits" = NA , "start_logits" = NA)
	#--- homogeneous
    if (type=="homo"){
		dfr1[1,"est"] <- FALSE	
		dfr1$est_parindex <- base::cumsum(dfr1$est) + base::max(dfr0$est_parindex)
		dfr1$start_logits <- tau_i_logits
	}		
	#--- equal phi parameters
	if (type=="unif"){
		phi0 <- base::rep(1/V,V)
		dfr1$start_logits <- probs2logits(phi0)
		dfr1$est <- FALSE
	}
	#--- set tau and phi parameters equal to each other
	if ( type == "equal"){
	    dfr1$est[1] <- FALSE
		dfr1$est_parindex <- dfr0$est_parindex		
		dfr1$start_logits <- dfr0$start_logits
	}
	#--- heterogeneous raters
	if (type == "hete"){
		dfr1[1,"est"] <- FALSE		
		dfr1a <- dfr1
		dfr1a$parname <- base::gsub( "phi_" , paste0("phi_" , base::colnames(y)[1] ),
				"_" , dfr1a$parname)		
		dfr1a$est_parindex <- base::cumsum(dfr1a$est) + base::max(dfr0$est_parindex)
		dfr1a$start_logits <- tau_i_logits	
		dfr1b <- dfr1
		dfr1b$parname <- base::gsub( "phi_" , paste0("phi_" , base::colnames(y)[2]) ,
				"_" , dfr1a$parname)		
		dfr1b$est_parindex <- base::cumsum(dfr1b$est) + base::max(dfr1a$est_parindex)
		dfr1b$start_logits <- tau_i_logits			
		dfr1 <- base::rbind( dfr1a , dfr1b)
	}
				

	dfr0 <- base::rbind( dfr0 , dfr1)
	
	#**** gamma parameters
	dfr1 <- base::data.frame( "pargroup" = "gamma" , "parnum" = NA , 
	            "parname" = base::paste0( "gamma_" , c(0,1) ) ,
				"est" = c(TRUE,FALSE) , "est_parindex" = 0 , 
				"rater" = 0  , "value" = NA , "value_logits" = NA , "start_logits" = NA)
	dfr1$est_parindex <- base::cumsum(dfr1$est) + base::max(dfr0$est_parindex)
	dfr1$start_logits <- base::c( stats::qlogis( gamma0 ) , 0 )
	dfr0 <- base::rbind( dfr0 , dfr1)
	dfr0[ ! dfr0$est , "est_parindex" ] <- 0	
	dfr0$parnum <- base::seq(1 , base::nrow(dfr0) )
		
	# extract starting values
	x0 <- dfr0[ dfr0$est ,]
	x0 <- x0[ ! base::duplicated( x0$est_parindex) , ]
	x <- x0[ x0$est , "start_logits" ] 

		#***********************************************************
		# optimization function
		l2rater_fct <- function(x){
			# reconstruct phi parameters
			ind <- dfr0[ ( dfr0$pargroup == "tau" )  & 
			             ( dfr0$est_parindex > 0 ) , "est_parindex" ]
			tau00 <- base::c( 0 , x[ind] )
			tau_probs <- logits2probs( tau00 )			
			# reconstruct phi parameters
			if ( type=="unif"){
				ind <- dfr0[ dfr0$pargroup == "phi" ,] 
				phi00 <- dfr0$start_logits
			}
			if ( type %in% base::c("homo") ){
  			    ind <- dfr0[ ( dfr0$pargroup == "phi" )  & 
			               ( dfr0$est_parindex > 0 ) , "est_parindex" ]			
				phi00 <- base::c(0,x[ind])
			}
			if ( type %in% base::c("equal") ){
  			    ind <- dfr0[ ( dfr0$pargroup == "phi" ) & ( dfr0$est ) , "est_parindex" ]
				phi00 <- base::c(0,x[ind])		
			}	
			if (type == "hete"){
  			    ind <- dfr0[ ( dfr0$pargroup == "phi" ) , "est_parindex" ]	
				L <- base::length(ind)
				L1 <- L / 2
				ind[ ind == 0] <- NA
				phi00 <- x[ ind ]
				phi00[ is.na(phi00) ] <- 0				
				phi_probs1 <- logits2probs( phi00[ base::seq(1,L1) ] )
				phi_probs2 <- logits2probs( phi00[ L1 + base::seq(1,L1) ] )
			}				
			if (type != "hete"){
			    phi_probs <- logits2probs( phi00 )
				phi_probs1 <- phi_probs
				phi_probs2 <- phi_probs
			}			
			# reconstruct gamma parameters
			ind <- dfr0[ ( dfr0$pargroup == "gamma" )  & 
			               ( dfr0$est_parindex > 0 ) , "est_parindex" ]
			gamma_probs <- logits2probs( base::c(x[ind] , 0 ) )
			# estimated probability
			est_prob <- ( y[,1] == y[,2]  ) * tau_probs[ y[,1] ] * gamma_probs[1] +
							phi_probs1[ y[,1] ] * phi_probs2[ y[,2] ] * ( 1 - gamma_probs[1] )
			ll <- - base::sum( w * base::log( est_prob + eps ) )
			base::return(ll)
		}
		#************************************************************	
	
	#***************
	# conduct optimization fitted model
	h1 <- stats::optim( x , l2rater_fct , method=method , hessian = TRUE , ... )
		
	#------ optimization independence model
		#***********************************************************
		# optimization function
		l2_indep <- function(x){
			phi_probs <- logits2probs( base::c(0, x ) )
			phi_probs1 <- phi_probs2 <- phi_probs
			est_prob <- phi_probs1[ y[,1] ] * phi_probs2[ y[,2] ]
			ll <- - base::sum( w * base::log( est_prob + eps ) )
			base::return(ll)
		}
		#************************************************************	
	h2 <- stats::optim( x[1:(V-1)] , l2_indep , method=method , hessian = FALSE , ... )
	
	#***************
	# update parameter table
	m1 <- dfr0$est_parindex
	m1[ m1 == 0] <- NA
	dfr0[  , "value_logits" ] <- h1$par[ m1 ]
	dfr0[ base::is.na(dfr0$value_logits) , "value_logits"] <- 0
	# tau parameters
	ind <- base::which(dfr0$pargroup == "tau")
	y <- dfr0[ ind , "value_logits" ]
	dfr0[ ind , "value"] <- logits2probs(y)
	
	ind <- base::which(dfr0$pargroup == "phi")
	if ( type != "hete"){
		y <- dfr0[ ind , "value_logits" ]
		dfr0[ ind , "value"] <- logits2probs(y)
	} else {
		y <- dfr0[ ind , "value_logits" ]
		y1 <- y[1:V]
		y2 <- y[V + 1:V]
		dfr0[ ind[1:V] , "value"] <- logits2probs(y1)
		dfr0[ ind[V + 1:V] , "value"] <- logits2probs(y2)														
	}
		
	ind <- base::which(dfr0$pargroup == "gamma")
	y <- dfr0[ ind , "value_logits" ]
	dfr0[ ind , "value"] <- logits2probs(y)
	dfr0[ dfr0$pargroup == "gamma" , "value"]
	# compute agreement by chance
	phi_probs <- dfr0[ dfr0$pargroup == "phi" , "value"]
	phi_probs1 <- phi_probs2 <- phi_probs
	if ( type == "hete"){
	   phi_probs1 <- phi_probs[1:V]
	   phi_probs2 <- phi_probs[V + 1:V]	   
	}
	
	# jj <- phi_probs1[ y0[,1] ] * phi_probs2[ y0[,2] ]
	agree_chance <- base::sum( phi_probs1 * phi_probs2 ) * 
				( 1 - dfr0[ dfr0$parname == "gamma_0" , "value" ] )

	#**** saturated model
	ll0 <- - base::sum( w * base::log( w / W ) )
	#******
	# collect model results
	model_output <- base::data.frame( "deviance" = 2*h1$value)
	model_output$npars <- base::max( dfr0$est_parindex ) 
	saturated_output <- base::data.frame("deviance" = 2*ll0 , "npars" = nrow(y0) - 1)
	
	#*****
	# independence output
	independence_output <- base::data.frame("deviance" = 2*h2$value , 
			"npars" = base::length(h2$par)  )
	
	#**** likelihood ratio test
	LRT_output <- base::data.frame("chisquare" = -2*(ll0 - h1$value) )
	LRT_output$df <- saturated_output$npars - model_output$npars
	LRT_output$p <- 1 - stats::pchisq( LRT_output$chisquare , 
							     df = LRT_output$df )
	#*** normed fit index (Clogg, xxxx)
	# see Uebersax (1990, Statistics in Medicine)
	L1 <- LRT_output$chisquare
	L0 <- -2*(ll0 - h2$value)
	NFI <- ( L0 - L1 ) / L0	
	
	#**** agreement statistics
	gamma0 <- dfr0[ dfr0$parname == "gamma_0" , "value" ] 
	# conditional probability of true agreement given observed agreement
	rel_agree <- gamma0 / ( gamma0 + agree_chance )
	#**** information criteria
	# total number of parameters
	dev <- model_output$deviance 
	ic <- base::data.frame( "dev" = dev)	
	ic$n <- W
	ic$Npars <- ic$np <- model_output$npars
    # AIC
    ic$AIC <- dev + 2*ic$np
	# AIC3
	ic$AIC3 <- dev + 3*ic$np
    # BIC
    ic$BIC <- dev + ( base::log(ic$n) )*ic$np
	# adjusted BIC 
	ic$aBIC <- dev + ( base::log( ( ic$n -2 ) / 24 ) )*ic$np
    # CAIC (consistent AIC)
    ic$CAIC <- dev + ( base::log(ic$n) + 1 )*ic$np
	# corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	#---------------
	# parameter output
	
	if ( type == "hete"){ NR <- 3} else { NR <- 2}
	dfr11 <- base::as.data.frame( base::matrix( 0 , nrow=2 , ncol=V+1 ) )
	base::colnames(dfr11) <- base::c( "parm" , base::paste0("Cat",values) )
	
	dfr11[1,-1] <- dfr0[ dfr0$pargroup == "tau" , "value"]
	dfr11[1,1] <- "tau"
	if ( type != "hete"){
		dfr11[2,-1] <- dfr0[ dfr0$pargroup == "phi" , "value"]
		dfr11[2,1] <- "phi"		
	} else {
		dfr11[2,-1] <- dfr0[ dfr0$pargroup == "phi" , "value"][ base::seq(1,V)]
		dfr11[2,1] <- base::paste0( "phi_" , base::colnames(y0)[1] )								 
		dfr11[3,-1] <- dfr0[ dfr0$pargroup == "phi" , "value"][ V + base::seq(1,V)]
		dfr11[3,1] <- base::paste0( "phi_" , base::colnames(y0)[2] )
	}
	
	#*** time
	s1 <- base::Sys.time()
	time0 <- base::data.frame("start"=s0 , "end"=s1)
	
	# output
	res <- base::list( 
	    "model_output" = model_output ,
		"saturated_output" = saturated_output , 
		"independence_output" = independence_output , 
		"LRT_output" = LRT_output , 
		"NFI" = NFI , 
	    "partable" = dfr0 ,
		"parmsummary" = dfr11 , 
		"agree_true" = gamma0 ,
		"agree_chance" = agree_chance ,
		"rel_agree" = rel_agree , 
		"optim_output" = h1 ,
		"nobs" = W , "type" = type ,
		"ic" = ic , 
		"loglike" = - model_output$deviance / 2 ,
		"npars" = model_output$df ,
		"y" = y0 , "w" = w , 
		"CALL" = CALL , "time" = time0 ,
		"description" = "Latent class model for 2 raters (Schuster & Smith, 2006)"
			)
	base::class(res) <- "lc2_agreement"
	base::return(res)
}
#########################################################		

