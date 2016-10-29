
#######################################################
# information criteria
immer_ic_hrm <- function( ic , summary.mcmcobj ){
		ic$n <- ic$N
		pars <- base::paste(summary.mcmcobj$parameter)		
		vars <- base::c("mu" , "sigma" , "a" , "b" ,  "phi" , "psi")
		VV <- base::length(vars)
		Npars <- NULL
		for (vv in 1:VV){
			# vv <- 1		
			ind <- base::which( base::substring( pars , 1 , base::nchar( vars[vv] ) ) == vars[vv] )
			Npars[ vars[vv] ] <- base::length(ind)
						}
		ic$Npars <- Npars
		ic$np <- base::sum(Npars)
		ic <- immer_IC_calc(ic)
		base::return(ic)
}
############################################################				