


##############################################################
# anova method immer_cml
anova_immer <- function( object , ... ){
    cl2 <- base::paste(match.call())[-1]
    if (length(list(object, ...)) != 2){ 
        base::stop("anova method can only be applied for comparison of two models.\n")		
	}
	objects <- base::list(object, ...)
	model1 <- objects[[1]]
	model2 <- objects[[2]]

	# define some necessary parameters
    model1$loglike <- model1$loglike
	model1$Npars <- model1$npars
    model2$loglike <- model2$loglike
	model2$Npars <- model2$npars
	# test
    dfr1 <- base::data.frame( "Model" = cl2[1] , 
		"loglike" = model1$loglike , 
		"Deviance" = -2*model1$loglike )
    dfr1$Npars <- base::sum(model1$Npars)
    dfr2 <- base::data.frame( "Model" = cl2[2] , 
		         "loglike" = model2$loglike , "Deviance" = -2*model2$loglike )
    dfr2$Npars <- base::sum(model2$Npars)
    dfr <- base::rbind( dfr1 , dfr2 )
    dfr <- dfr[ base::order( dfr$Npars ), ]
    dfr$Chisq <- NA
    dfr$df <- NA
    dfr$p <- NA
    dfr[1,"Chisq"] <- dfr[1,"Deviance"] - dfr[2,"Deviance"]
    dfr[1,"df"] <- base::abs( dfr[1,"Npars"] - dfr[2,"Npars"] )
    dfr[ 1, "p" ] <- base::round( 1 - stats::pchisq( dfr[1,"Chisq"] , df= dfr[1,"df"] ) , 5 )
    for ( vv in 2:( base::ncol(dfr))){ dfr[,vv] <- base::round( dfr[,vv] , 5 ) }
    base::print( dfr )
    base::invisible(dfr)
}
##############################################################

anova.immer_cml <- anova_immer
# anova.immer_HRM <- IRT.anova.immer_HRM
anova.immer_HRM <- anova_immer
anova.lc2_agreement <- anova_immer