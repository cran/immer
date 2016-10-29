

###############################################################
# log-likelihood lc2_agreement
logLik.lc2_agreement <- function (object, ...) {
	# extract log-likelihood
	out <- object$loglike
    # number of parameters
    base::attr(out, "df") <- object$model_output$npars
	# extract number of observations
    base::attr(out, "nobs") <- object$nobs
    base::class(out) <- "logLik"
    base::return(out)
}
#################################################################