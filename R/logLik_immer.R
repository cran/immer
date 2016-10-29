
###############################################################
# log-likelihood immer_HRM
logLik.immer_HRM <- function (object, ...) {
	# extract log-likelihood
	out <- object$like
    # number of parameters
    base::attr(out, "df") <- base::sum(object$ic$Npars)
	# extract number of observations
    base::attr(out, "nobs") <- object$ic$N
    base::class(out) <- "logLik"
    base::return(out)
}
#################################################################



###############################################################
# log-likelihood immer_cml
logLik.immer_cml <- function (object, ...) {
	# extract log-likelihood
	out <- object$loglike
    # number of parameters
    base::attr(out, "df") <- base::sum(object$npars)
	# extract number of observations
    base::attr(out, "nobs") <- object$N
    base::class(out) <- "logLik"
    base::return(out)
}
#################################################################