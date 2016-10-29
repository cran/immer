
########################################################
# lpcm initial parameters
# This function is mainly copied from the
# pcmodel function from the psychotools package
lpcm_inits <- function( dat , weights , maxK , b_const , W ,
		irtmodel , normalization ){
	I <- base::ncol(dat)
	m <- I
	oj_vec <- base::sapply( 1:I , FUN = function(ii){
				base::seq( 0 , maxK[ii] ) } , simplify=FALSE )
	y <- dat
	ctot <- base::vector("list", length = m)
    for (j in base::seq_len(m)){
		ctot[[j]] <- base::as.vector( base::tapply( weights, factor(y[, j], levels = oj_vec[[j]]), sum))
	}
    start <- base::lapply(ctot, function(x){
				- base::cumsum(diff.default( base::log(x))) } ) 
			# delta_jk = log(x_j(k-1) - log(x_jk), beta_jk = cumsum_k=1^k(delta_jk)
    start <- base::unlist(start)
    start[ is.na(start) ] <- 0
	par_init <- ( start - b_const ) %*% W
    par_init <- par_init[1,]
	if ( ( irtmodel=="PCM2") & ( normalization == "sum") ){ 
			par_init <- 0 * par_init
	}
    base::return(par_init)
}
########################################################