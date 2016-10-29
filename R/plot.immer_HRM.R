
plot.immer_HRM <- function( x , ... ){
	base::class(x) <- "mcmc.sirt"
	graphics::plot( x , ... )
}
	
