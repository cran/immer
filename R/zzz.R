#  zzz.R
#
# This function is simply copied from mice package.

# on attach immer
.onAttach <- function(libname,pkgname){
  d <- utils::packageDescription("immer")
  d1 <- d$Version
  base::packageStartupMessage(
		base::paste("* " ,d$Package," " , d1 ," (",d$Date,")",sep="")  )
	}
	
	
version <- function(pkg="immer"){
  lib <- base::dirname(system.file(package = pkg))
  d <- utils::packageDescription(pkg)
  base::return( base::paste(d$Package,d$Version,d$Date,lib) )
}

# .First.lib <- function(lib, pkg){
#          library.dynam("sirt", package = pkg, lib.loc = lib)
#          return(invisible(0))
#        } 