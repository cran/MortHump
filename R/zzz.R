
.onAttach <- function(libname, pkgname){
	suppressWarnings(descr <- utils::packageDescription("MortHump"))
# 	if(utils::packageVersion("MortHump")$minor %% 2 == 0) {
# 		state <- "stable"
# 	}
# 	else {
#  		state <- "development"
#  	}
  state <- ""
	if(!is.null(descr$Built)){
		builtDate <- paste(" (Built: ", strsplit(strsplit(descr$Built, ";")[[1]][3], " ")[[1]][2], ")", sep="")
	}else{
		builtDate <- ""
	}
	packageStartupMessage("This is MortHump ", state, " version ", descr$Version, builtDate)
	packageStartupMessage("\nTo cite MortHump in publications please use:")
	packageStartupMessage("Remund, Camarda and Riffe. Analyzing the Young Adult Mortality Hump in R with MortHump.")
	packageStartupMessage("MPIDR Technical Report TR-2018-3. (2018).")
}
