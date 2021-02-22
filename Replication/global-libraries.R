####################################
# global libraries used everywhere #
####################################

mran.date <- "2020-11-30"
options(repos=paste0("https://cran.microsoft.com/snapshot/",mran.date,"/"))



pkgTest <- function(x)
{
	if (!require(x,character.only = TRUE))
	{
		install.packages(x,dep=TRUE)
		if(!require(x,character.only = TRUE)) stop("Package not found")
	}
	return("OK")
}

global.libraries <- c("dplyr","haven","foreign","ggplot2","ggrepel","gridExtra",'modeest")

results <- sapply(as.list(global.libraries), pkgTest)
