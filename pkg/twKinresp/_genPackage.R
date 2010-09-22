#------------ generating package tw.DEMC
# inlcuding generation of Rd files from inline-docs
# install.packages("inlinedocs")

library(RUnit)
library(R.methodsS3)
library(nlme)
#library(twNlme)
library(lmtest)
#library(snowfall)
library(inlinedocs)

#source(file.path("R","inlinedocs.R"))	# does not work on install time
tmp <- sapply(Sys.glob(file.path("R","*.R")), source)
#tmp <- sapply(Sys.glob(file.path("R","multiTemp.R")), sfSource)
#library(twMisc)
#data( list=twStripFileExt(basename(Sys.glob(file.path("data","*.RData")))))
data(respWutzler10)

#mtrace(twUtest)
#(res <- twUtestF(fitKinrespExperiment,divertOutputFile=NULL))

library(twMisc)
(res <- twUtestF())

tmp.f <- function(){
	library(inlinedocs)
	unlink( file.path("man","*.Rd") )	
	package.skeleton.dx(".")
}

tmp.f <- function(){
	unloadNamespace("twKinresp")
	library(twKinresp)
	?twKinresp
}

tmp.f <- function(){
	mtrace(package.skeleton.dx)
	mtrace(extract.docs.chunk)
	mtrace(modify.Rd.file)
	mtrace(modify.Rd.file,F)
}
#R CMD check --no-vignettes --no-latex --no-install twKinresp
#R CMD check --no-vignettes --no-latex --no-codoc twKinresp
#R CMD INSTALL --html twKinresp
#R CMD INSTALL --build twKinresp
