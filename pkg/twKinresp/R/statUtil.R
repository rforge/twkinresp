# from package twNlme

.twCalcr2 <- function(
	### Calculate coefficient of determination (r-square).
	ypred	##<< variable 1, i.e. predicted values
	,y		##<< variable 2, i.e. observed values, must have the same length as ypred
	,weights=1	##<< weights, also same length as ypred
){
	bo <- !is.na(ypred) & !is.na(y) 
	1 - sum(((ypred-y)[bo])^2*weights) / sum( (y[bo] - mean(y[bo]))^2*weights )	
}

setMethodS3("fixef","gnls", function(
	### Fixed coefficients of either nlme or gnls (S3 method).
	object	##<< model, whose fixed coefficients should be accessed
	,...
){
	# fixef.gnls
	##details<< 
	## To support similar access to fixed coefficients of nlme and gnls models.
	## Because coef(nlme) will give more details than fixed.effects 
	coef(object)
})

.twVarFuncCoef <- function(
	### Extracting the estimated coefficients of the variance functions, such as power in varPower
	tmp.m	##<< results of gnls or nlme
){
	tmp.mslength = length(tmp.m$modelStruct)
	coef(tmp.m$modelStruct[[tmp.mslength]], uncons = FALSE, allCoef = TRUE )	
}

#support similar access to the correlation matrix of gnls and nlme
#extractCorrMatrix <- function(x, ...){ UseMethod("extractCorrMatrix") }

setMethodS3("extractCorrMatrix","gnls", function( 
	### Access to correlation matrix of gnls.
	modelfit	##<< result of gnls
	,...
){
	summary(modelfit)$corBeta
})

setMethodS3("extractCorrMatrix","nlme", function( 
	### Access to correlation matrix of nlme.
	modelfit	##<< result of nlme
	,...
){
	summary(modelfit)$corFixed
})

.extractRandomCovMatrix <- function(
	### Extract the Covariance Matrix for the Random effects.
	tmp.fit		##<< result of gnls or nlme
){
	pdMatrix(tmp.fit$modelStruct$reStruct[[1]])*tmp.fit$sigma^2	
}

.extractNlmeFullVar <- function(
	### Extract var = VarFixed + VarRand from nlme or gnls.
	tmp.fit	##<< result of gnls or nlme
){
	diag(tmp.fit$varFix) + diag( .extractRandomCovMatrix(tmp.fit) )
}

setMethodS3("confint","nlme", function( 
	### confidence intervals for the parameters of nlme fit, population level. 
	object
	,...
){
	#confint.nlme
	#tmp.se <- sqrt.extractNlmeFullVarr(tmp.fit))
	tmp.se <- sqrt(diag(object$varFix)) #only interested in the population interval
	tmp.df <- object$dims$N + 
			-object$dims$ncol[object$dims$Q + 1]+ # number of parameters
			-length(coef(object[["modelStruct"]]))+
			-1
	#tmp.t <- qt(0.975, tmp.fit$fixDF$X )
	tmp.t <- qt(0.975, tmp.df )
	tmp <- cbind( fixed.effects(object)  -tmp.t*tmp.se, fixed.effects(object)+tmp.t*tmp.se )
	dimnames(tmp)[[2]] <- c("2.5%","97.5%")
	tmp
})






