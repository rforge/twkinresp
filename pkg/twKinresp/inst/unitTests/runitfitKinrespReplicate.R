.setUp <- function(){
	data(respWutzler10)
	.setUpDf <- within( list(),{
			rder <- subset(respWutzler10, suite=="Face" & experiment==3 & replicate==2 )
			rder.e <- rder[1:which.max(rder$resp),]
		})
	attach(.setUpDf)
}

.tearDown <- function(){
	detach()
}

test.fitKinrespBetaReplicate <- function(){
	res2 <- fitKinrespBetaReplicate(rder.e$time, rder.e$resp, weights=varPower(fixed=0.5))
	summary(res2)
	#mtrace(.coefKinrespBeta)
	resBeta <- coefKinrespBeta(coef(res2))
	checkEquals( paste("beta",0:2,sep=""), names(resBeta))
	checkTrue( all(resBeta > 0) )
	checkEqualsNumeric( fitted(res2), modelKinrespBeta(rder.e$time,resBeta) )
	checkEquals( coef(res2), .coefKinrespBetaLogStart(resBeta) )
	
	#mtrace(calcKinrespCoef)
	res3 <- calcKinrespCoef( coef(res2) )
	#regression test from former result
	checkEqualsNumeric( c(mumax=0.20, r0=0.02, x0=118.48), res3, tol=0.01 )
	
	res4 <- coefKinresp(coef(res2))
	checkEquals( res3,  res4)
}

.tmp.f <- function(){
	plot(res2)
	plot( resp~time, data=rder.e)
	lines( fitted(res2)~rder.e$time)
	lines( modelKinrespBeta(rder.e$time,resBeta)~rder.e$time, col="red")
	
	lines( fitted(res10)~rder.e$time, col="blue")
}

test.kinrespReplicate <- function(){
	#mtrace(kinrespReplicate)
	res10 <- fitKinrespReplicate(rder.e, weights=varPower(fixed=0.5)) 
	summary(res10)
	resParms <- coefKinresp(coef(res10))
	checkTrue( all(c("mumax","r0","x0") %in% names(resParms)) )
	checkTrue( all(resParms > 0) )
	checkEquals( coef(res10), coefKinrespNormStart(resParms) )
	confintKinresp(confint(res10))
	#mtrace(kinrespParDist.gnls)
	(pars <- kinrespParDist(res10))
	iPar="x0"
	xGrid <- seq( pars[iPar,"cf025"]*0.8, pars[iPar,"cf975"]*1.2, length.out=80)
	fx <- dlnorm(xGrid, meanlog=pars[iPar,"mu"],sdlog=pars[iPar,"sigma"])
	plot( fx ~ xGrid, type="l" )
	abline(v=pars[iPar,c("mle","median","mean","cf025","cf975")], col=c("red","green","blue","gray","gray"))
	iPar="r0"
	xGrid <- seq( pars[iPar,"cf025"]*0.8, pars[iPar,"cf975"]*1.2, length.out=80)
	fx <- dlogitnorm(xGrid, mean=pars[iPar,"mu"],sd=pars[iPar,"sigma"])
	plot( fx ~ xGrid, type="l" )
	abline(v=pars[iPar,c("mle","median","mean","cf025","cf975")], col=c("red","green","blue","gray","gray"))
}

