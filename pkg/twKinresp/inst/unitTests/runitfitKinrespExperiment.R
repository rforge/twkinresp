.setUp <- function(){
	data(respWutzler10)
	.setUpDf <- within( list(),{
			rde <- subset(respWutzler10, suite=="Face" & experiment==3 )
		})
	attach(.setUpDf)
}

.tearDown <- function(){
	detach()
}

test.fitKinrespExperiment <- function(){
	plot( resp ~ time, data=rde, col=as.numeric(rde$rep) )
	#mtrace(kinrespGrowthphaseExperiment)
	res4 <- kinrespGrowthphaseExperiment(rde, weights=varPower(fixed=0.5) )
	# extract the dataset of unlimited growth
	rde.e <- getUnlimitedGrowthData(res4)
	# get the initial values (microbial parameters of replicate fits)
	coefRep <- coefKinresp(res4)
	# fit the mixed model
	#mtrace(fitKinrespExperiment)
	res5 <- fitKinrespExperiment( rde.e, coefRep )
	summary(res5$model)
	coefKinresp(fixef(res5$model))
	confintKinresp(confint(res5$model))
	#mtrace(momentsLogitnorm)
	(pars <- kinrespParDist(res5$model))
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


