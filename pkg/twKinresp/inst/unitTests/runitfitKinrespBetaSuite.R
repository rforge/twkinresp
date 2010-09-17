.setUp <- function(){
	data(respWutzler10)
}

.tearDown <- function(){
}

scratch.kinrespGrowthphase <- function(){
	## get the 3 experiments from two different series (to constrain processing time)
	#rd <- subset(respWutzler10, experiment %in% levels(respWutzler10$experiment)[15:17] )
	#rd <- subset(respWutzler10, experiment %in% levels(respWutzler10$experiment)[9:11] )
	rd <- subset(respWutzler10, suite=="Face" )
	#mtrace(kinrespGrowthphase)
	#weights=varPower(fixed=0.5)
	#res1 <- kinrespGrowthphase(rd, weights=varPower(fixed=0.5))
	res2 <- kinrespGrowthphase(rd, weights=varPower(fixed=0.5))
	res3 <- kinrespGrowthphase(rd, weights=varPower(fixed=0.3))
	res4 <- kinrespGrowthphase(rd, weights=varPower(fixed=0.8))
	
	checkTrue( all(c("n","resRep","errors") %in% names(res1) ))
	checkEquals( length(res1$resRep), nrow(res1$n) )  # regression test (from former run)

	
	rds.e <- getUnlimitedGrowthData(res1)
	repFits <- coefKinrespBeta(res1)
	rds.e <- getUnlimitedGrowthData(res2)
	repFits <- coefKinrespBeta(res2)
	rds.e <- getUnlimitedGrowthData(res3)
	repFits <- coefKinrespBeta(res3)
	rds.e <- getUnlimitedGrowthData(res4)
	#mtrace(coefKinrespBeta.kinrespList)
	repFits <- coefKinrespBeta(res4, rds.e)	#0.47
	#coefList(res2)

	
	rd <- subset(respWutzler10, suite=="Fal" )
	res2b <- kinrespGrowthphase(rd, weights=varPower(fixed=0.5))
	rds.e <- getUnlimitedGrowthData(res2b)
	repFits <- coefKinrespBeta(res2b, rds.e)
	#mtrace(fitKinrespBetaSuite)
	tmp <- fitKinrespBetaSuite( rds.e, repFits, weights=varPower(0.5))	#0.34

	rd <- subset(respWutzler10, suite=="Pushchino" )
	res3b <- kinrespGrowthphase(rd, weights=varPower(fixed=0.5))
	rds.e <- getUnlimitedGrowthData(res3b)
	#mtrace(coefKinrespBeta.kinrespList)
	repFits <- coefKinrespBeta(res3b, rds.e)
	#mtrace(fitKinrespBetaSuite)
	tmp <- fitKinrespBetaSuite( rds.e, repFits, weights=varPower(0.5))
	# does not converge

	rd <- subset(respWutzler10, suite=="Pushchino" )
	res4b <- kinrespGrowthphase(rd, weights=varPower(fixed=0.3))
	rds.e <- getUnlimitedGrowthData(res4b)
	#mtrace(coefKinrespBeta.kinrespList)
	repFits <- coefKinrespBeta(res4b, rds.e)
	#mtrace(fitKinrespBetaSuite)
	tmp <- fitKinrespBetaSuite( rds.e, repFits, weights=varPower(0.3))
	#does not converge
	
	# tracing error of fit==NULL
	# options(error=recover) # find the prolemnatic experiment
	rd <- subset(respWutzler10, suite=="Pushchino" & experiment==42 )
	res3b <- kinrespGrowthphase(rd, weights=varPower(fixed=0.5))
	plotKinrespDiagnostics( res3b$resRep[["Pushchino_42_3"]] )
	#mtrace(plotKinrespDiagnostics.kinresp)
	#mtrace(.plotRepResidualArray)
	plotKinrespDiagnostics( res3b$resRep[["Pushchino_42_3"]], rangeN=c(33,40), nCol=4 )
	plotKinrespDiagnostics( res3b$resRep[["Pushchino_42_2"]] )
	rds.e <- getUnlimitedGrowthData(res3b)
	#mtrace(coefKinrespBeta.kinrespList)
	repFits <- coefKinrespBeta(res3b, rds.e)
	
}

