.setUp <- function(){
	data(respWutzler10)
}

.tearDown <- function(){
}

test.kinrespGrowthphaseReplicate <- function(){
	rder <- subset(respWutzler10, suite=="Face" & experiment==3 & replicate==2 )
	#mtrace(kinrespReplicate)
	res2 <- kinrespGrowthphaseReplicate(rder, weights=varPower(fixed=0.5)) 
}

