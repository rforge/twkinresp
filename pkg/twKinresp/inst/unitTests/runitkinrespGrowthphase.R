.setUp <- function(){
	data(respWutzler10)
}

.tearDown <- function(){
}

test.kinrespGrowthphase <- function(){
	# get the 3 experiments from two different series (to constrain processing time)
	rd <- subset(respWutzler10, experiment %in% levels(respWutzler10$experiment)[15:17] )
	#mtrace(kinrespGrowthphase)
	res1 <- kinrespGrowthphase(rd, weights=varPower(fixed=0.5))
	checkTrue( all(c("n","resRep","errors") %in% names(res1) ))
	checkEquals( length(res1$resRep), nrow(res1$n) )  # regression test (from former run)
	#checkEquals( 19, length(res1$resRep) )  # regression test (from former run)
	checkEquals( 9, length(res1$resRep) )  # regression test (from former run)
	
	suppressWarnings(dir.create("tmp"))
	windows(pointsize=10,width=6.4,height=6.4, record=TRUE)
	plotKinrespDiagnostics(res1, plotFileBasename = file.path("tmp/testKinrespGrowthphase") )
	
	rd.e <- getUnlimitedGrowthData(res1)
	tmp.n <- list()
	tmp.n[ as.character(res1$n[1,"ser"]) ] <- res1$n[1,"n"]-1
	rd.e2 <- getUnlimitedGrowthData(res1, n=tmp.n)
	checkEquals(nrow(rd.e)-1, nrow(rd.e2) )	#one point less
}

