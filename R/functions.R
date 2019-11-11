#-----------------------------------------------------------------------------------------------------------
# All functions required for analysis
#-----------------------------------------------------------------------------------------------------------
require(dplyr)
require(combinat)
require(DEoptimR)
require(LaplacesDemon)
#-----------------------------------------------------------------------------------------------------------
checkCountsFormat <- function(counts){

	errors <- 0	

	# check is a data.frame
	if(!is.data.frame(counts)){
		warning('Raw count data must be a dataframe')
		errors <- errors + 1
		}
	# Check just one row
	if(nrow(counts)!=1){
		warning('Raw count data must comprise a single row dataframe')
		errors <- errors + 1
		}
	# check integers
	if(!is.integer(as.matrix(counts))){
		warning('Raw count data must comprise only integers')
		errors <- errors + 1
		}
	# check col names are unique
	if(sum(duplicated(names(counts)))!=0){
		warning('Column names of raw count data must be unique')
		errors <- errors + 1
		}
	# check col names are capital letters
	letters <- unique(strsplit(paste(names(counts),collapse=''),split='')[[1]])
	if(sum(!letters%in%LETTERS)!=0){
		warning("Column names of raw count data must comprise only capital letters, eg 'EFG' or 'E' ")
		errors <- errors + 1
		}
	if(errors==0)res <- 'OK'
	if(errors!=0)res <- paste(errors, 'errors')
return(res)}
#-----------------------------------------------------------------------------------------------------------
uniqueClassMaker <- function(x){
	# memory demands can be reduced hugely by retaining a multiclass that ALWAYS (every row of 'counts') have counts combined 
	# x: string of class names, comprising capital letters from A. Eg c('ABCD','A','AB','B',...)

	classes.individual <- sort(unique(strsplit(paste(x,collapse=''),split='')[[1]]))
	keep <- c()
	for(n in 1:length(classes.individual)){
		member <- grepl(classes.individual[n],x) * nchar(x)
		member[member==0] <- NA
		keep <- c(keep, which(member==min(member,na.rm=T)))
		}
return(x[unique(keep)])}
#-----------------------------------------------------------------------------------------------------------
allArrangements <- function(counts){

	if(checkCountsFormat(counts)!='OK')stop()
	classes <- uniqueClassMaker(names(counts))
	m1 <- matrix(0,1,length(classes)); colnames(m1) <- classes

	for(n in 1:length(counts)){

		# overheads
		class.counts <- counts[n]

		# separate into the classes found in 'classes'
		member <- c()
		for(c in 1:length(classes))member[c] <- (grepl(classes[c],names(class.counts)))
		class.names <- classes[member]

		# all possible arrangements of class counts, for the n th (multi)-class
		x <- t(xsimplex(length(class.names),as.integer(class.counts)))
		colnames(x) <- class.names
		m2 <- matrix(0,nrow(x),ncol(m1)); colnames(m2) <- classes
		m2[,colnames(x)] <- x

		# save memory, store as integers
		storage.mode(m2) <- 'integer'

		# all possible combinations of m1 and m2
		multi.1 <- do.call("rbind", rep(list(m1), nrow(m2)))
		multi.2 <- m2[rep(1:nrow(m2),each=nrow(m1)),]
		both <- multi.1 + multi.2

		# keep unique rows, and iterate
		m1 <- distinct(as.data.frame(both))
		}
return(m1)}
#-----------------------------------------------------------------------------------------------------------
randomArrangement <- function(counts){
	# Helper function required by GOF()
	# Generates one possible random realisation (conversion of multi-class counts into single class counts), under the assumption that the 
	# archaeologist's original multiclass assumption represents an equal probability of belonging to each of the individual classes. 
	# I.e, a count belonging to multiclass 'AB' was believed to belong to either class 'A' or 'B' with an equal probability. 

	# counts: a one row data.frame satisfying the requirements of checkCountsFormat(). I.e, handles just one dataset at a time.
	
	# ensure only one row
	if(nrow(counts)>1)warning('randomArrangement() must be handed a one-row data.frame. Only the first row has been used.')
	counts <- counts[1,]

	# scrap any with a count of zero
	counts <- counts[,counts>0,drop=F]

	# separate into uncertain and certain calls, only uncertains need sampling
	cert.counts <- counts[nchar(names(counts))==1]
	uncert.counts <- counts[nchar(names(counts))>1]

	# blank df of counts
	letters <- LETTERS[1:9]
	all.counts <- as.data.frame(matrix(0,1,length(letters))); names(all.counts) <- letters

	# add summary of random assignment to each uncertain class
	for(n in 1:length(uncert.counts)){
		rand <- sample(strsplit(names(uncert.counts[n]),'')[[1]],as.numeric(uncert.counts[n]),replace=T)
		tb <- table(rand)
		all.counts[,names(tb)] <- all.counts[,names(tb)] + as.integer(tb)
		}

	# add the certain counts
	all.counts[,names(cert.counts)] <- all.counts[,names(cert.counts)] + as.integer(cert.counts)
return(all.counts)}
#-----------------------------------------------------------------------------------------------------------
combineClasses <- function(comb, model){
	# Helper function required by GOF()
	# Aggregates individual classes to match the model classes. 
	# Column names of 'comb' must be single capital letters that occur within the column names of 'model'

	# comb: data.frame of counts, in single classes only. Eg, output by randomArrangement()
	# model: data.frame of model age class probabilities, with headers

	class.agg <- names(model)
	comb.agg <- matrix(,nrow(comb),length(class.agg))
	for(n in 1:length(class.agg)){
		v <- strsplit(class.agg[n],'')[[1]]
		comb.agg[,n] <- rowSums(comb[,v,drop=F])
		}

	# keep only the unique combinations
	if(nrow(comb.agg)>1)comb.agg <- comb.agg[!duplicated(comb.agg),,drop=F]
	comb.agg <- as.data.frame(comb.agg); colnames(comb.agg) <- class.agg
return(comb.agg)}
#--------------------------------------------------------------------------------------------------------
GOF <- function(counts, model, N = 1e+04){

	if(checkCountsFormat(counts)!='OK')stop()
	all.p <- numeric(N)
	for(n in 1:N){

		# random assignment of multi-class counts into single classes
		r <- randomArrangement(counts)

		# aggregate single classes to match the model classes
		r <- combineClasses(r,model)

		# chi-squared p-value
		x <- as.integer(r)
		p <- as.numeric(model)
		
		# remove model zero probability classes (meaningless for Chisq)
		i <- abs(p)>0.00000001
		x <- x[i]
		p <- p[i]
		
		all.p[n] <- suppressWarnings(chisq.test(x=x,p=p)$p.value)

		if(n%in%seq(0,N,by=1000))print(paste(n,'of',N,'samples completed'))
		}
	res <- formatC(mean(all.p), format = "e", digits = 1)

return(res)}
#-------------------------------------------------------------------------------------------------------
convertClassAges <- function(class.ages, multi.classes){
	# helper function to assist in computational efficiency, since class.ages can be combined where observed data permit
	# Checks various requirements first.
	
	# ensure multi.classes only comprise sequential classes
	for(n in 1:length(multi.classes)){
		char <- (strsplit(multi.classes[n],split='')[[1]])
		i <- which(LETTERS%in%char)
		if(length(i)>1){
			problem <- sum(diff(i)!=1)
			if(problem!=0)stop(paste('Multiclass',multi.classes[n],'is not valid, must comprise sequential letters'))
			}
		}

	# ensure multi.classes only comprises capital letters
	mc <- strsplit(paste(multi.classes,collapse=''),split='')[[1]]
	if(sum(!mc%in%LETTERS)!=0)stop('multi.classes must comprise only capital letters')

	# ensure class ages only comprises capital letters
	ca <- strsplit(paste(names(class.ages),collapse=''),split='')[[1]]
	if(sum(!ca%in%LETTERS)!=0)stop('names of class.ages must comprise only capital letters')

	# ensure class.ages and multi.class ages comprise the same letters
	match <- identical(sort(unique(mc)),  sort(unique(ca)))
	if(!match)stop('names of class.ages and multi.classes dont match')

	# convert the class.ages to match the multi.classes
	N <- length(multi.classes)
	new <- as.data.frame(matrix(,1,N)); names(new) <- multi.classes
	for(n in 1:N){
		i <- strsplit(multi.classes[n],split='')[[1]][1] == names(class.ages)
		new[n] <- class.ages[i]
		}
return(new)}
#-------------------------------------------------------------------------------------------------------
gammaLoglik <- function(x, class.ages, shape, mean){
	# log likelihood under any Gamma distribution 
	# converts the continuous PDF of any Gamma distribution into discrete probabilities (a multinomial probability distribution) for each age class
	# column names of 'x' and 'class.ages' must match, ie the classes (or multi-classes) must be the same

	# x: data.frame of all possible counts in each age class for which likelihoods are to be calculated and summed. Ie, the output of allArrangements()
	# class.ages: a one row data.frame of the starting age of each age class or multiclass
	# These must be sequential, such that the start of the (n+1)th class equals the end of the nth class.
	# The final age class will automatically include all counts above this. The first starting age therefore must be zero. 
	# shape, mean: gamma shape parameters

	if(!identical(names(x),names(class.ages)))stop("The class names of 'x' and 'class.names' must match")
	if(class.ages[1]!=0)stop('The first value must be zero to ensure all possible ages are encompassed')
	rate <- shape/mean
	p <- diff(c(pgamma(q=as.numeric(class.ages),shape=shape,rate=rate),1))
	loglik <- log(sum(apply(x, 1, dmultinom, prob=p, log=FALSE)))

	# in the extreme, the probabilities can shrink to almost zero (loglik becomes -Inf). This causes trouble downstream in MCMC acceptance ratios
	if(is.infinite(loglik))loglik <- -99999
return(loglik)}
#-----------------------------------------------------------------------------------------------------------
proposalFunction <- function(params, prop){
	# helper function for MCMC
	# params: vector or two values representing the 'shape' and 'mean' parameters of the gamma distribution
	# prop: hyperparameter controlling the size of the next random jump. May require tuning: smaller jumps result in a higher acceptance ratio (AR). Aim for AR between 0.2 and 0.6
	
	N <- length(params)
	jumps <- rep(prop,N) 
	moves <- rnorm(N,0,jumps)
	new.params <- abs(params + moves)
return(new.params)}
#-----------------------------------------------------------------------------------------------------------
mcmc <- function(counts, class.ages = NULL, N=30000, burn = 2000, thin = 5, prop = 0.4, plot.chain = TRUE){

	if(checkCountsFormat(counts)!='OK')stop()
	if(is.null(class.ages))class.ages <- data.frame(A=0,B=1/6,C=1/2,D=1,E=2,F=3,G=4,H=6,I=8)
	if(class.ages[1]!=0)stop('The first value must be zero to ensure all possible ages are encompassed')
	aa <- allArrangements(counts)
	ca <- convertClassAges(class.ages,names(aa))

	# initiate random Gamma parameters, between 0 and 5
	params <- runif(2,0,5)
	all.params <- matrix(,N,2)

	# house keeping progress
	progress <- round(seq(0,N,length.out=6))

	# mcmc loop
	accepted <- rep(0,N)
	for(n in 1:N){
		all.params[n,] <- params
		llik <- gammaLoglik(aa,ca,params[1],params[2])
		prop.params <- proposalFunction(params,prop)
		prop.llik <- gammaLoglik(aa,ca,prop.params[1],prop.params[2])
		ratio <- min(exp(prop.llik-llik),1)
		move <- sample(c(T,F),size=1,prob=c(ratio,1-ratio))
		if(move){
			params <- prop.params
			accepted[n] <- 1
			}
		if(n%in%progress)print(paste(n,'of',N,'samples completed'))
		}

	# remove burnin
	N.sub <- all.params[burn:N,]

	# thinning
	i <- seq(burn,N,by=thin)
	res <- all.params[i,]
	res <- as.data.frame(res); names(res) <- c('shape','mean')

	# print acceptance rate
	print(paste('acceptance rate = ',round(sum(accepted[burn:N])/nrow(N.sub),3)))

	# look at chain
	if(plot.chain){
		par(mfrow=c(2,2))
		plot(all.params[,1],type='l',xlab='Samples in chain',ylab='Gamma shape parameter',main='Entire chain')
		plot(all.params[,2],type='l',xlab='Samples in chain',ylab='Gamma mean parameter',main='Entire chain')
		hist(res[,1],breaks=20,xlab='shape',main='Burn in removed and thinned',ylab='count')
		hist(res[,2],breaks=20,xlab='mean',main='Burn in removed and thinned',ylab='count')	
		}

return(res)}
#-----------------------------------------------------------------------------------------------------------
gammaSearch <- function(counts, class.ages, trace){
	#  DEoptimR search, used by gammaMLparameters() and gammaLogMLE()

	# counts: a one row data.frame satisfying the requirements of checkCountsFormat(). I.e, handles just one dataset at a time.
	# class.ages: a one row data.frame of the starting age of each age class. See gammaLoglik()

	if(checkCountsFormat(counts)!='OK')stop()
	if(class.ages[1]!=0)stop('The first value must be zero to ensure all possible ages are encompassed')
	aa <- allArrangements(counts)
	class.ages <- convertClassAges(class.ages,names(aa))

	fn <- function(x,aa,class.ages){-gammaLoglik(aa,class.ages,x[1],x[2])}
	res <- JDEoptim(lower=c(0,0),upper=c(20,20),fn=fn,aa=aa,class.ages=class.ages,tol=1e-7,NP=10,trace=trace)
return(res)}
#-----------------------------------------------------------------------------------------------------------
gammaMLparameters <- function(counts, class.ages = NULL, trace = FALSE){
	# Search for the Maximum Likelihood Gamma parameters
	if(is.null(class.ages))class.ages <- data.frame(A=0,B=1/6,C=1/2,D=1,E=2,F=3,G=4,H=6,I=8)
	res <- gammaSearch(counts,class.ages,trace)
	x <- data.frame(shape=round(res$par[1],3),mean=round(res$par[2],3))
return(x)}
#-----------------------------------------------------------------------------------------------------------
gammaLogMLE <- function(counts, class.ages = NULL, trace = FALSE){
	# Search for the log Maximum Likelihood Estimate
	if(is.null(class.ages))class.ages <- data.frame(A=0,B=1/6,C=1/2,D=1,E=2,F=3,G=4,H=6,I=8)
	res <- gammaSearch(counts,class.ages,trace)
return(-res$value)}
#-----------------------------------------------------------------------------------------------------------
ageClassSearch <- function(counts, N, I, C){
	# Search used by both ageClassMLparameters() and ageClassLogMLE()
	# N: pop size
	# I: iterations
	# C: cooling: annealing schedule

	if(checkCountsFormat(counts)!='OK')stop()
	aa <- allArrangements(counts)
	R <- ncol(aa)

	# first iteration
	liks <- numeric(N)
	pars <- rdirichlet(N,rep(1,R))
	for(n in 1:N)liks[n] <- log(sum(apply(aa, 1, dmultinom, prob=as.numeric(pars[n,]), log=FALSE)))
	keep.pars <- pars[liks==max(liks),,drop=F]

	# subsequent iterations
	for(i in 1:I){
		liks <- numeric(N+1)
		pars <- rbind(keep.pars,rdirichlet(N,rep(1,R)+keep.pars*(C^i)))
		for(n in 1:(N+1))liks[n] <- log(sum(apply(aa, 1, dmultinom, prob=as.numeric(pars[n,]), log=FALSE)))
		keep.pars <- pars[liks==max(liks),,drop=F]	
		}

	pars <- as.data.frame(matrix(round(keep.pars,4),1,R));names(pars) <- names(aa) 
return(list(pars=pars,liks=liks))}
#-----------------------------------------------------------------------------------------------------------
ageClassMLparameters <- function(counts, N = 200, I = 40, C = 2){
	# Search for the Maximum Likelihood multinomial parameters (probabilities)
	res <- ageClassSearch(counts, N, I, C)$pars
return(res)}
#-----------------------------------------------------------------------------------------------------------
ageClassLogMLE <- function(counts, N = 200, I = 40, C = 2){
	# Search for the log Maximum Likelihood Estimate under the best fitting Multinomial
	liks <- ageClassSearch(counts, N, I, C)$liks
	res <- max(liks)
return(res)}
#-----------------------------------------------------------------------------------------------------------
