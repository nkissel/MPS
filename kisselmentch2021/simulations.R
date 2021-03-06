#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

#################################################################################
# Loading snowfall and looking for available nodes, then initializing:
#################################################################################
Sys.setenv(OMPI_MCA_btl="tcp,self")
library(snowfall)
machines <- scan("node_list.txt", what="")
machines
nmach = length(machines)
sfInit(parallel=TRUE, type='MPI',cpus=nmach, socketHosts=machines)

#################################################################################
# Loading Other necessary R Libraries:
#################################################################################
sfLibrary(glmnet)
sfLibrary(MASS)

#################################################################################
#  Support functions
#################################################################################
rmvnorm <- function(n,mean,sigma){
	ev<-eigen(sigma,symmetric=TRUE)
	Rp<-t(matrix(unlist(ev[2]),nrow=dim(sigma)[1],ncol=dim(sigma)[2],byrow=FALSE)%*%(t(matrix(unlist(ev[2]),nrow=dim(sigma)[1],ncol=dim(sigma)[2],byrow=FALSE))*sqrt(pmax(as.vector(unlist(ev[1])),0))))
	retval<-matrix(rnorm(n*ncol(sigma)),nrow=n,byrow=TRUE)%*%Rp
	retval<-sweep(retval,2,mean,"+")
	colnames(retval)<-names(mean)
	return(retval)
}

boot.select.gen <- function(resp.name, c.vars, myframe, r){
	B <- r*ncol(myframe)*1.2
	chosen1<-rep(NA,B)
	a.vars<-colnames(myframe)[which(!colnames(myframe)%in%c(c.vars,resp.name))]
	var.ind<-length(c.vars)+2
	for(i in 1:B){
		samp<-sample(1:dim(myframe)[1],floor(sqrt(dim(myframe)[1])),replace=FALSE)
		train<-myframe[samp,]
		if(is.null(c.vars)){
			lower<-lm(as.formula(paste(resp.name,"~","1",sep="")),train)
		}else{
			lower<-lm(as.formula(paste(resp.name, paste(c.vars, collapse=" + "), sep=" ~ ")),data=train)
		}
		upper<-lm(as.formula(paste(resp.name, paste(c(c.vars,a.vars), collapse=" + "), sep=" ~ ")),data=train)
		mystep<-step(lower,scope=list(lower=lower,upper=upper),direction="forward",steps=1,trace=0,k=0)
		chosen1[i]<-names(coef(mystep))[var.ind]
		if(max(table(chosen1)) >= r) {
			break
		}
	}
	rt <- table(c(chosen1, a.vars))[order(table(c(chosen1, a.vars)),decreasing=TRUE)]
	rt <- rt - 1

	return(rt)
}

boot.select.gen.fast <- function(resp.name, a.vars, c.vars, myframe, r, curr_counts){
	B <- r*ncol(myframe)*1.2
	chosen1<-rep(NA,B)
	curr_chosen <- rep(names(curr_counts), times = curr_counts)
	chosen1[seq_along(curr_chosen)] <- curr_chosen
	var.ind<-length(c.vars)+2
	for(i in (length(curr_chosen) + 1):B){
		samp<-sample(1:dim(myframe)[1],floor(sqrt(dim(myframe)[1])),replace=FALSE)
		train<-myframe[samp,]
		if(is.null(c.vars)){
			lower<-lm(as.formula(paste(resp.name,"~","1",sep="")),train)
		}else{
			lower<-lm(as.formula(paste(resp.name, paste(c.vars, collapse=" + "), sep=" ~ ")),data=train)
		}
		upper<-lm(as.formula(paste(resp.name, paste(c(c.vars,a.vars), collapse=" + "), sep=" ~ ")),data=train)
		mystep<-step(lower,scope=list(lower=lower,upper=upper),direction="forward",steps=1,trace=0,k=0)
		chosen1[i]<-names(coef(mystep))[var.ind]
		if(max(table(chosen1)) >= r) {
			break
		}
	}
	rt <- table(c(chosen1, a.vars))[order(table(c(chosen1, a.vars)),decreasing=TRUE)]
	rt <- rt - 1

	return(rt)
}

find_worstcases_chen86 <- function(r, k, P) {
	N <- r * k
	nsim <- 1e4
	checker <- TRUE
	r_prime <- r
	probs <- rep(1/k, k)
	selected <- matrix(nrow = nsim, ncol = k)
	f <- function(filler, N, probs) {
		temp <- rmultinom(N * 1.1, 1, probs)
		applied_temp <- apply(temp, 1, cumsum)
		found <- function(vect) {
			which(vect == r)
		}
		v <- apply(applied_temp, 2, FUN = found)
		col_stop <- min(unlist(v))
		return(rowSums(temp[,1:col_stop, drop = F]))
	}
	selected <- sapply(rep(NA, nsim), f, N, probs)
	if(k == 1) {
		selected <- matrix(selected, nrow = 1, ncol = nsim)
	}
	for(i in 1:(r+1)) {
		r_prime <- r - i + 1
		if(dim(selected)[1] == k) { # unsure if a change happened between R-4.0 and R-3.
			sel_prop <- mean(selected[1,] >= r_prime)
		} else {
			sel_prop <- mean(selected[,1] >= r_prime)
		}
		if(sel_prop >= P) {
			break
		}
	}
	r_prime
	return(list(r_prime = r_prime, P = sel_prop))
}

naRemover<-function(myMat){
	y<-dim(myMat)[2]
	for(i in 1:(y-1)){
		tv<-ifelse(is.na(myMat[,i]),"REP",myMat[,i])
		myRle<-rle(tv)
		toRep<-na.omit(ifelse(myRle$values!="REP"&c(myRle$values[-1],"HOLDER")!="REP",myRle$lengths,
													ifelse(myRle$values!="REP"&c(myRle$values[-1],"HOLDER")=="REP",c(myRle$lengths[-1],1)+myRle$lengths,NA)))
		myMat[,i]<-rep(myRle$values[myRle$values!="REP"],toRep)
	}
	return(myMat)
}

wc_updater_finder <- function(r, k, P) {
	# mypath <- "~/mps/check/worst_case_table.Rds"
	mypath <- "/ihome/lmentch/nkissel/mps/worst_case_table.Rds"
	worst_cases <- readRDS(mypath)
	where <- (worst_cases$r %in% r) & (worst_cases$k %in% k) & (worst_cases$P %in% P)
	is_in <- sum(where)

	if(is_in == 1) {
		return(list(r_prime = worst_cases$r_prime[where], P = P))
	} else {
		curr_rows <- nrow(worst_cases)
		t <- find_worstcases_chen86(r = r, k = k, P = P)
		worst_cases[curr_rows + 1, 1] <- r
		worst_cases[curr_rows + 1, 2] <- k
		worst_cases[curr_rows + 1, 3] <- P
		worst_cases[curr_rows + 1, 4] <- t$r_prime
		saveRDS(worst_cases, file = "./worst_case_table.Rds", version = 2)
		return(list(r_prime = t$r_prime, P = P))
	}
}

full.select.gen<-function(myframe, resp.name, depth, r,
													Pstar, condense = TRUE){
	#params for selector
	desir_prob <- Pstar
	p <- ncol(myframe)

	############################################ Data Read #############################################
	var.list<-""
	vars<-colnames(myframe)[which(colnames(myframe)!=resp.name)]

	######################################### Depth 1 select ###########################################
	select1<-boot.select.gen(resp.name,NULL,myframe,r)

	var1.list.all<-names(select1)
	var1.c.list<-select1

	p <- length(select1)
	worst_case <- wc_updater_finder(r, p, desir_prob)
	# worst_case <- find_worstcases_chen86(r, p, desir_prob)
	r_prime <- worst_case$r_prime
	var1.list <- var1.list.all[which(select1 >= r_prime )]
	c1.list <- select1[which(select1 >= r_prime )]

	inner.mat <- matrix(NA, nrow=length(var1.list),ncol=1)
	inner.mat.c <- matrix(NA, nrow=length(c1.list),ncol=1)
	inner.mat[,1] <- var1.list
	inner.mat.c[,1] <- c1.list
	print(inner.mat)

	##################################### To full depth selection #######################################
	if(depth > 1){
		for(i in 2:depth){
			running.vars <- running.cs <- NA #these are the list of ALL vars for a given depth
			running.lens<-0 #records the length of each var2.list
			for(j in 1:dim(inner.mat)[1]){ #loops through each collection of parents and finds child selections
				select2 <- boot.select.gen(resp.name,(inner.mat[j,]),myframe,r)
				var2.list.all<-names(select2)
				var2.c.list <- select2
				p <- length(select2)
				worst_case <- wc_updater_finder(r, p, desir_prob)
				# worst_case <- find_worstcases_chen86(r, p, desir_prob)
				r_prime <- worst_case$r_prime
				var2.list <- var2.list.all[which(select2 >= r_prime )]
				c2.list <- select2[which(select2 >= r_prime )]

				running.vars <- c(running.vars,var2.list)
				running.cs <- c(running.cs, c2.list)
				running.lens <- c(running.lens,length(var2.list))
			}
			running.vars <- running.vars[-1]
			running.cs <- running.cs[-1]
			running.lens <- running.lens[-1]
			inner.mat2 <- matrix(NA,nrow=length(running.vars),ncol=i)
			inner.mat.c2 <- matrix(NA,nrow=length(running.cs),ncol=i)
			ind <- cumsum(running.lens) + 1
			inner.mat2[c(1, ind[-length(ind)]), 1:(i-1)] <- inner.mat
			inner.mat2[, i] <- running.vars
			inner.mat <- naRemover(inner.mat2)
			inner.mat.c2[c(1, ind[-length(ind)]), 1:(i-1)] <- inner.mat.c
			inner.mat.c2[, i] <- running.cs
			inner.mat.c <- inner.mat.c2
			if(condense){
				inner.mat3<-t(apply(inner.mat,1,sort))
				inner.mat.c<-inner.mat.c[!duplicated(inner.mat3),,drop=FALSE]
				inner.mat<-inner.mat[!duplicated(inner.mat3),,drop=FALSE]
			}
			print(inner.mat)
		}
	}
	return(list(mps = inner.mat, counts = inner.mat.c))
}

full.select.gen.fast<-function(myframe, resp.name, depth, r, Pstar, condense = TRUE){
	#params for selector
	desir_prob <- Pstar
	p <- ncol(myframe) - 1

	############################################ Data Read #############################################
	var.list<-""
	vars <- colnames(myframe)[which(colnames(myframe)!=resp.name)]

	######################################### Depth 1 select ###########################################
	intermed_r <- 9
	intermed_P <- 0.95
	select1.temp <- boot.select.gen.fast(resp.name, vars, NULL, myframe, intermed_r, NULL)
	p <- length(select1.temp)
	var1.list.all <- names(select1.temp)
	var1.c.list <- select1.temp

	worst_case <- wc_updater_finder(intermed_r, p, intermed_P)
	r_prime <- worst_case$r_prime
	var1.list.temp <- var1.list.all[which(select1.temp >= r_prime )]
	c1.list.temp <- select1.temp[which(select1.temp >= r_prime )]

	select1.temp2 <- boot.select.gen.fast(resp.name, var1.list.temp, NULL, myframe, r, c1.list.temp)
	select1.temp2 <- select1.temp2[order(names(select1.temp2))]
	select1 <- select1.temp2
	select1 <- rev(sort(select1))
	worst_case <- wc_updater_finder(r, length(select1), desir_prob)
	r_prime <- worst_case$r_prime

	var1.list.all<-names(select1)
	var1.c.list <- select1
	var1.list <- var1.list.all[which(select1 >= r_prime )]
	c1.list <- select1[which(select1 >= r_prime )]

	# print(var1.c.list)



	inner.mat <- matrix(NA, nrow=length(var1.list),ncol=1)
	inner.mat.c <- matrix(NA, nrow=length(c1.list),ncol=1)
	inner.mat[,1] <- var1.list
	inner.mat.c[,1] <- c1.list
	print(inner.mat)

	##################################### To full depth selection #######################################
	if(depth > 1){
		for(i in 2:depth){
			running.vars <- running.cs <- NA #these are the list of ALL vars for a given depth
			running.lens<-0 #records the length of each var2.list
			for(j in 1:dim(inner.mat)[1]){ #loops through each collection of parents and finds child selections
				c.vars <- inner.mat[j,]
				a.vars <- colnames(myframe)[which(colnames(myframe)!=c(c.vars,resp.name))]
				select2.temp <- boot.select.gen.fast(resp.name, a.vars, c.vars, myframe, intermed_r, NULL)
				p <- length(select2.temp)
				var2.list.all <- names(select2.temp)
				var2.c.list <- select2.temp

				worst_case <- wc_updater_finder(intermed_r, p, intermed_P)
				r_prime <- worst_case$r_prime
				var2.list.temp <- var2.list.all[which(select2.temp >= r_prime )]
				c2.list.temp <- select2.temp[which(select2.temp >= r_prime )]

				select2.temp2 <- boot.select.gen.fast(resp.name, var2.list.temp, c.vars, myframe, r, c2.list.temp)
				select2.temp2 <- select2.temp2[order(names(select2.temp2))]
				select2 <- select2.temp2
				select2 <- rev(sort(select2))
				worst_case <- wc_updater_finder(r, length(select2), desir_prob)
				r_prime <- worst_case$r_prime

				var2.list.all<-names(select2)
				var2.c.list <- select2
				var2.list <- var2.list.all[which(select2 >= r_prime )]
				c2.list <- select2[which(select2 >= r_prime )]

				# print(var2.c.list)


				running.vars <- c(running.vars,var2.list)
				running.cs <- c(running.cs, c2.list)
				running.lens <- c(running.lens,length(var2.list))
			}
			running.vars <- running.vars[-1]
			running.cs <- running.cs[-1]
			running.lens <- running.lens[-1]
			inner.mat2 <- matrix(NA,nrow=length(running.vars),ncol=i)
			inner.mat.c2 <- matrix(NA,nrow=length(running.cs),ncol=i)
			ind <- cumsum(running.lens) + 1
			inner.mat2[c(1, ind[-length(ind)]), 1:(i-1)] <- inner.mat
			inner.mat2[, i] <- running.vars
			inner.mat <- naRemover(inner.mat2)
			inner.mat.c2[c(1, ind[-length(ind)]), 1:(i-1)] <- inner.mat.c
			inner.mat.c2[, i] <- running.cs
			inner.mat.c <- inner.mat.c2
			if(condense){
				inner.mat3<-t(apply(inner.mat,1,sort))
				inner.mat.c<-inner.mat.c[!duplicated(inner.mat3),,drop=FALSE]
				inner.mat<-inner.mat[!duplicated(inner.mat3),,drop=FALSE]
			}
			print(inner.mat)
		}
	}
	return(list(mps = inner.mat, counts = inner.mat.c))
}


CleanErrMat <- function(err.mat) {
	# This is a subroutine called in the main function "CVPerm"
	# It removes identical columns in the cross-validated loss matrix
	n2 <- nrow(err.mat)
	M <- ncol(err.mat)
	err.mean <- apply(err.mat, 2, mean)
	err.mat.center <- err.mat - matrix(err.mean, nrow = n2, ncol = M, byrow = T)

	res <- list()
	list.length <- 0
	remaining.index <- 1:M

	while (length(remaining.index) > 0) {
		m <- remaining.index[1]
		err.diff.center <- err.mat.center[,m] - err.mat.center
		err.mean.diff <- err.mean[m] - err.mean
		sd.vec <- apply(err.diff.center, 2, sd)
		J.1 <- which(sd.vec==0 & err.mean.diff > 0)
		J.2 <- which(sd.vec==0 & err.mean.diff == 0)
		J.3 <- which(sd.vec==0 & err.mean.diff < 0)
		if (length(J.1) == 0) {
			list.length <- list.length+1
			res[[list.length]] <- J.2
		}
		remaining.index <- setdiff(remaining.index, union(J.2,J.3))
	}
	return(list(err.mat.c = err.mat[,sapply(res,function(x){x[1]})], ind=res))
}

CVPerm <- function(err.mat, B = 200, screen=T, alpha.s=0.005) {
	# This is the main function of cross-validation with confidence.
	# Requires subroutine "CleanErrMat"
	# Input: "err.mat" is a matrix where each row corresponds to a data point
	#            and each column corresponds to a candidate model.
	#            err.mat[i,j] records the cross-validated loss evaluated at
	#            the i'th data point and the j'th candidate model.
	#            For example, in least square linear regression this is the
	#            cross-validated squared residual
	#        "B" is the number of bootstrap samples
	#        "screen" is an indicator if pre-screening is used for the test
	#        "alpha.s" is the threshold used in pre-screening
	clean.res <- CleanErrMat(err.mat)
	err.mat.c <- clean.res$err.mat.c
	err.mat.c.ind <- clean.res$ind
	n2 <- nrow(err.mat.c)
	M <- ncol(err.mat.c)
	err.mean <- apply(err.mat.c, 2, mean)
	err.mat.center <- err.mat.c - matrix(err.mean, nrow = n2, ncol = M, byrow = T)
	sign.mat <- matrix(rnorm(n2*B),ncol = B)
	norm.quantile <- qnorm(1-alpha.s/(M-1))
	screen.th <- ifelse(norm.quantile^2>=n2,Inf,norm.quantile / sqrt(1-norm.quantile^2/n2))
	sgmb.p.val <- sapply(1:M, function(m){
		err.diff.center <- err.mat.center[,m] - err.mat.center[, setdiff(1:M, m),drop=F]
		err.mean.diff <- err.mean[m] - err.mean[setdiff(1:M, m)]
		sd.vec <- apply(err.diff.center, 2, sd)
		err.mean.diff.scale <- err.mean.diff / sd.vec
		test.stat <- sqrt(n2)*max(err.mean.diff.scale)
		if (test.stat >= screen.th) {
			return(alpha.s)
		}
		if (screen) {
			J.screen <- which(sqrt(n2) * err.mean.diff.scale > -2*screen.th)
			if (length(J.screen)==0) {return(1-alpha.s)}
		} else {
			J.screen <- 1:(M-1)
		}
		err.diff.center.scale <- err.diff.center[,J.screen] /
			matrix(sd.vec[J.screen], nrow = n2, ncol = length(J.screen), byrow=T)
		test.stat.vec <- sapply(1:B,
														function(ib){max(apply(err.diff.center.scale*sign.mat[,ib],2,mean))})
		return(mean(test.stat.vec>max(err.mean.diff.scale[J.screen])))
	})
	res <- rep(0,ncol(err.mat))
	for (i in 1:length(err.mat.c.ind)) {
		res[err.mat.c.ind[[i]]] <- sgmb.p.val[i]
	}
	return(res)
}

CVC <- function(X, Y, type = "stepwise", max.steps=NULL, lambda=NULL, nlambda=NULL, ols.models=NULL, n.fold = 5, B = 200, screen=T, alpha.s = 0.005) {
	# This is the main wrapper function that applies CVC to various linear regression methods
	# Input: "X", "Y" are the covariate and response, in matrix form
	#        "type": "OLS", "stepwise", "lasso", "lars"
	#            if type=="OLS", then a list of candidate models "ols.models" shall be provided,
	#            see "low_dim_example.R" for example
	#        "max.steps" is the maximum number of steps used in lars package, which is used to implement forward stepwise
	#        "lambda" is the vector of lambda values used in lasso, implemented by "glmnet" package
	#        "nlambda": instead of providing a sequence of lambda, one can also just provide a total number of lambda values,
	#            the sequence of lambda will then be chosen by glmnet.
	#        "ols.models": the list of candidate ols models when type=="OLS", see "low_dim_example.R" for example
	#        "n.fold" is the number of folds
	#        "B", "screen" are parameters passed to the main CVC function "CVPerm"
	n <- nrow(X)
	p <- ncol(X)
	fold.size <- floor(n/n.fold)
	n <- n.fold * fold.size
	fold.ind <- matrix(sample.int(n),nrow = n.fold)
	if (type=="ols") {
		M <- length(ols.models)
		cv.res <- list(test.err=matrix(0,ncol=M,nrow=0), beta=array(0,c(p,M,n.fold)))
		for (i.fold in 1:n.fold) {
			test.ind <- fold.ind[i.fold,]
			train.ind <- setdiff(1:n, test.ind)
			X.train <- X[train.ind,]
			Y.train <- Y[train.ind]
			X.test <- X[test.ind,]
			Y.test <- Y[test.ind]
			beta.est <- sapply(ols.models,
												 function(x){est.coef <- rep(0,p)
												 if (length(x)>0) {
												 	est.coef[x] <- lm(Y.train ~ X.train[,x]-1)$coef
												 	est.coef[is.na(est.coef)] <- 0
												 }
												 return(est.coef)})
			predict.vals <- X.test %*% beta.est
			test.err <- predict.vals - Y.test
			cv.res$test.err <- rbind(cv.res$test.err, test.err)
			cv.res$beta[,,i.fold] <- beta.est
		}
	}

	test.err.1 <- cv.res$test.err[1:fold.size,]
	p.vals.1 <- CVPerm(test.err.1^2, B = B, screen=screen, alpha.s=alpha.s)
	test.err.mean.1 <- apply(test.err.1, 2, mean)
	p.vals.c <- CVPerm(cv.res$test.err^2, B = B, screen=screen, alpha.s=alpha.s)
	return(list(test.err.c = cv.res$test.err, p.vals.c = p.vals.c, test.err.1 = test.err.1, p.vals.1 = p.vals.1,
							beta.1 = cv.res$beta[,,1], beta = cv.res$beta))
}

mps_cvc_func <- function(x, y, ntrue, p) {
	picked_mods <- list()
	for(i in 1:ntrue) {
		sel_mods <- list()
		for(j in seq_along(picked_mods)) {
			a_vars_inds <- (1:p)[!((1:p) %in% picked_mods[[j]])]
			my_ols.models_1 <- matrix(nrow = length(a_vars_inds), ncol = i)
			my_ols.models_1[,1:(i-1)] <- rep(picked_mods[[j]], each = length(a_vars_inds))
			my_ols.models_1[,i] <- a_vars_inds
			my_ols.models_1 <- split(t(my_ols.models_1), rep(1:nrow(my_ols.models_1), each = ncol(my_ols.models_1)))
			mycvc_temp <- CVC(x, as.matrix(y), type = "ols", ols.models = my_ols.models_1, n.fold = as.numeric(args[9]), B = 100)
			picked_temp <- which(mycvc_temp$p.vals.c > 0.05)
			sel_mods <- c(sel_mods, my_ols.models_1[picked_temp])
			names(sel_mods) <- NULL
		}
		if(i == 1) {
			my_ols.models_1 <- as.list(1:p)
			mycvc_temp <- CVC(x, as.matrix(y), type = "ols", ols.models = my_ols.models_1, n.fold = as.numeric(args[9]), B = 100)
			picked_temp <- which(mycvc_temp$p.vals.c > 0.05)
			sel_mods <- c(sel_mods, my_ols.models_1[picked_temp])
			names(sel_mods) <- NULL
		}
		picked_mods <- sel_mods
		picked_mods_mat <- do.call(rbind, picked_mods)
		picked_mods <- picked_mods[!duplicated(t(apply(picked_mods_mat, 1, sort)))]
	}
	return(picked_mods)
}


##################################################################################################################################################################
#  Main function
##################################################################################################################################################################
distributed_func <- function(){

	go_fast <- as.logical(args[8])
	do_cvc <- TRUE

	tfmat0<-tfmat1<-tfmat2<-tfmat3<-tfmat4<-tfmat5<-NA

	B <- 100
	p <- as.numeric(args[1])
	n <- as.numeric(args[2])
	r <- as.numeric(args[3])
	SNR <- as.numeric(args[4])

	Pstar <- as.numeric(args[7])
	ntrue <- 5

	mysigma<-diag(p)
	noise_type <- "toep"

	rho <- as.numeric(args[5])
	if(noise_type == "toep") {
		for(j in 1:p){
			for(k in 1:p){
				mysigma[j,k]<-rho^abs(j-k)
			}
		}
	}

	n1<-rmvnorm(n=n,mean=rep(0,p),sigma=mysigma)
	colnames(n1)<-paste("x",1:p,sep="")
	# coefs <- rep(1, ntrue)
	coefs <- seq(10, 0.5, len = ntrue)
	s.covar <- 1:ntrue
	# s.covar <- round(seq(1, p, len = ntrue))
	betas <- rep(0,p)
	betas[s.covar] <- coefs
	betas <- matrix(betas,ncol=1)
	signal <- (n1%*%betas)
	eps <- rnorm(n,0,sqrt(var(signal)/SNR))
	y <- signal+eps

	n1.test <- rmvnorm(n = 100000, mean = rep(0, p), sigma = mysigma)
	colnames(n1.test)<-paste("x",1:p,sep="")
	signal.test <- (n1.test %*% betas)
	eps.test <- rnorm(100000, 0, sqrt(var(signal) / SNR))
	y.test <- signal.test + eps.test

	model<-glmnet(x = n1, y = y)
	mycond<-which(colSums(ifelse(as.matrix(coef(model)[-1,])==0,0,1))<=ntrue)
	mysnum<-names(which.max(mycond))
	red.coefs<-coef(model)[-1,which(colnames(coef(model))==mysnum)]
	myranked<-rank(-abs(red.coefs))
	myrankedOr<-myranked[order(myranked)]
	lasso.sel<-	names(myrankedOr)[1:ntrue]
	tfmat0<-lasso.sel
	pmat <- predict(model, n1.test)
	tfacc0 <- mean((pmat[,colnames(pmat) == mysnum] - y.test)^2)

	mycv <- cv.glmnet(x = n1, y = y)
	lasso_cved <- glmnet(x = n1, y = y, lambda = mycv$lambda.min)
	red.coefs <- coef(lasso_cved)[-1,]
	lasso.sel <- names(which(red.coefs != 0))
	tfmat0_cv <- lasso.sel
	tfacc0_cv <- mean((predict(lasso_cved, n1.test) - y.test)^2)


	mymat<-matrix(NA,nrow=B,ncol=floor(sqrt(.8*p)))
	for(j in 1:B){
		samp<-sample(1:n, floor(n/2))
		train.n1<-n1[samp,]
		train.y<-y[samp]
		model<-glmnet(x=train.n1,y=train.y)
		mysnum<-names(which.max(colSums(coef(model)[-1,which(colSums(ifelse(as.matrix(coef(model)[-1,])==0,0,1))<=sqrt(.8*p)),drop=FALSE]^2)))
		select.vect<-names(na.omit(ifelse(coef(model)[-1,which(colnames(coef(model)[-1,])==mysnum)]==0,NA,1)))
		if(length(select.vect)>1){
			mymat[j,1:length(select.vect)]<-select.vect
		}else{
			mymat[j,1]<-NA
		}
	}
	tb<-table(mymat)
	myrank<-rank(-tb[order(tb,decreasing=TRUE)],ties.method="random")
	myrankOrdered<-myrank[order(myrank,decreasing=FALSE)] #ordered by rank
	tfmat1<-names(myrankOrdered)[1:ntrue]

	mymat<-matrix(NA,nrow=B,ncol=floor(sqrt(.8*p)))
	for(j in 1:B){
		samp<-sample(1:n, floor(sqrt(n)))
		train.n1<-n1[samp,]
		train.y<-y[samp]
		model<-glmnet(x=train.n1,y=train.y)
		mysnum<-names(which.max(colSums(coef(model)[-1,which(colSums(ifelse(as.matrix(coef(model)[-1,])==0,0,1))<=sqrt(.8*p)),drop=FALSE]^2)))
		select.vect<-names(na.omit(ifelse(coef(model)[-1,which(colnames(coef(model)[-1,])==mysnum)]==0,NA,1)))
		if(length(select.vect)>1){
			mymat[j,1:length(select.vect)]<-select.vect
		}else{
			mymat[j,1]<-NA
		}
	}
	tb<-table(mymat)
	myrank<-rank(-tb[order(tb,decreasing=TRUE)],ties.method="random")
	myrankOrdered<-myrank[order(myrank,decreasing=FALSE)] #ordered by rank
	tfmat1.sub<-names(myrankOrdered)[1:ntrue]

	mydata<-data.frame(n1,y)
	lower<-lm(as.formula(paste("y","~","1",sep="")),mydata)
	upper<-lm(as.formula(paste("y","~",".",sep="")),mydata)
	mystep<-step(lower,scope=list(lower=lower,upper=upper),direction="forward",steps=ntrue,trace=0,k=0)
	tfmat2<-names(coef(mystep))[-1]
	#
	mypicks<-c()
	for(j in 1:B){
		samp<-sample(1:n,n,replace=TRUE)
		train<-mydata[samp,]
		lower<-lm(as.formula(paste("y","~","1",sep="")),train)
		upper<-lm(as.formula(paste("y","~",".",sep="")),train)
		mystep<-step(lower,scope=list(lower=lower,upper=upper),direction="forward",steps=ntrue,trace=0,k=0)
		steps<-names(coef(mystep))[-1]
		mypicks<-c(mypicks,steps)
	}
	tfmat3<-names(table(mypicks)[order(table(mypicks),decreasing=T)][1:ntrue])


	mypicks<-c()
	for(j in 1:B){
		samp<-sample(1:(n),floor(sqrt(n)),replace=FALSE)
		train<-mydata[samp,]
		lower<-lm(as.formula(paste("y","~","1",sep="")),train)
		upper<-lm(as.formula(paste("y","~",".",sep="")),train)
		mystep<-step(lower,scope=list(lower=lower,upper=upper),direction="forward",steps=ntrue,trace=0,k=0)
		steps<-names(coef(mystep))[-1]
		mypicks<-c(mypicks,steps)
	}
	tfmat3.sub<-names(table(mypicks)[order(table(mypicks),decreasing=T)][1:ntrue])

	full_cvc_res <- mps_cvc_res <- cvc_time <- mps_cvc_time <- NULL
	if(do_cvc) {
		mymat <- combn(1:p, ntrue)
		my_ols.models <- split(mymat, rep(1:ncol(mymat), each = nrow(mymat)))
		pt <- proc.time()
		mycvc <- CVC(n1, as.matrix(y), type = "ols", ols.models = my_ols.models, n.fold = as.numeric(args[9]), B = 100)
		cvc_time <- proc.time() - pt
		full_cvc_res <- matrix(paste0("x", do.call(rbind, my_ols.models[which(mycvc$p.vals.c > 0.05)])), ncol = ntrue)

		pt <- proc.time()
		picked_mods_ls <- mps_cvc_func(n1, as.matrix(y), ntrue, p)
		mps_cvc_time <- proc.time() - pt
		mps_cvc_res <- matrix(paste0("x", do.call(rbind, picked_mods_ls)), ncol = ntrue)
	}


	mse <- function(y,y.pred){return(mean((y-y.pred)^2))}
	pt <- proc.time()
	if(go_fast) {
		full_sel <- full.select.gen.fast(myframe = data.frame(y,n1), resp.name = "y", depth = ntrue, r = r,
																		 Pstar = Pstar, condense = TRUE)
	} else {
		full_sel <- full.select.gen(myframe = data.frame(y,n1), resp.name = "y", depth = ntrue, r = r,
																Pstar = Pstar, condense = TRUE)
	}
	mps_time <- proc.time() - pt
	myc <- full_sel$mps
	myc_counts <- full_sel$counts
	tfmat6 <- myc[1,]
	tfmat7 <- myc


	tfacc1<-tfacc2<-tfacc3<-tfacc4<-tfacc5<-tfacc6<-tfacc7<-tfacc_cvc<-tfacc_mps_cvc<-NA

	rte0 <- NULL

	mod<-lm(as.formula(paste("y", paste(c(na.omit(tfmat1)), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
	pred<-predict(mod,data.frame(n1.test))
	tfacc1<-mean((pred-y.test)^2)
	lmcoefs1<-(coefficients(mod)[-1])
	lmbetas1<-rep(0,p)
	lmbetas1[as.numeric(gsub("[a-zA-Z ]", "", na.omit(tfmat1)))]<-lmcoefs1
	lmbetas1<-matrix(lmbetas1,ncol=1)
	rte1<-(t(lmbetas1-betas)%*%(mysigma)%*%(lmbetas1-betas)+(var(signal)/SNR))/(var(signal)/SNR)

	mod<-lm(as.formula(paste("y", paste(c(na.omit(tfmat1.sub)), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
	pred<-predict(mod,data.frame(n1.test))
	tfacc1.sub<-mean((pred-y.test)^2)
	lmcoefs1<-(coefficients(mod)[-1])
	lmbetas1<-rep(0,p)
	lmbetas1[as.numeric(gsub("[a-zA-Z ]", "", na.omit(tfmat1.sub)))]<-lmcoefs1
	lmbetas1<-matrix(lmbetas1,ncol=1)
	rte1.sub<-(t(lmbetas1-betas)%*%(mysigma)%*%(lmbetas1-betas)+(var(signal)/SNR))/(var(signal)/SNR)

	mod<-lm(as.formula(paste("y", paste(c(na.omit(tfmat2)), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
	pred<-predict(mod,data.frame(n1.test))
	tfacc2<-mean((pred-y.test)^2)
	lmcoefs1<-(coefficients(mod)[-1])
	lmbetas1<-rep(0,p)
	lmbetas1[as.numeric(gsub("[a-zA-Z ]", "", na.omit(tfmat2)))]<-lmcoefs1
	lmbetas1<-matrix(lmbetas1,ncol=1)
	rte2<-(t(lmbetas1-betas)%*%(mysigma)%*%(lmbetas1-betas)+(var(signal)/SNR))/(var(signal)/SNR)

	mod<-lm(as.formula(paste("y", paste(c(na.omit(tfmat3)), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
	pred<-predict(mod,data.frame(n1.test))
	tfacc3<-mean((pred-y.test)^2)
	lmcoefs1<-(coefficients(mod)[-1])
	lmbetas1<-rep(0,p)
	lmbetas1[as.numeric(gsub("[a-zA-Z ]", "", na.omit(tfmat3)))]<-lmcoefs1
	lmbetas1<-matrix(lmbetas1,ncol=1)
	rte3<-(t(lmbetas1-betas)%*%(mysigma)%*%(lmbetas1-betas)+(var(signal)/SNR))/(var(signal)/SNR)

	mod<-lm(as.formula(paste("y", paste(c(na.omit(tfmat3.sub)), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
	pred<-predict(mod,data.frame(n1.test))
	tfacc3.sub<-mean((pred-y.test)^2)
	lmcoefs1<-(coefficients(mod)[-1])
	lmbetas1<-rep(0,p)
	lmbetas1[as.numeric(gsub("[a-zA-Z ]", "", na.omit(tfmat3.sub)))]<-lmcoefs1
	lmbetas1<-matrix(lmbetas1,ncol=1)
	rte3.sub<-(t(lmbetas1-betas)%*%(mysigma)%*%(lmbetas1-betas)+(var(signal)/SNR))/(var(signal)/SNR)

	rte_cvc<-NULL
	if(do_cvc){
		for(i in 1:dim(full_cvc_res)[1]){
			mod<-lm(as.formula(paste("y", paste(c(na.omit(full_cvc_res[i,])), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
			pred<-predict(mod,data.frame(n1.test))
			tfacc_cvc[i]<-mean((pred-y.test)^2)
			lmcoefs1<-(coefficients(mod)[-1])
			lmbetas1<-rep(0,p)
			lmbetas1[as.numeric(gsub("[a-zA-Z ]", "", na.omit(full_cvc_res[i,])))]<-lmcoefs1
			lmbetas1<-matrix(lmbetas1,ncol=1)
			rte_cvc[i]<-(t(lmbetas1-betas)%*%(mysigma)%*%(lmbetas1-betas)+(var(signal)/SNR))/(var(signal)/SNR)
		}
	}

	rte_mps_cvc<-NULL
	if(do_cvc){
		for(i in 1:dim(mps_cvc_res)[1]){
			mod<-lm(as.formula(paste("y", paste(c(na.omit(mps_cvc_res[i,])), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
			pred<-predict(mod,data.frame(n1.test))
			tfacc_mps_cvc[i]<-mean((pred-y.test)^2)
			lmcoefs1<-(coefficients(mod)[-1])
			lmbetas1<-rep(0,p)
			lmbetas1[as.numeric(gsub("[a-zA-Z ]", "", na.omit(mps_cvc_res[i,])))]<-lmcoefs1
			lmbetas1<-matrix(lmbetas1,ncol=1)
			rte_mps_cvc[i]<-(t(lmbetas1-betas)%*%(mysigma)%*%(lmbetas1-betas)+(var(signal)/SNR))/(var(signal)/SNR)
		}
	}


	mod<-lm(as.formula(paste("y", paste(c(na.omit(tfmat6)), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
	pred<-predict(mod,data.frame(n1.test))
	tfacc6<-mean((pred-y.test)^2)
	lmcoefs1<-(coefficients(mod)[-1])
	lmbetas1<-rep(0,p)
	lmbetas1[as.numeric(gsub("[a-zA-Z ]", "", na.omit(tfmat6)))]<-lmcoefs1
	lmbetas1<-matrix(lmbetas1,ncol=1)
	rte6<-(t(lmbetas1-betas)%*%(mysigma)%*%(lmbetas1-betas)+(var(signal)/SNR))/(var(signal)/SNR)

	rte7<-NULL
	for(i in 1:dim(tfmat7)[1]){
		mod<-lm(as.formula(paste("y", paste(c(na.omit(tfmat7[i,])), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
		pred<-predict(mod,data.frame(n1.test))
		tfacc7[i]<-mean((pred-y.test)^2)
		lmcoefs1<-(coefficients(mod)[-1])
		lmbetas1<-rep(0,p)
		lmbetas1[as.numeric(gsub("[a-zA-Z ]", "", na.omit(tfmat7[i,])))]<-lmcoefs1
		lmbetas1<-matrix(lmbetas1,ncol=1)
		rte7[i]<-(t(lmbetas1-betas)%*%(mysigma)%*%(lmbetas1-betas)+(var(signal)/SNR))/(var(signal)/SNR)
	}

	mod<-lm(as.formula(paste("y", paste(c(colnames(n1)[s.covar]), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
	pred<-predict(mod,data.frame(n1.test))
	tacc<-mean((pred-y.test)^2)
	pred_true <- n1.test %*% betas
	tacc_true <- mean((pred_true - y.test)^2)
	lmcoefs1<-(coefficients(mod)[-1])
	lmbetas1<-rep(0,p)
	lmbetas1[as.numeric(gsub("[a-zA-Z ]", "", colnames(n1)[s.covar]))]<-lmcoefs1
	lmbetas1<-matrix(lmbetas1,ncol=1)
	trte<-(t(lmbetas1-betas)%*%(mysigma)%*%(lmbetas1-betas)+(var(signal)/SNR))/(var(signal)/SNR)

	return(list(truecov = colnames(n1)[s.covar], lasso = tfmat0, lasso_cv = tfmat0_cv,
							ss=tfmat1, ss.sub = tfmat1.sub, fs = tfmat2, ssfs = tfmat3,
							ssfs.sub = tfmat3.sub, cvc = full_cvc_res, mps_cvc = mps_cvc_res,
							fss.sub = tfmat6, pi.sub = tfmat7, pi.counts = myc_counts,
							trueacc = tacc_true, trueacc_emp = tacc, lasso_acc = tfacc0, lasso_acc_cv = tfacc0_cv,
							ss_acc = tfacc1, ss_acc.sub = tfacc1.sub, fs_acc = tfacc2, ssfs_acc = tfacc3,
							ssfs_acc.sub = tfacc3.sub, tfacc_cvc = tfacc_cvc, tfacc_mps_cvc = tfacc_mps_cvc,
							fss_acc.sub = tfacc6, pi_acc.sub = tfacc7, lasso_rte = rte0,
							ss_rte = rte1, ss_rte.sub=rte1.sub, fs_rte = rte2, ssfs_rte = rte3,
							ssfs_rte.sub = rte3.sub, rte_cvc = rte_cvc, rte_mps_cvc = rte_mps_cvc, fss_rte.sub = rte6,
							pi_rte.sub = rte7, t_rte=trte,
							cvc_time = cvc_time, mps_time = mps_time, mps_cvc_time = mps_cvc_time))
}

#################################################################################
#  Exporting the data and functions needed by the worker nodes:
#################################################################################
sfExport("args")
sfExport("rmvnorm")
sfExport("boot.select.gen.fast")
sfExport("boot.select.gen")
sfExport("naRemover")
sfExport("full.select.gen.fast")
sfExport("full.select.gen")
sfExport("find_worstcases_chen86")
sfExport("wc_updater_finder")
sfExport("distributed_func")
sfExport("mps_cvc_func")
sfExport("CVPerm")
sfExport("CVC")
sfExport("CleanErrMat")


#################################################################################
#  Creating a wrapper function
#################################################################################
wrapper <- function(nsim){
	result <- distributed_func()
	return(result)
}


#################################################################################
#  Running the function nSim times:
#################################################################################
nSim <- as.numeric(args[6])
FinalResult <- sfLapply(1:nSim, wrapper)


#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(FinalResult, file = paste0('tibs_b3_nsim', args[6],
																'_p', args[1],'_n', args[2],'_r',
																args[3],'_SNR',args[4],'_rho',args[5],"_ps",args[7],'_fast',args[8],'.RData'), version = 2)


#################################################################################
#  Close all connections:
#################################################################################
sfStop()
