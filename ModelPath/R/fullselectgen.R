
f<-function(y,y.pred){return(mean((y-y.pred)^2))} #mean squared error

######################################################################################################
######################################################################################################
##################################  Full Selection Procedure Gen  ####################################
######################################################################################################
######################################################################################################
#' Perform MPS
#'
#' Function that performs Model Path Selection (MPS)
#' @param myframe data.frame holding explanatory and response data
#' @param resp.name matrix of MPS data
#' @param depth maximum depth of paths
#' @param r maximum cell count
#' @param Pstar inclusion probability bound
#' @param f loss function
#' @param type 'reg' or 'class' for regression and classification. Not needed if \code{f} is specified.
#' @param model name of modeling method; must be name of the function called in R
#' @param fun.args list of arguments that are passed to modeling function
#' @param pred.args list of arguments that are passed to model predict function
#' @param condense logical for whether or not paths should be merged
#' @param save.table whether you want to save a table for faster computation
#' @param table.path file path to where table is saved; if not specified, will default to working directory
#' @param trace T/F for whether intermediary output is printed
#' @keywords matrix
#' @examples set.seed(200)
#' library(MASS)
#' n <- 1000
#' p <- 10
#' x <- mvrnorm(n, rep(0, p), diag(p))
#' signal <- rowSums(x[, 1:3])
#' noise <- rnorm(n, 0, 1/4)
#' y <- signal+noise
#' mydata <- data.frame(x, y)
#' mps1 <- full.select.gen(myframe = mydata, resp.name = 'y', depth = 3, r = 50,
#'                         model = 'lm', condense = FALSE)
#' mps1
#' build.tree(mps1) #graph
#'
#' #merged paths
#' mps2<-full.select.gen(myframe=mydata,resp.name='y',depth=3,
#'    r=100,model='lm',condense=TRUE)
#' mps2
#' build.tree(mps2) #graph
#' @export

full.select.gen <- function(
		myframe, resp.name, depth, r, Pstar=0.95, f, type, model, fun.args,
		pred.args, condense = F, save.table = F, table.path = NULL, trace=F) {
	####################################### input info and stops #######################################
	if (missing(depth)) depth <- 3
	if (missing(myframe)) stop("Need data.frame")
	if (missing(type)) type <- "reg"
	if ((type != "reg") && (type != "class")) {
		stop("Need to specify type as \"class\" or \"reg\"")
	}
	if (missing(r)) r <- 100
	if (missing(f)) {
		f <- function(y, y.pred) return(mean((y - y.pred)^2)) #mean squared error
		if(type=="class"){
			f <- function(y, y.pred) return(sum(y != y.pred) / length(y)) #misclassification rate
		}
	}
	if (!is.function(f)) stop("f must be a function")
	if (missing(model)) model <- "lm"
	if (missing(fun.args)) fun.args <- list()
	if (missing(pred.args)) pred.args <- list()

	############################################ Data Read #############################################
	var.list <- ""
	vars <- colnames(myframe)[which(colnames(myframe) != resp.name)]

	######################################### Depth 1 select ###########################################
	select1 <- resample_select(
		resp.name, NULL, myframe, type, r, f, model, fun.args, pred.args)
	var1.list.all <- names(select1)
	var1.c.list <- select1

	p <- length(select1)
	desir_prob <- Pstar
	worst_case <- wc_updater_finder(
		r, p, desir_prob, save.table=save.table, table.path=table.path)
	# worst_case <- find_worstcases_chen86(r, p, desir_prob)
	r_prime <- worst_case$r_prime
	var1.list <- var1.list.all[which(select1 >= r_prime )]
	c1.list <- select1[which(select1 >= r_prime )]

	inner.mat <- matrix(NA, nrow=length(var1.list),ncol=1)
	inner.mat.c <- matrix(NA, nrow=length(c1.list),ncol=1)
	inner.mat[,1] <- var1.list
	inner.mat.c[,1] <- c1.list
	if(trace) print(inner.mat)

	##################################### To full depth selection #######################################
	if (depth > 1){
		for (i in 2:depth) {
			running.vars <- running.cs <- NA #these are the list of ALL vars for a given depth
			running.lens <- 0 #records the length of each var2.list
			for (j in seq_len(nrow(inner.mat))) { #loops through each collection of parents and finds child selections
				select2 <- resample_select(
					resp.name, inner.mat[j,], myframe, type, r, f, model, fun.args, pred.args)
				var2.list.all <- names(select2)
				var2.c.list <- select2
				p <- length(select2)
				worst_case <- wc_updater_finder(
					r, p, desir_prob, save.table = save.table, table.path = table.path)
				# worst_case <- find_worstcases_chen86(r, p, desir_prob)
				r_prime <- worst_case$r_prime
				var2.list <- var2.list.all[which(select2 >= r_prime )]
				c2.list <- select2[which(select2 >= r_prime )]

				running.vars <- c(running.vars, var2.list)
				running.cs <- c(running.cs, c2.list)
				running.lens <- c(running.lens, length(var2.list))
			}
			running.vars <- running.vars[-1]
			running.cs <- running.cs[-1]
			running.lens <- running.lens[-1]
			inner.mat2 <- matrix(NA, nrow = length(running.vars), ncol = i)
			inner.mat.c2 <- matrix(NA, nrow = length(running.cs), ncol = i)
			ind <- cumsum(running.lens) + 1
			inner.mat2[c(1, ind[-length(ind)]), 1:(i-1)] <- inner.mat
			inner.mat2[, i] <- running.vars
			inner.mat <- naRemover(inner.mat2)
			inner.mat.c2[c(1, ind[-length(ind)]), 1:(i-1)] <- inner.mat.c
			inner.mat.c2[, i] <- running.cs
			inner.mat.c <- inner.mat.c2
			if(condense){
				inner.mat3 <- t(apply(inner.mat,1,sort))
				inner.mat.c <- inner.mat.c[!duplicated(inner.mat3),,drop=FALSE]
				inner.mat <- inner.mat[!duplicated(inner.mat3),,drop=FALSE]
			}
			if(trace) print(inner.mat)
		}
	}
	return(inner.mat)
}
