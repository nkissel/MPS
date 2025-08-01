######################################################################################################
# General selection proportion generator, defined for use in full.select
######################################################################################################
#' Perform a Step of FSS
#'
#' Internal function that performs 1 step of forward stability selection. Returns a table of selection counts.
#' @param resp.name matrix of MPS data
#' @param c.vars matrix of MPS data
#' @param myframe data.frame holding explanatory and response data
#' @param type 'reg' or 'class' for regression and classification. Not needed if `f` is specified.
#' @param r maximum cell count
#' @param f loss function
#' @param model name of modeling method; must be name of the function called in R
#' @param fun.args list of arguments that are passed to modeling function
#' @param pred.args list of arguments that are passed to model predict function
#' @keywords matrix
#' @examples library(MASS)
#' n<-100
#' p<-10
#' x<-mvrnorm(n,rep(0,p),diag(p))
#' signal<-rowSums(x[,1:3])
#' noise<-rnorm(n,0,1/4)
#' y<-signal+noise
#' mydata<-data.frame(x,y)
#' resample_select(resp.name='y',c.vars=NULL,myframe=mydata,r=100,model='lm')
#' @export
resample_select<-function(resp.name,c.vars,myframe,type,r,f,model,fun.args,pred.args,subsampr){#c.vars are chosen, a.vars are available
	if(missing(resp.name)) stop("Need response variable name")
	if(missing(c.vars)) c.vars<-NULL
	if(!is.character(c.vars)&&!is.null(c.vars)) stop("c.vars must be of type character or NULL")
	if(missing(myframe)) stop("Need data.frame")
	if(missing(type)) type<-"reg"
	if(type!="reg"&&type!="class") stop("Need to specify type as \"class\" or \"reg\"")
	if(missing(r)) r<-100
	if(missing(f)){
		f<-function(y,y.pred){return(mean((y-y.pred)^2))} #mean squared error
		if(type=="class"){
			f<-function(y,y.pred){return(sum(y!=y.pred)/length(y))} #misclassification rate
		}
	}
	if(!is.function(f)) stop("f must be a function")
	if(missing(model)) model<-"lm"
	my.try<-try(match.fun(model),silent=TRUE)
	if("try-error" %in% class(my.try)) stop(paste("You do not have package with \"",
																								model,"\" as a function installed.",sep=""))
	if(missing(fun.args)) fun.args<-list()
	if(missing(pred.args)) pred.args<-list()
	if(missing(subsampr)) subsampr<-function(x) sqrt(x)

	B <- r*ncol(myframe)*1.2
	chosen1<-rep(NA,B)
	a.vars<-colnames(myframe)[which(!colnames(myframe)%in%c(c.vars,resp.name))]
	for(i in 1:B){
		samp<-sample(1:dim(myframe)[1],floor(subsampr(dim(myframe)[1])),replace=FALSE) ### CHECK this is bootstrap at a rate ALSO below
		train<-myframe[samp,]
		acc.temp<-rep(1,length(a.vars))
		for(j in 1:length(a.vars)){
			fun.new.args<-c(fun.args,
											list(formula=as.formula(paste(resp.name,paste(c(c.vars,a.vars[j]),collapse=" + "),sep=" ~ ")),data=train))
			mod<-do.call(model,fun.new.args)
			pred.new.args<-c(pred.args,list(object=mod,newdata=train))
			pred.temp<-do.call("predict",pred.new.args)
			acc.temp[j]<-f(pred.temp,train[,which(colnames(train)==resp.name)])
		}

		if(length(which(acc.temp==min(acc.temp)))>1){
			chosen1[i]<-a.vars[sample(which(acc.temp==min(acc.temp)),1)]
		}
		if(length(which(acc.temp==min(acc.temp)))==1){
			chosen1[i]<-a.vars[which.min(acc.temp)]
		}
		if(max(table(chosen1)) >= r) {
			break
		}
	}
	rt <- table(c(chosen1, a.vars))[order(table(c(chosen1, a.vars)),decreasing=TRUE)]
	rt <- rt - 1
	return(rt)
}
