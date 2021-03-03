# imports: none!
# suggests: whatever function/package you use to form models

######################################################################################################
# Read Tree Label, likely only has use for build.tree function
######################################################################################################
#' Tree Processing Function
#'
#' This is an internal function used for building MPS graph.
#' @param var.mat.na matrix of MPS data
#' @keywords matrix
#' @export
read.tree.label<-function(var.mat.na){
  # The interworkings of this method are actually somewhat complex bc use of RLE and the
  #   existance of so many NA. Rest assured I'll eventually get around to documenting
  #   all of this better, but for the time being, I won't add anything
  
  ##################################### data input and read ##########################################
  if(is.null(dim(var.mat.na))) {
    endPoint<-length(var.mat.na)
  } else {
    endPoint<-dim(var.mat.na)[2]
  }
  if(is.null(dim(var.mat.na))) {
    txt<-(var.mat.na[endPoint])
  } else {
    txt<-(var.mat.na[,endPoint])
  }
  if(is.null(dim(var.mat.na))) {
    v1<-ifelse(is.na(var.mat.na[endPoint-1]),FALSE,TRUE)
  } else {
    v1<-ifelse(is.na(var.mat.na[,endPoint-1]),FALSE,TRUE)
  }
  if(dim(var.mat.na)[2]>1){
    c<-0
    same<-0
    for(i in 1:length(v1)){
      if(v1[i]){
        c<-c+1
      }
      same[i]<-c
    }
    rle1<-rle(same)
    txt2<-rep(NA,length(rle(same)$values))
    c<-1
    for(i in 1:length(txt2)){
      txt2[i]<-paste("(",paste(txt[c:(c+rle(same)$lengths[i]-1)],collapse=","),")",sep="")
      c<-c+rle(same)$lengths[i]
    }
    if(endPoint>2){
      for(h in (endPoint-2):1){
        if(is.null(dim(var.mat.na))) v2<-ifelse(is.na(var.mat.na[h]),FALSE,TRUE) else v2<-ifelse(is.na(var.mat.na[,h]),FALSE,TRUE)
        c<-0
        same2<-0
        for(i in 1:length(v2)){
          if(v2[i]){
            c<-c+1
          }
          same2[i]<-c
        }
        rle2<-rle(same2)
        
        txt3<-rep(NA,length(rle2$values))
        rle.summer<-(rle1$lengths[1])
        if(length(rle1$lengths)>1){
          for(i in 2:length(rle1$lengths)){
            rle.summer[i]<-rle.summer[i-1]+(rle1$lengths[i])
          }
        }
        
        adjuster<-0
        c<-0
        for(i in 1:length(txt3)){
          txt3[i]<-paste("(",paste(txt2[(c+1):((which((rle.summer-adjuster)==rle2$lengths[i])))],collapse=","),")",sep="")
          c<-(which((rle.summer-adjuster)==rle2$lengths[i]))
          adjuster<-adjuster+rle2$lengths[i]
        }
        txt2<-txt3
        rle1<-rle2
      }
    }else{
      txt3<-txt2
    }
  }else{
    txt3<-paste(as.vector(var.mat.na))
  }
  final<-paste("(",paste(txt3,collapse=","), ");",sep="") # no-var at top
  return(final)
}


######################################################################################################
######################################################################################################
###################################  Build Tree Graphics Function  ###################################
######################################################################################################
######################################################################################################
#' Tree Graphic Function
#'
#' Builds MPS graphic in graphics window
#' @param var.mat.na matrix of MPS data
#' @param cex text size
#' @param alt integer that vertically staggers terminal nodes (max of 3)
#' @param shape1 either \code{'trad'} or \code{'rad'} for traditional trees and radial trees
#' @param shape2 either \code{'right'} or \code{'tri'} for rectangular and triangular;
#'  only matters for \code{shape1='trad'}
#' @param col node color
#' @keywords tree
#' @examples mat_input<-c('Var1','Var1','Var1','Var1','Var2','Var2',
#'    'Var3','Var3','Var4','Var5','Var4','Var5')
#' mat1<-matrix(mat_input,nrow=4)
#' build.tree(mat1) #base example
#' build.tree(mat1,alt=2) #different alt
#' build.tree(mat1,shape1="rad") #radial tree
#' build.tree(mat1,shape2="tri") #triangular tree (only for shape1)
#' build.tree(mat1,col="red") #node color
#' @export
build.tree<-function(var.mat.na,cex=.9, alt=1, shape1="trad", shape2="right", lwd=1, col){
  ####################################### input info and stops #######################################
  if(!is.numeric(alt)) alt<-1
  if(shape1=="rad") shape2<-NULL
  if(shape1!="trad"&&shape1!="rad") stop("invalid shape")
  if(shape2!="right"&&!is.null(shape2)) shape2<-"tri" #only matters for shape1="right"
  if(missing(col)) col<-"lightgoldenrod"
  
  ####################################### begin read and setup #######################################
  var.mat.na<-naAdder(var.mat.na)
  final<-read.tree.label(var.mat.na)
  par(mar=c(0,0,0,0))
  #plot(NA,NA,xlim=c(-1.1,1.1),ylim=c(-1.1,1.1), axes = TRUE)
  final<-gsub("_","FILLER",final)
  final<-gsub("\\.","FILLERPERIOD",final)
  str<-unlist(strsplit(final, "(?<=\\pP)(?=\\PP)|(?<=\\PP)(?=\\pP)", perl=T))
  maxDepth<-dim(var.mat.na)[2]
  x.depth<-rep(NA,length(na.omit(as.vector(var.mat.na))))
  current<-0
  str<-str[-length(str)]
  shift.right<-1
  shift<-1
  
  ################################## Node-information matrix building ##################################
  for (i in 1:length(str)){
    shift.up.down<-0
    temp.str<-unlist(strsplit(str[i],"(?<=\\pP)(?=\\pP)",perl=T))
    current.var<-sum(!is.na(x.depth))+1
    if(i==1){
      shift.up.down<-1:length(unlist(strsplit(str[i],"(?<=\\pP)(?=\\pP)",perl=T)))
      x.depth[1:length(shift.up.down)]<-shift.up.down
      shift[1:length(shift.up.down)]<-1
      current<-maxDepth+1
    }else{
      if((","%in%temp.str)&(length(temp.str)>1)){
        shift.right<-shift.right+1
        shift.up.down<-length(which(!temp.str==","))/2
        current<-current-shift.up.down
        x.depth[current.var:(current.var+(shift.up.down)-1)]<-current:(maxDepth)
        shift[current.var:(current.var+(shift.up.down)-1)]<-shift.right
        current<-(maxDepth)
      }else{
        if(!","%in%temp.str){
          x.depth[current.var]<-(maxDepth+1)
          shift[current.var]<-shift.right
          current<-(maxDepth+1)
        }
      }
    }
  }
  x.depth.true<-na.omit(x.depth)
  tot.mat<-(matrix(ncol=6,nrow=length(x.depth.true)))
  tot.mat[,1]<-c(1:length(x.depth.true))
  tot.mat[,2]<-as.numeric(x.depth.true)
  tot.mat[,3]<-as.numeric(shift)
  tot.mat[,6] #pointer to father variable
  
  ################################## Terminal node coordinates ##################################
  if(shape1=="rad"){
    r<-seq(from=0,to=1,len=maxDepth+2+alt-1)
    r<-r[-1]
    theta<-seq(from=0,to=2*pi,len=length(which(x.depth.true==maxDepth+1))+1)
    theta<-theta[-length(theta)]
  }
  if(shape1=="trad"){
    xTerm<-seq(-.95,.95,len=length(which(x.depth.true==maxDepth+1)))
    if(dim(var.mat.na)[1]==2&&dim(var.mat.na)[2]==1) xTerm<-seq(-.4,.4,len=length(which(x.depth.true==maxDepth+1)))
    if(length(which(x.depth.true==maxDepth+1))==1) xTerm<-0
    y.heights<-seq(to=-1,from=1,len=maxDepth+2)
    y.heights<-y.heights[-length(y.heights)]
    if(alt==2){
      y.heights<-c(y.heights,
                   y.heights[length(y.heights)]-strheight("1","user")*2*cex)
    }
    if(alt==3){
      y.heights<-c(y.heights,
                   y.heights[length(y.heights)]-strheight("1","user")*2*cex,
                   y.heights[length(y.heights)]-strheight("1","user")*4*cex)
    }
  }
  
  ################################## Internal node coordinates ##################################
  for(i in unique(rev(x.depth.true))){
    if(i==maxDepth+1){
      if(shape1=="trad"){
        tot.mat[which(tot.mat[,2]==maxDepth+1),4]<-xTerm
        tot.mat[which(tot.mat[,2]==maxDepth+1),5]<-y.heights[maxDepth+1]
      }
      if(shape1=="rad"){
        tot.mat[which(tot.mat[,2]==maxDepth+1),4]<-theta
        wch <- which(tot.mat[,2]==maxDepth+1)
        tot.mat[wch[seq_along(wch) %% alt == 0],5]<-r[i]
        tot.mat[wch[seq_along(wch) %% alt == 1],5]<-r[i + 1]
        tot.mat[wch[seq_along(wch) %% alt == 2],5]<-r[i + 2]
      }
      for(s in unique(tot.mat[which(tot.mat[,2]==i),3])){
        cond1<-which((tot.mat[,2]==(i-1))&(tot.mat[,3]==s))
        cond2<-which((tot.mat[,3]-s)==min(tot.mat[cond1,3]-s))
        final.cond<-cond2[cond2%in%cond1]
        pointer<-tot.mat[final.cond,1]
        tot.mat[which((tot.mat[,2]==i)&(tot.mat[,3]==s)),6]<-pointer
      }
    }else{
      for(s in unique(tot.mat[which(tot.mat[,2]==i),3])){
        if(i!=1){
          cond1<-which((tot.mat[,2]==(i-1))&(tot.mat[,3]<=s))
          cond1<-max(cond1)
          if(0%in%(tot.mat[cond1,3]-s)){
            cond2<-which((tot.mat[,3]-s)==0)
          }else{
            cond2<-which((tot.mat[,3]-s)==min(tot.mat[cond1,3]-s))
          }
          final.cond<-cond2[cond2%in%cond1]
          pointer<-tot.mat[final.cond,1]
          tot.mat[which((tot.mat[,2]==i)&(tot.mat[,3]==s)),6]<-pointer
        }
        if(i==1){
          tot.mat[which(tot.mat[,2]==1),6]<-0
        }
        
        if(shape1=="trad"){
          c.row<-tot.mat[which((tot.mat[,2]==i)&(tot.mat[,3]==s)),] #current row being editted
          x.val<-mean(tot.mat[which(tot.mat[,6]==c.row[1]),4])
          y.val<-y.heights[i]
          tot.mat[which((tot.mat[,2]==i)&(tot.mat[,3]==s)),4]<-x.val
          tot.mat[which((tot.mat[,2]==i)&(tot.mat[,3]==s)),5]<-y.val
        }
        if(shape1=="rad"){
          c.row<-tot.mat[which((tot.mat[,2]==i)&(tot.mat[,3]==s)),] #current row being editted
          theta.val<-mean(tot.mat[which(tot.mat[,6]==c.row[1]),4])
          r.val<-r[i]
          tot.mat[which((tot.mat[,2]==i)&(tot.mat[,3]==s)),4]<-theta.val
          tot.mat[which((tot.mat[,2]==i)&(tot.mat[,3]==s)),5]<-r.val
        }
      }
    }
  }
  
  ############################# Terminal node height alternator ##############################
  if(alt>1&&shape1=="trad"){
    indeces.alt<-which(tot.mat[,2]==maxDepth+1)
    if(alt==2){
      if(length(indeces.alt)%%2==0){
        tot.mat[indeces.alt,5]<-c(y.heights[maxDepth+1],y.heights[maxDepth+2])
      }else{
        tot.mat[indeces.alt,5]<-c(rep(c(y.heights[maxDepth+1],
                                        y.heights[maxDepth+2]),floor(length(indeces.alt)/2)),
                                  y.heights[maxDepth+1])
      }
    }
    if(alt==3){
      if(length(indeces.alt)%%3==0){
        tot.mat[indeces.alt,5]<-c(y.heights[maxDepth+1],y.heights[maxDepth+2],y.heights[maxDepth+3])
      }
      if(length(indeces.alt)%%3==1){
        tot.mat[indeces.alt,5]<-c(rep(c(y.heights[maxDepth+1],
                                        y.heights[maxDepth+2],
                                        y.heights[maxDepth+3]),floor(length(indeces.alt)/3)),
                                  y.heights[maxDepth+1])
      }
      if(length(indeces.alt)%%3==2){
        tot.mat[indeces.alt,5]<-c(rep(c(y.heights[maxDepth+1],
                                        y.heights[maxDepth+2],
                                        y.heights[maxDepth+3]),floor(length(indeces.alt)/3)),
                                  y.heights[maxDepth+1],y.heights[maxDepth+2])
      }
    }
  }
  
  ######################################## Plotting lines #########################################
  plot(NA,NA,xlim=c(-1.1,1.1),ylim=c(-1.1,1.1), axes = FALSE)
  if(shape1=="trad"){
    for(k in dim(tot.mat)[1]:1){
      if(tot.mat[k,6]>1){
        if(shape2=="right"){
          segments(tot.mat[k,4],tot.mat[k,5],
                   tot.mat[k,4],
                   tot.mat[which(tot.mat[,1]==tot.mat[k,6]),5],col=if(k<maxDepth+2)"black"else"black", lwd=lwd) #vertical line
          segments(tot.mat[k,4],tot.mat[which(tot.mat[,1]==tot.mat[k,6]),5],
                   tot.mat[which(tot.mat[,1]==tot.mat[k,6]),4],
                   tot.mat[which(tot.mat[,1]==tot.mat[k,6]),5],col=if(k<maxDepth+2)"black"else"black", lwd=lwd) #horizontal line
        }else{
          segments(tot.mat[k,4],tot.mat[k,5],
                   tot.mat[which(tot.mat[,1]==tot.mat[k,6]),4],
                   tot.mat[which(tot.mat[,1]==tot.mat[k,6]),5],col=if(k<maxDepth+2)"black"else"black", lwd=lwd)
        }
      }
    }
  }
  if(shape1=="rad"){
    for(k in dim(tot.mat)[1]:1){
      if(tot.mat[k,6]>1){
        new.theta<-tot.mat[which(tot.mat[,1]==tot.mat[k,6]),4]
        new.r<-tot.mat[which(tot.mat[,1]==tot.mat[k,6]),5]
        segments(tot.mat[k,5]*cos(tot.mat[k,4]),tot.mat[k,5]*sin(tot.mat[k,4]),
                 new.r*cos(tot.mat[k,4]),new.r*sin(tot.mat[k,4]),col=if(k<maxDepth+2)"black"else"black", lwd=lwd) #radial line
        lines(new.r*cos(seq(tot.mat[k,4],new.theta,len=100)),
              new.r*sin(seq(tot.mat[k,4],new.theta,len=100)),col=if(k<maxDepth+2)"black"else"black", lwd=lwd)
      }
    }
  }
  
  #################################### Plotting words and boxes #####################################
  mywords<-c(na.omit(as.vector(t(var.mat.na))))
  mywords<-gsub("FILLER","_",mywords)
  mywords<-gsub("FILLERPERIOD","\\.",mywords)
  
  if(shape1=="trad"){
    rect(tot.mat[-1,4]-strwidth(mywords,"user")*cex/1.9,
         tot.mat[-1,5]-strheight(mywords,"user")*cex/1.1,
         tot.mat[-1,4]+strwidth(mywords,"user")*cex/1.9,
         tot.mat[-1,5]+strheight(mywords,"user")*cex/1.1,col=col, border = col, lwd=lwd)
    lx<-tot.mat[-1,4]
    ly<-tot.mat[-1,5]
  }
  if(shape1=="rad"){
    rect(tot.mat[-1,5]*cos(tot.mat[-1,4])-strwidth(mywords,"user")*cex/1.9,
         tot.mat[-1,5]*sin(tot.mat[-1,4])-strheight(mywords,"user")*cex/1.1,
         tot.mat[-1,5]*cos(tot.mat[-1,4])+strwidth(mywords,"user")*cex/1.9,
         tot.mat[-1,5]*sin(tot.mat[-1,4])+strheight(mywords,"user")*cex/1.1,col=col, border = col, lwd=lwd)
    lx<-tot.mat[-1,5]*cos(tot.mat[-1,4])
    ly<-tot.mat[-1,5]*sin(tot.mat[-1,4])
  }
  text(lx,ly,mywords,cex=cex,col="black")
}

######################################################################################################
# naAdder, adds NAs, likely only has use or applicability for full.select method
######################################################################################################
#' Add NAs to MPS Graph Matrix
#'
#' Internal function that adds \code{NA}s to matrix that holds MPS graph data
#' @param myMat matrix of MPS data
#' @keywords matrix
#' @examples mat_input<-c('Var1','Var1','Var1','Var1','Var2','Var2',
#'    'Var3','Var3','Var4','Var5','Var4','Var5')
#' mat1<-matrix(mat_input,nrow=4)
#' naAdder(mat1)
#' @export
naAdder<-function(myMat){
  x<-dim(myMat)[1]
  y<-dim(myMat)[2]
  sticker<-rep(FALSE,x)
  if(y>1){
    for(i in 1:(y-1)){
      myRle<-rle(myMat[,i])
      toNA<-cumsum(myRle$lengths)+1
      toNA<-c(1,toNA[-length(toNA)])
      indNA<-rep(FALSE,length(myMat[,i]))
      indNA[toNA]<-TRUE
      sticker<-sticker|(indNA)
      myMat[!sticker,i]<-NA
    }
  }
  return(myMat)
}

######################################################################################################
# naRemover, fills in NAs, likely only has use or applicability for full.select method
######################################################################################################
#' Remove NAs to MPS Graph Matrix
#'
#' Internal function that removes \code{NA}s to matrix that holds MPS graph data
#' @param myMat matrix of MPS data
#' @keywords matrix
#' @examples mat_input<-c('Var1',NA,NA,NA,'Var2',NA,
#'    'Var3',NA,'Var4','Var5','Var4','Var5')
#' mat1<-matrix(mat_input,nrow=4)
#' naRemover(mat1)
#' @export
naRemover<-function(myMat){
  y<-dim(myMat)[2]
  for(i in 1:(y-1)){
    tv<-ifelse(is.na(myMat[,i]),"REP",myMat[,i])
    myRle<-rle(tv)
    toRep<-na.omit(ifelse(myRle$values!="REP"&c(myRle$values[-1],"HOLDER")!="REP",myRle$lengths,
                          ifelse(myRle$values!="REP"&c(myRle$values[-1],"HOLDER")=="REP",
                                 c(myRle$lengths[-1],1)+myRle$lengths,NA)))
    myMat[,i]<-rep(myRle$values[myRle$values!="REP"],toRep)
  }
  return(myMat)
}

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
#' boot.select.gen(resp.name='y',c.vars=NULL,myframe=mydata,B=100,model='lm')
#' @export
boot.select.gen<-function(resp.name,c.vars,myframe,type,r,f,model,fun.args,pred.args){#c.vars are chosen, a.vars are available
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
  
  B <- r*ncol(myframe)*1.2
  chosen1<-rep(NA,B)
  a.vars<-colnames(myframe)[which(!colnames(myframe)%in%c(c.vars,resp.name))]
  for(i in 1:B){
    samp<-sample(1:dim(myframe)[1],floor(sqrt(dim(myframe)[1])),replace=FALSE) ### CHECK this is bootstrap at a rate ALSO below
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

f<-function(y,y.pred){return(mean((y-y.pred)^2))} #mean squared error


######################################################################################################
#
######################################################################################################
#' Find set likely to contain the most probable category
#'
#' Function that that finds the set likely to contain the most probable multinomial category
#' @param r maximum cell count
#' @param k number of cells
#' @param P desired inclusion probability
#' @examples find_worstcases_chen86(100. 5, 0.80)
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

######################################################################################################
#
######################################################################################################
#' Update table
#'
#' Function that that updates a table containing info generated by `find_worstcases_chen86()`
#' @param r maximum cell count
#' @param k number of cells
#' @param P desired inclusion probability
#' @examples find_worstcases_chen86(100. 5, 0.80)
wc_updater_finder <- function(r, k, P) {
  mypath <- "~/mps/speed_test/worst_case_table.Rds"
  # mypath <- "/ihome/lmentch/nkissel/mps/speed_test/worst_case_table.Rds"
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
    saveRDS(worst_cases, file = "~/mps/speed_test/worst_case_table.Rds", version = 2)
    return(list(r_prime = t$r_prime, P = P))
  }
}

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
#' @keywords matrix
#' @examples set.seed(200)
#' library(MASS)
#' n<-1000
#' p<-10
#' x<-mvrnorm(n,rep(0,p),diag(p))
#' signal<-rowSums(x[,1:3])
#' noise<-rnorm(n,0,1/4)
#' y<-signal+noise
#' mydata<-data.frame(x,y)
#' mps1<-full.select.gen(myframe=mydata,resp.name='y',depth=3,
#'    B=100,model='lm',condense=FALSE)
#' mps1
#' build.tree(mps1) #graph
#'
#' #merged paths
#' mps2<-full.select.gen(myframe=mydata,resp.name='y',depth=3,
#'    r=100,model='lm',condense=TRUE)
#' mps2
#' build.tree(mps2) #graph
#' @export
full.select.gen<-function(myframe,resp.name,depth,r,Pstar,f,type,model,fun.args,pred.args,condense=FALSE){
  ####################################### input info and stops #######################################
  if(missing(depth)) depth<-3
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
  if(missing(fun.args)) fun.args<-list()
  if(missing(pred.args)) pred.args<-list()
  
  ############################################ Data Read #############################################
  var.list<-""
  vars<-colnames(myframe)[which(colnames(myframe)!=resp.name)]
  
  ######################################### Depth 1 select ###########################################
  select1<-boot.select.gen(resp.name,NULL,myframe,type,r,f,model,fun.args,pred.args)
  var1.list.all<-names(select1)
  var1.c.list<-select1
  
  p <- length(select1)
  desir_prob <- Pstar
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
        select2 <- boot.select.gen(resp.name, inner.mat[j,],myframe,type,r,f,model,fun.args,pred.args)
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
  return(inner.mat)
}

