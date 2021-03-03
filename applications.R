################################################################################
################################### Diabetes ###################################
################################################################################
library(lars)
library(ggplot2)
source("../R/myPackage.R")

data(diabetes)
xnm <- colnames(diabetes$x)
mydata <- cbind(diabetes$x, diabetes$y)
colnames(mydata) <- c(xnm, "y")
mydata <- data.frame(mydata)
tosq <- as.matrix(mydata[,-c(11)])
newd_2 <- poly(tosq, degree = 2, raw = TRUE)
nm <- strsplit(colnames(newd_2), "[.]")
for(i in seq_along(nm)) {
	which1 <- which(nm[[i]] == "1")
	if(length(which1) > 0) {
		colnames(newd_2)[i] <- paste0(xnm[which1], collapse = "_")
	} else {
		which2 <- which(nm[[i]] == "2")
		colnames(newd_2)[i] <- paste0(xnm[which2], "_", xnm[which2])
	}
	# newd_2[,i] <- (newd_2[,i] - mean(newd_2[,i]))/sd(newd_2[,i])
}
newd_2 <- newd_2[, which(colnames(newd_2) != "sex_sex")]
final_d <- cbind(newd_2, data.frame(y = diabetes$y))
final_d <- as.data.frame(final_d)
final_d$y <- (final_d$y - mean(final_d$y))/sd(final_d$y)

# number of steps of fwd sel
set.seed(1)
sss <- sample(1:nrow(final_d), 300)
final_d_sub <- final_d[sss, ]
resp.name <- 'y'
nfolds <- 5
mse <- matrix(nrow = 5, ncol = ncol(final_d_sub) - 1)
inds <- sample(rep_len(1:nfolds, nrow(final_d_sub)), nrow(final_d_sub))
for(i in 1:nfolds) {
	train <- final_d_sub[inds != i, ]
	test <- final_d_sub[inds == i, ]
	lower<-lm(as.formula(paste(resp.name,"~","1",sep="")),train)
	upper<-lm(as.formula(paste(resp.name,"~",".",sep="")),data=train)
	mystep<-step(lower,scope=list(lower=lower,upper=upper),direction="forward",steps=(ncol(final_d_sub) - 1),trace=0,k=0)
	pth <- names(colSums(attr(mystep$terms, 'factors') == 1))
	mse_temp <- NULL
	for(j in 1:(ncol(final_d_sub) - 1)) {
		lm_temp <- lm(as.formula(paste(resp.name,"~",paste0(pth[1:j], collapse='+'),sep="")),data=train)
		mse_temp[j] <- mean((predict(lm_temp, test) - test$y)^2)
	}
	mse[i,] <- mse_temp
}
which.min(colMeans(mse))
mse[which.min(colMeans(mse))]

# fwd selection
lower<-lm(as.formula(paste(resp.name,"~","1",sep="")),final_d_sub)
upper<-lm(as.formula(paste(resp.name,"~",".",sep="")),data=final_d_sub)
mystep<-step(lower,scope=list(lower=lower,upper=upper),direction="forward",steps=which.min(colMeans(mse)),trace=0,k=0)
fs_err <- mean((predict(mystep, final_d[-sss, ]) - final_d[-sss, ]$y)^2)

set.seed(1)
full.sel1 <- full.select.gen(myframe = final_d_sub, resp.name = "y",
														 r = 100, depth = 7, Pstar = 0.95, condense = TRUE)
# tree plot
build.tree(full.sel1, shape1 = "rad", cex = 0.55, alt = 3)

mse_mps <- NULL
for(i in 1:nrow(full.sel1)) {
	pth <- full.sel1[i,]
	lm_temp <- lm(as.formula(paste(resp.name,"~",paste0(pth, collapse='+'),sep="")),data=train)
	mse_mps[i] <- mean((predict(lm_temp, test) - test$y)^2)
}


set.seed(1)
library(glmnet)
cvgl <- cv.glmnet(as.matrix(final_d_sub[,colnames(final_d_sub) != 'y']),
									final_d_sub[,colnames(final_d_sub) == 'y'], nfolds = 5)
las <- glmnet(as.matrix(final_d_sub[,colnames(final_d_sub) != 'y']),
							final_d_sub[,colnames(final_d_sub) == 'y'], lambda = cvgl$lambda.min)
sum(as.matrix(las$beta)!=0)
las_err <- mean((predict(las, as.matrix(final_d[-sss, colnames(final_d) != 'y'])) - final_d[-sss,]$y)^2)

sca <- 1.5
library(ggplot2)
plt <- ggplot(data.frame(mse_mps), aes(x = mse_mps, y = "MPS")) +
	geom_boxplot() +
	theme_bw() +
	labs(x = "Test error", y="") +
	geom_vline(data = data.frame(y="MPS", x = fs_err), aes(xintercept=x), color="black", lty = 3, size = 1) +
	geom_vline(data = data.frame(y="MPS", x = las_err), aes(xintercept=x), color="black", lty = 5, size = 1) +
	theme(axis.text = element_text(size = 6 * sca),
				axis.title = element_text(size = 8 * sca),
				plot.title = element_text(size = 8* sca),
				legend.text = element_text(size = 7 * sca),
				strip.text = element_text(size = 6 * sca),
				strip.background = element_rect(),
				panel.grid.minor = element_blank())
plt

################################################################################
################################ Breast Cancer #################################
################################################################################
library(parallel)
mydata <- as.data.frame(read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.data", header = F))
mydata <- mydata[,-1]
colnames(mydata)
mydata$V11 <- ifelse(mydata$V11 == 2, 0, 1)
prev <- mean(mydata$V11)
colnames(mydata) <- c("thickness", "size_unif", "shape_unif", "margin_ad", "epi_size",
											"br_nuclei", "bl_chrm", "norm_nucl", "mitoses", "class")
o_name <-c("thickness", "size_unif", "shape_unif", "margin_ad", "epi_size",
					 "br_nuclei", "bl_chrm", "norm_nucl", "mitoses")
n_name <-c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9")

library(rpart)
class_and_error <- function(y,y.pred){
	return(sum(y != ifelse(y.pred > prev, 1, 0))/length(y))
}
set.seed(999)
full.sel1 <- full.select.gen(myframe = mydata, resp.name = "class",
														 r = 200, type = "class", depth = 3, Pstar = 0.75, condense = TRUE,
														 model = "glm", fun.args = list(family = binomial()),
														 pred.args = list(type = "response"), f = class_and_error)

set.seed(999)
full.sel2<-full.select.gen(myframe = mydata, resp.name = "class",
													 r = 200, type = "class", depth = 3, Pstar = 0.75, condense = TRUE,
													 model = "rpart")

conv_x <- function(m) {
	ret <- ifelse(m == o_name[1], n_name[1],
								ifelse(m == o_name[2], n_name[2],
											 ifelse(m == o_name[3], n_name[3],
											 			 ifelse(m == o_name[4], n_name[4],
											 			 			 ifelse(m == o_name[5], n_name[5],
											 			 			 			 ifelse(m == o_name[6], n_name[6],
											 			 			 			 			 ifelse(m == o_name[7], n_name[7],
											 			 			 			 			 			 ifelse(m == o_name[8], n_name[8], n_name[9]))))))))
	return(ret)
}

# pdf("breastcancer_paths.pdf", width = 12*1.5, height = 4*1.5)
sep_s <- which(!is.na(naAdder(full.sel1)[,1]))
layout(matrix(c(1,2,3,1,4,5,6,7,8), 3, 3, byrow = TRUE),
			 widths=c(1/49,24/49,24/49), heights=c(1/3,1/3,1/3))
par(mar = rep(0, 4))
plot.new()
rect(0,0,1,1, col = "lightgrey", border = 'lightgrey')
text(0.5,0.5,labels = "Logistic Regression",srt = 90,cex = 2, adj =c(0.5,0.5))
build.tree(conv_x(full.sel1)[sep_s[1]:(sep_s[2] - 1),], alt = 1, cex = 1.6, lwd = 2)
build.tree(conv_x(full.sel1)[sep_s[2]:(sep_s[3] - 1),], alt = 1, cex = 1.6, lwd = 2)
build.tree(conv_x(full.sel1)[sep_s[3]:(sep_s[4] - 1),], alt = 1, cex = 1.6, lwd = 2)
build.tree(conv_x(full.sel1)[sep_s[4]:nrow(full.sel1),], alt = 1, cex = 1.6, lwd = 2)
plot.new()
rect(0,0,1,1, col = "lightgrey", border = 'lightgrey')
text(0.5,0.5,labels = "Trees",srt = 90,cex = 2, adj =c(0.5,0.5))
build.tree(conv_x(full.sel2), alt = 1, cex = 1.6, lwd = 2)
# dev.off()

