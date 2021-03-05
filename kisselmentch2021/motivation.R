library(ggplot2)
library(MASS)
library(glmnet)

set.seed(1)
n <- 500
p <- 18
SNR <- 20

# covariance matrix
mysigma <- diag(p)
mysigma[1:3,1:3] <- mysigma[4:6,4:6] <- mysigma[7:9, 7:9] <- 0.90
mysigma[diag(p) == 1] <- 1

x <- mvrnorm(n, rep(0, p), mysigma) # covariates

cf <- c(3, 2, 1) #coefficients

# generate response
signal <- cf[1]*rowSums(x[,1:3]) + cf[2]*rowSums(x[,4:6]) + cf[3]*rowSums(x[,7:9])
eps <- rnorm(n, 0, sqrt(var(signal) / SNR))
y <- signal + eps
train <- data.frame(x,y)

ls1 <- glmnet(x, y)
myl <- c(ls1$lambda[1] + 10, ls1$lambda, 0)

# SS lasso
B <- 100
ss_mat <- matrix(nrow = B, ncol = p)
ss_mat_ls <- matrix(nrow = B, ncol = p)
for(i in 1:B) {
  tr <- train[sample(1:n, n/2, replace = F),]
  lower <- lm(y ~ 1, tr)
  upper <- lm(y ~ ., tr)
  mystep <- step(lower, scope = list(lower = lower, upper = upper),
                 direction="forward", steps = p,trace = 0,k = 0)
  ls_temp <- glmnet(as.matrix(tr[, -(p+1)]), as.vector(tr[, (p+1)]), lambda = myl)
  b <- t(as.matrix(ls_temp$beta))
  first_nonzero <- function(x){min(which(x != 0))}
  loc_add <- apply(b, 2, first_nonzero)
  ss_mat_ls[i, ] <- loc_add
  ss_mat[i, ] <- names(coef(mystep)[-1])
}

var_names <- names(train[,1:p])
colnames(ss_mat_ls) <- var_names
var_counts <- matrix(0, ncol = p, nrow = p)
for(i in 1:p) {
  tb <- table(ss_mat[,1:i])
  var_counts[i,match(names(tb), var_names)] <- tb
}
colnames(var_counts) <- var_names
rownames(var_counts) <- 1:p

var_counts_ls <- matrix(0, nrow = length(myl), ncol = p)
for(i in 1:nrow(var_counts_ls)) {
  has_been_picked <- function(x){length(which(x <= i))}
  var_counts_ls[i, ] <- apply(ss_mat_ls, 2, has_been_picked)
}
colnames(var_counts_ls) <- var_names

var_prop <- var_counts / B
var_prop_ls <- var_counts_ls / B
types <- as.factor(c(rep("X1-X3", p*3), rep("X4-X6", p*3), rep("X7-X9", p*3),
                     rep("X10-X18", p*p - p*9) ))
types <- factor(types, levels = c("X1-X3", "X4-X6", "X7-X9", "X10-X18"))
types_ls <- as.factor(c(rep("X1-X3", length(myl)*3), rep("X4-X6", length(myl)*3),
                        rep("X7-X9", length(myl)*3),
                        rep("X10-X18", length(myl)*p - length(myl)*9) ))
types_ls <- factor(types_ls, levels = c("X1-X3", "X4-X6", "X7-X9", "X10-X18"))
to_plot <- data.frame(props = as.vector(var_prop), reg = rep(1:p, p),
                      var = rep(var_names, each = p), type = types)
to_plot_ls <- data.frame(props = as.vector(var_prop_ls), reg = 1 / (rep(myl, p) + 1),
                         var = rep(var_names, each = length(myl)), type = types_ls)


p0 <- ggplot(to_plot_ls,aes(x = reg, y = props,group = var, col = type)) + 
  geom_line() + labs(col = "Variable group", y = "Selection proportion",
                     x = "Regularization") + theme_bw()

# generate test data
x <- mvrnorm(10000, rep(0, p), mysigma)
signal <- cf[1]*rowSums(x[,1:3]) + cf[2]*rowSums(x[,4:6]) + cf[3]*rowSums(x[,7:9])
eps <- rnorm(10000, 0, sqrt(var(signal) / SNR))
y <- signal + eps
test <- data.frame(x,y)

# get all 3 group model errors
errs_3grp <- rep(0, 27)
ctr <- 1
for(i in 1:3) {
  for(j in 4:6) {
    for(k in 7:9) {
      train_temp <- data.frame(train[,i], train[,j],train[,k], train$y)
      colnames(train_temp) <- c(paste0("X", i), paste0("X", j), paste0("X", k), "y")
      
      test_temp <- data.frame(test[,i], test[,j],test[,k], test$y)
      colnames(test_temp) <- c(paste0("X", i), paste0("X", j), paste0("X", k), "y")
      
      lm_temp <- lm(y ~ ., train_temp)
      errs_3grp[ctr] <- mean((test_temp$y - predict(lm_temp, test_temp))^2)
      ctr <- ctr + 1
    }
  }
}

# get 3 var 2 group model errors (grp1 & grp2)
ctr <- 1
opts <- combn(1:6, 3)
errs_2grp <- rep(0, ncol(opts))
err123 <- NA
for(i in 1:ncol(opts)) {
  train_temp <- data.frame(train[,opts[,i][1]], train[,opts[,i][2]],train[,opts[,i][3]], train$y)
  colnames(train_temp) <- c(paste0("X", opts[,i][1]), paste0("X", opts[,i][2]),
                            paste0("X", opts[,i][3]), "y")
  
  test_temp <- data.frame(test[,opts[,i][1]], test[,opts[,i][2]],test[,opts[,i][3]], test$y)
  colnames(test_temp) <- c(paste0("X", opts[,i][1]), paste0("X", opts[,i][2]),
                           paste0("X", opts[,i][3]), "y")
  
  lm_temp <- lm(y ~ ., train_temp)
  errs_2grp[ctr] <- mean((test_temp$y - predict(lm_temp, test_temp))^2)
  ctr <- ctr + 1
}

p1 <- ggplot(to_plot, aes(x = reg, y = props, group = var, col = type)) + 
  geom_line() + 
  labs(col = "Variable group", y = "Selection proportion", x = "Variables selected") + 
  theme_bw()
p2 <- ggplot(data.frame(errs_3grp), aes(x = errs_3grp)) + 
  geom_histogram(fill = "gray", col = "black", bins = 30) + 
  geom_vline(aes(xintercept = min(errs_2grp, na.rm = T)), lty = 2) +
  labs(y = "Count", x = "Mean squared error") + 
  theme_bw()

p0 # lasso stability paths
p1 # fwd sel stability paths
p2 # 3 group errors & best 2 grp
