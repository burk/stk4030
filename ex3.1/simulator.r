library(MASS)

exact.b <- c(3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0);
exact.sd <- 5;
exact.p <- 15;

exact.generateX <- function(N) {
  mu <- rep(0, exact.p);
  sigma <- array(0.8, c(exact.p, exact.p)) + diag(0.2, exact.p);
  return(mvrnorm(n=N, mu=mu, sigma));
}

exact.generateY <- function(X) {
  N <- dim(X)[1];
  
  Y <- cbind(1, X) %*% exact.b + rnorm(N, 0, exact.sd);
  return(Y)
}

beta.res <- c();
pred.res <- c();

for (N in c(20, 100, 1000)) {
  X.train <- exact.generateX(N);
  X.test  <- exact.generateX(1000);

  nrep <- 1000
  beta.errs <- matrix(, nrow=length(exact.b), ncol=nrep);
  pred.errs <- matrix(, nrow=length(X.test), ncol=nrep);
  for (rep in 1:nrep) {
    Y.train <- exact.generateY(X.train);
    Y.test  <- exact.generateY(X.test);

    lmod <- lm(Y~., data.frame(Y=Y.train, X=X.train));
    beta.errs[, rep] <- lmod$coef - exact.b;
    
    Y.pred <- predict(lmod, data.frame(X=X.test));
    pred.errs[, rep] <- Y.pred - Y.test;
  }
  beta.bias <- apply(beta.errs, 1, mean);
  beta.var  <- apply(beta.errs, 1, var);
  beta.mse  <- apply(beta.errs^2, 1, mean);
  
  tmp <- cbind(beta.bias, beta.var, beta.mse);
  colnames(tmp) <- c('beta bias', 'beta var', 'beta MSE');
  beta.res <- cbind(beta.res, tmp);
  
  pred.bias <- mean(as.vector(pred.errs));
  pred.var  <- var(as.vector(pred.errs));
  pred.mse  <- mean(as.vector(pred.errs)^2);
  
  tmp <- cbind(pred.bias, pred.var, pred.mse);
  colnames(tmp) <- c('pred bias', 'pred var', 'pred MSE');
  pred.res <- cbind(pred.res, tmp);
}
