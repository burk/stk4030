library(MASS)
library(leaps)

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

kfold.cv <- function(k, X, Y, fn, parms) {
  N <- dim(X.train)[1];
  
  per.fold <- floor(N / k);
  
  ind <- 1:N;
  ind <- sample(ind);
  
  res <- c();
  for (i in 1:k) {
    bot <- per.fold * (i-1) + 1;
    top <- bot + per.fold - 1;
    
    tmp <- fn(X[-ind[bot:top],], Y[-ind[bot:top],], X[ind[bot:top],], Y[ind[bot:top],], parms);
    
    res <- cbind(res, tmp);
  }
  
  return(res);
}

cv.leap <- function(X, Y, testX, testY, parms) {
  res <- leaps(X, Y, nbest=1);

  errsq <- c();
  
  for (incl in parms) {
    best <- res$which[incl, ];
    lmod <- lm(Y~., data.frame(X=X[, best], Y=Y));
    pred <- predict(lmod, data.frame(X=testX[, best]));
    errsq <- rbind(errsq, as.array(pred - testY)^2);
  }
  
  return(errsq);
}

cv.ridge <- function(X, Y, testX, testY, parms) {
  errsq <- c();
  
  lam <- parms;
           
  for (l in lam) {
    lmod <- lm.ridge(Y~., data.frame(X=X, Y=Y), lambda=l);

    coef <- coef(lmod);
    pred <- cbind(1, testX) %*% coef;
    
    errsq <- cbind(errsq, as.array(pred - testY)^2);
  }
  
  return(t(errsq));
}

betadf <- data.frame();
preddf <- data.frame();

for (experiment in 1:3) {
  niter <- c(1000, 500, 200)[experiment];
  ntrain <- c(20, 100, 500)[experiment];
  
  for (iter in 1:niter) {
    cat("Iteration", iter, "\n");
    
    # Ordinary Least Squares regression
    X.train <- exact.generateX(ntrain);
    Y.train <- exact.generateY(X.train);
    
    X.test <- exact.generateX(1000);
    Y.test <- exact.generateY(X.test);
    
    lmod <- lm(Y~., data.frame(X=X.train, Y=Y.train));
    pred.ols <- predict(lmod, data.frame(X=X.test));
    pred.ols.err <- pred.ols - Y.test;
    beta.ols.err <- unname(lmod$coef - exact.b);
    betadf <- rbind(betadf, data.frame(niter=niter, ntrain=ntrain, type="ols", err=beta.ols.err, beta=0:15));
    preddf <- rbind(preddf, data.frame(niter=niter, ntrain=ntrain, type="ols", err=pred.ols.err));
    
    # Best subset regression
    parms <- 1:dim(X.train)[2];
    out <- kfold.cv(10, X.train, Y.train, cv.leap, parms);
    leap.rmse <- sqrt(apply(out, 1, mean));
    best <- which.min(leap.rmse);
    cat("Best subset size", best, "\n");
    bmod <- leaps(X.train, Y.train, nbest=1);
    include <- bmod$which[best,];
    
    lbmod <- lm(Y~., data.frame(X=X.train[, include], Y=Y.train))
    pred.best <- predict(lbmod, data.frame(X=X.test[, include]));
    pred.best.err <- pred.best - Y.test;
    beta <- rep(0, exact.p + 1);
    beta[c(TRUE, include)] <- lbmod$coef
    beta.best.err <- unname(beta - exact.b);
    betadf <- rbind(betadf, data.frame(niter=niter, ntrain=ntrain, type="best", err=beta.best.err, beta=0:15));
    preddf <- rbind(preddf, data.frame(niter=niter, ntrain=ntrain, type="best", err=pred.best.err));
    
    # Ridge regression
    parms <- c(50, 25, 20, 17.5, 15, 12.5, 10, 7.5, 5, 2.5, 1);
    out <- kfold.cv(10, X.train, Y.train, cv.ridge, parms);
    ridge.rmse <- sqrt(apply(out, 1, mean));
    lambda <- parms[which.min(ridge.rmse)];
    cat("Ridge lambda parameter", lambda, "\n")
    
    rmod <- lm.ridge(Y~., data.frame(X=X.train, Y=Y.train), lambda=lambda);
    coef <- coef(rmod);
    pred.ridge <- cbind(1, X.test) %*% coef;
    pred.ridge.err <- pred.ridge - Y.test;
    beta.ridge.err <- unname(coef - exact.b);
    betadf <- rbind(betadf, data.frame(niter=niter, ntrain=ntrain, type="ridge", err=beta.ridge.err, beta=0:15));
    preddf <- rbind(preddf, data.frame(niter=niter, ntrain=ntrain, type="ridge", err=pred.ridge.err));
  }
}

aggregate(err~beta+type+ntrain, data=betadf, FUN=mean);
aggregate(err~beta+type+ntrain, data=betadf, FUN=var);
aggregate(err~beta+type+ntrain, data=betadf, FUN=function(x) mean(x^2));

aggregate(err~type+ntrain, data=preddf, FUN=mean);
aggregate(err~type+ntrain, data=preddf, FUN=var);
aggregate(err~type+ntrain, data=preddf, FUN=function(x) mean(x^2));

