## dependencies
library(MASS)
library(gam)
library(pROC)
library(MatchIt)
library(WeightIt)

##########################
## utility functions
##########################


## generate correlation matrix
CreateCorrelationMatrix <- function(rho, p) {
  aux1 <- matrix(rep(1:p, p), p, p)
  aux2 <- matrix(rep(1:p, each = p), p, p) 
  rho^abs(aux1 - aux2)
}


GenerateData <- function(n,
                         n.feat,
                         beta.XC = rep(1, n.feat),
                         beta.XY = rep(1, n.feat),
                         rho.X = 0.5,
                         var.C = 1,
                         var.Y = 1,
                         cov.YC = 0) {
  sigma.CC <- var.C
  sigma.YC <- cov.YC
  sigma.YY <- var.Y
  Sig.CY <- matrix(c(sigma.CC, sigma.YC, 
                     sigma.YC, sigma.YY), nr = 2, nc = 2, byrow = TRUE)
  U.CY <- mvrnorm(n, mu = rep(0, 2), Sigma = Sig.CY)
  U.C <- U.CY[, 1]
  U.Y <- U.CY[, 2]
  Sig.X <- CreateCorrelationMatrix(rho.X, n.feat)
  err.X <- mvrnorm(n, mu = rep(0, n.feat), Sigma = Sig.X)
  
  C <- U.C
  Y <- U.Y
  X <- matrix(NA, n, n.feat)
  colnames(X) <- paste("X", seq(n.feat), sep = "")
  for (i in seq(n.feat)) {
    X[, i] <- beta.XC[i] * C + beta.XY[i] * Y + err.X[, i]
  }
  dat <- data.frame(Y, C, X)
  
  dat
}


## simpler implementation 
## (does not work when the response is a factor)
##
CausalityAwareDataLm.0 <- function(dat, 
                                   idx.train,
                                   idx.test,
                                   response.name, 
                                   confounder.name, 
                                   feature.names) {
  n.feat <- length(feature.names)
  rhs.formula <- paste(" ~ ", response.name, " + ", confounder.name, sep = "")
  for (i in seq(n.feat)) {
    my.formula <- as.formula(paste(feature.names[i], rhs.formula, sep = ""))
    fit <- lm(my.formula, data = dat[idx.train,])
    conf.coef.train <- fit$coefficients[3]
    dat[idx.train, feature.names[i]] <- dat[idx.train, feature.names[i]] - dat[idx.train, confounder.name] * conf.coef.train
    dat[idx.test, feature.names[i]] <- dat[idx.test, feature.names[i]] - dat[idx.test, confounder.name] * conf.coef.train
  }
  
  dat
}


## more general implementation
## (response can be numeric or categorical)
## 
CausalityAwareDataLm <- function(dat, 
                                 idx.train,
                                 idx.test,
                                 response.name, 
                                 confounder.name, 
                                 feature.names) {
  n.feat <- length(feature.names)
  rhs.formula <- paste(" ~ ", response.name, " + ", confounder.name, sep = "")
  confounder.mean.train <- mean(dat[idx.train, confounder.name])
  confounder.mean.test <- mean(dat[idx.test, confounder.name])
  for (i in seq(n.feat)) {
    my.formula <- as.formula(paste(feature.names[i], rhs.formula, sep = ""))
    fit <- lm(my.formula, data = dat[idx.train,])
    conf.coef.train <- fit$coefficients[confounder.name]
    aux.train <- predict(fit, newdata = dat[idx.train, c(response.name, confounder.name)], type = "terms")
    dat[idx.train, feature.names[i]] <- dat[idx.train, feature.names[i]] - aux.train[, -1] - conf.coef.train * confounder.mean.train
    aux.test <- predict(fit, newdata = dat[idx.test, c(response.name, confounder.name)], type = "terms") 
    dat[idx.test, feature.names[i]] <- dat[idx.test, feature.names[i]] - aux.test[, -1] - conf.coef.train * confounder.mean.test
  }
  
  dat
}


## compute MSE
LinearRegresssionMSE <- function(dat,
                                 idx.train,
                                 idx.test,
                                 response.name,
                                 feature.names) {
  dat.train <- dat[idx.train,]
  dat.test <- dat[idx.test,]
  my.formula <- as.formula(paste(response.name, " ~ ", paste(feature.names, collapse = "+"), sep = ""))
  fit <- lm(my.formula, data = dat.train)
  y.test.hat <- predict(fit, newdata = dat.test[, feature.names, drop = FALSE])
  y.test <- dat.test[, response.name]
  mse <- mean((y.test.hat - y.test)^2)
  mse
}


FitAnchorRegression <- function(dat, 
                                anchor.names, 
                                input.names, 
                                output.name, 
                                gamma.seq) {
  n <- nrow(dat)
  A <- as.matrix(dat[, anchor.names, drop = FALSE])
  aux <- solve(t(A) %*% A)
  Pi.A <- A %*% aux %*% t(A) 
  n.vals <- length(gamma.seq)
  n.inputs <- length(input.names)
  out <- matrix(NA, n.inputs + 1, n.vals)
  colnames(out) <- gamma.seq
  rownames(out) <- c("intercept", input.names)
  I <- diag(n)
  dat.output <- as.matrix(dat[, output.name])
  dat.input <- as.matrix(dat[, input.names])
  for (i in seq(n.vals)) {
    #cat("gamma = ", gamma.seq[i], "\n")
    W <- I - (1 - sqrt(gamma.seq[i])) * Pi.A
    output.tilde <- W %*% dat.output 
    input.tilde <- W %*% dat.input
    dat.tilde <- data.frame(output.tilde, input.tilde)
    fit <- lm(output.tilde ~ ., data = dat.tilde)
    out[, i] <- fit$coefficients
  }
  
  out
}


ComputeAnchorMSE <- function(dat.tr,
                             dat.ts,
                             anchor.names = "A", 
                             input.names = "X", 
                             output.name = "Y", 
                             gamma.seq = c(0, 1, 1000)) {
  ## fit anchor regression to training data
  b.gamma <- FitAnchorRegression(dat.tr, 
                                 anchor.names, 
                                 input.names, 
                                 output.name, 
                                 gamma.seq)
  
  ## generate model matrix for the test set
  rhs.formula <- paste0("~ ", paste(input.names, collapse = " + "))
  model.matrix.ts <- model.matrix(as.formula(rhs.formula), data = dat.ts)
  
  ## compute predictions
  pred.anchor <- model.matrix.ts %*% b.gamma 
  
  ## compute MSEs
  MSE <- (dat.ts[, output.name] - pred.anchor)^2
  MSE <- apply(MSE, 2, mean)
  
  list(b.gamma = b.gamma, MSE = MSE)
}



## run the experiments
RunShiftExperiments1 <- function(n, 
                                 n.sim, 
                                 n.feat,
                                 response.name, 
                                 confounder.name,
                                 var.C.vals, 
                                 var.Y.vals, 
                                 cov.YC.vals,
                                 beta.XC.range,
                                 beta.XY.range,
                                 rho.X.range,
                                 my.seed,
                                 gamma.seq) {
  n.test <- length(cov.YC.vals)
  dat.test.list <- vector(mode = "list", length = n.test)
  feature.names <- paste("X", seq(n.feat), sep = "")
  MSE.na <- matrix(NA, n.sim, n.test, 
                   dimnames = list(NULL, paste("test set", seq(n.test), sep = " ")))
  MSE.ca <- MSE.na
  MSE.b1 <- MSE.na
  
  n.gamma <- length(gamma.seq)
  MSE.ar <- array(NA, dim = c(n.sim, n.test, n.gamma),
                  dimnames = list(NULL, 
                                  paste("test set", seq(n.test), sep = " "),
                                  gamma.seq))
  
  idx.train <- seq(n)
  idx.test <- seq(n + 1, 2 * n, by = 1)
  
  if (!is.null(my.seed)) {
    set.seed(my.seed)
  }
  for (i in seq(n.sim)) {
    cat(i, "\n")
    
    beta.XC <- runif(n.feat, beta.XC.range[1], beta.XC.range[2])
    beta.XY <- runif(n.feat, beta.XY.range[1], beta.XY.range[2])
    rho.X <- runif(1, rho.X.range[1], rho.X.range[2])
    
    dat.train <- GenerateData(n,
                              n.feat,
                              beta.XC,
                              beta.XY,
                              rho.X,
                              var.C = var.C.vals[1],
                              var.Y = var.Y.vals[1],
                              cov.YC = cov.YC.vals[1])
    dat.train <- data.frame(scale(dat.train, center = TRUE, scale = FALSE))
    
    for (j in seq(n.test)) {
      dat.test <- GenerateData(n,
                               n.feat,
                               beta.XC,
                               beta.XY,
                               rho.X,
                               var.C = var.C.vals[j],
                               var.Y = var.Y.vals[j],
                               cov.YC = cov.YC.vals[j])
      dat.test <- data.frame(scale(dat.test, center = TRUE, scale = FALSE))
      
      ## no adjustmnet
      ## confounded training and test sets
      dat.na <- rbind(dat.train, dat.test)
      
      ## causality-aware
      ## process confounded training and test sets
      dat.ca <- CausalityAwareDataLm(dat.na, 
                                     idx.train, 
                                     idx.test, 
                                     response.name, 
                                     confounder.name,
                                     feature.names)
      
      
      ## baseline 1
      ## unconfounded training (no effect of C on X) but confounded test data
      dat.b1 <- dat.ca
      dat.b1[idx.test, feature.names] <- dat.na[idx.test, feature.names]
      
      MSE.na[i, j] <- LinearRegresssionMSE(dat.na,
                                           idx.train,
                                           idx.test,
                                           response.name = "Y",
                                           feature.name = feature.names) 
      
      MSE.ca[i, j] <- LinearRegresssionMSE(dat.ca,
                                           idx.train,
                                           idx.test,
                                           response.name = "Y",
                                           feature.name = feature.names) 
      
      MSE.b1[i, j] <- LinearRegresssionMSE(dat.b1,
                                           idx.train,
                                           idx.test,
                                           response.name = "Y",
                                           feature.name = feature.names)
      
      aux <- ComputeAnchorMSE(dat.tr = dat.train, 
                              dat.ts = dat.test,
                              anchor.names = "C", 
                              input.names = feature.names, 
                              output.name = "Y", 
                              gamma.seq = gamma.seq)
      MSE.ar[i, j, 1:n.gamma] <- aux$MSE
    }
  }
  
  list(MSE.ca = MSE.ca, 
       MSE.na = MSE.na,
       MSE.b1 = MSE.b1,
       MSE.ar = MSE.ar)
}  




RunShiftExperimentsFailureMode <- function(n, 
                                           n.sim, 
                                           response.name, 
                                           confounder.name,
                                           beta.YC.range,
                                           beta.XC.range,
                                           beta.XY.range,
                                           my.seed,
                                           gamma.seq,
                                           beta.XC.test.values) {
  GenerateDataAR <- function(n,
                             beta.YC = 1,
                             beta.XC = 1,
                             beta.XY = 1) {
    C <- rnorm(n)
    Y <- beta.YC * C + rnorm(n)
    X <- beta.XC * C + beta.XY * Y + rnorm(n)
    dat <- data.frame(C, X, Y)
    
    dat
  }
  
  n.gamma <- length(gamma.seq)
  n.test <- length(beta.XC.test.values)
  dat.test.list <- vector(mode = "list", length = n.test)
  feature.names <- "X"
  MSE.na <- matrix(NA, n.sim, n.test, 
                   dimnames = list(NULL, 
                                   paste("test set", seq(n.test), sep = " ")))
  MSE.ca <- MSE.na
  MSE.b1 <- MSE.na
  n.gamma <- length(gamma.seq)
  MSE.ar <- array(NA, dim = c(n.sim, n.test, n.gamma),
                  dimnames = list(NULL, 
                                  paste("test set", seq(n.test), sep = " "),
                                  gamma.seq))
  
  idx.train <- seq(n)
  idx.test <- seq(n + 1, 2 * n, by = 1)
  
  if (!is.null(my.seed)) {
    set.seed(my.seed)
  }
  for (i in seq(n.sim)) {
    cat(i, "\n")
    
    beta.YC <- runif(1, beta.YC.range[1], beta.YC.range[2])
    beta.XC <- runif(1, beta.XC.range[1], beta.XC.range[2])
    beta.XY <- runif(1, beta.XY.range[1], beta.XY.range[2])
    
    dat.train <- GenerateDataAR(n,
                                beta.YC,
                                beta.XC,
                                beta.XY)
    dat.train$X <- dat.train$X - mean(dat.train$X)
    dat.train$Y <- dat.train$Y - mean(dat.train$Y)
    
    for (j in seq(n.test)) {
      dat.test <- GenerateDataAR(n,
                                 beta.YC,
                                 beta.XC = beta.XC.test.values[j],
                                 beta.XY)
      dat.test$X <- dat.test$X - mean(dat.test$X)
      dat.test$Y <- dat.test$Y - mean(dat.test$Y)
      
      ## no adjustmnet
      ## confounded training and test sets
      dat.na <- rbind(dat.train, dat.test)
      
      ## causality-aware
      ## process confounded training and test sets
      dat.ca <- CausalityAwareDataLm(dat.na, 
                                     idx.train, 
                                     idx.test, 
                                     response.name, 
                                     confounder.name,
                                     feature.names)
      
      ## baseline 1
      ## unconfounded training (no effect of C on X) but confounded test data
      dat.b1 <- dat.ca
      dat.b1[idx.test, feature.names] <- dat.na[idx.test, feature.names]
      
      
      MSE.na[i, j] <- LinearRegresssionMSE(dat.na,
                                           idx.train,
                                           idx.test,
                                           response.name = "Y",
                                           feature.name = feature.names) 
      
      MSE.ca[i, j] <- LinearRegresssionMSE(dat.ca,
                                           idx.train,
                                           idx.test,
                                           response.name = "Y",
                                           feature.name = feature.names) 
      
      MSE.b1[i, j] <- LinearRegresssionMSE(dat.b1,
                                           idx.train,
                                           idx.test,
                                           response.name = "Y",
                                           feature.name = feature.names)
      
      ## anchor regression
      
      aux <- ComputeAnchorMSE(dat.tr = dat.train, 
                              dat.ts = dat.test,
                              anchor.names = "C", 
                              input.names = feature.names, 
                              output.name = "Y", 
                              gamma.seq = gamma.seq)
      MSE.ar[i, j, 1:n.gamma] <- aux$MSE
    }
  }
  
  list(MSE.ca = MSE.ca, 
       MSE.na = MSE.na,
       MSE.b1 = MSE.b1,
       MSE.ar = MSE.ar)
}  





#################################################################
## additional functions for linear vs additive model comparison
#################################################################

CausalityAwareDataGam <- function(dat, 
                                  idx.train,
                                  idx.test,
                                  response.name, 
                                  confounder.name, 
                                  feature.names) {
  n.feat <- length(feature.names)
  if (is.factor(dat[, response.name])) {
    respName <- response.name
  }
  else {
    respName <- paste("s(", response.name, ")", sep = "")
  }
  confName <- paste("s(", confounder.name, ")", sep = "")
  rhs.formula <- paste(" ~ ", respName, " + ", confName, sep = "")
  confounder.mean.train <- mean(dat[idx.train, confounder.name])
  confounder.mean.test <- mean(dat[idx.test, confounder.name])
  for (i in seq(n.feat)) {
    my.formula <- as.formula(paste(feature.names[i], rhs.formula, sep = ""))
    fit <- gam(my.formula, data = dat[idx.train,])
    conf.coef.train <- fit$coefficients[confName]
    aux.train <- predict(fit, newdata = dat[idx.train, c(response.name, confounder.name)], type = "terms")
    dat[idx.train, feature.names[i]] <- dat[idx.train, feature.names[i]] - aux.train[, confName] - conf.coef.train * confounder.mean.train
    aux.test <- predict(fit, newdata = dat[idx.test, c(response.name, confounder.name)], type = "terms") 
    dat[idx.test, feature.names[i]] <- dat[idx.test, feature.names[i]] - aux.test[, confName] - conf.coef.train * confounder.mean.test
  }
  
  dat
}


CreateCorrelationMatrix <- function(rho, p) {
  aux1 <- matrix(rep(1:p, p), p, p)
  aux2 <- matrix(rep(1:p, each = p), p, p) 
  rho^abs(aux1 - aux2)
}


GenerateDataL <- function(n,
                          n.feat,
                          rho.YC,
                          beta.XC,
                          beta.XY,
                          rho.X) {
  Sig <- CreateCorrelationMatrix(rho.YC, 2)
  aux <- mvrnorm(n, mu = c(0, 0), Sigma = Sig)
  C <- aux[, 1]
  Y <- aux[, 2]
  Sig.X <- CreateCorrelationMatrix(rho.X, n.feat)
  err.X <- mvrnorm(n, mu = rep(0, n.feat), Sigma = Sig.X)
  X <- matrix(NA, n, n.feat)
  colnames(X) <- paste("X", seq(n.feat), sep = "")
  for (i in seq(n.feat)) {
    X[, i] <- beta.XC[i] * C + beta.XY[i] * Y + err.X[, i]
  }
  
  data.frame(C, Y, X)
}


GenerateDataNL <- function(n,
                           n.feat,
                           rho.YC,
                           beta.XC,
                           beta.XY,
                           rho.X) {
  Sig <- CreateCorrelationMatrix(rho.YC, 2)
  aux <- mvrnorm(n, mu = c(0, 0), Sigma = Sig)
  C <- aux[, 1]
  Y <- aux[, 2]
  Sig.X <- CreateCorrelationMatrix(rho.X, n.feat)
  err.X <- mvrnorm(n, mu = rep(0, n.feat), Sigma = Sig.X)
  X <- matrix(NA, n, n.feat)
  colnames(X) <- paste("X", seq(n.feat), sep = "")
  for (i in seq(n.feat)) {
    X[, i] <- pmax(beta.XC[i] * C, 0) + pmax(beta.XY[i] * Y, 0) + err.X[, i]
  }
  
  data.frame(C, Y, X)
}



TrainLm <- function(dat.train, 
                    response.name,
                    feature.names) {
  my.formula <- paste0(response.name, " ~ ", paste(feature.names, collapse = " + "))
  lm.model <- lm(as.formula(my.formula), data = dat.train)
  
  lm.model
}


TrainGam <- function(dat.train, 
                     response.name,
                     feature.names) {
  featNames <- paste("s(", feature.names, ")", sep = "")
  my.formula <- paste(response.name, " ~ ", paste(featNames, collapse = " + "), sep = "")
  gam.model <- gam(as.formula(my.formula), data = dat.train)
  
  gam.model
}



RunShiftExperiments2 <- function(n,
                           n.runs,
                           n.feat,
                           rho.YC.train,
                           rho.YC.test,
                           beta.XC.range,
                           beta.XY.range,
                           generate.non.linear.data = TRUE,
                           my.seed = NULL) {
  n.test <- length(rho.YC.test)
  MSEs.lm <- matrix(NA, n.runs, n.test)
  colnames(MSEs.lm) <- paste("test set ", seq(n.test), sep = " ")
  MSEs.lm.ca <- MSEs.lm
  MSEs.gam <- MSEs.lm
  MSEs.gam.ca <- MSEs.lm
  response.name <- "Y"
  confounder.name <- "C"
  feature.names <- paste0("X", seq(n.feat))
  idx.train <- seq(n)
  idx.test <- seq(n + 1, 2 * n)
  
  if (!is.null(my.seed)) {
    set.seed(my.seed)
  }
  for (i in seq(n.runs)) {
    cat("run ", i, "\n")
    aux <- sign(rho.YC.train)
    rho.X <- runif(1, -0.8, 0.8)
    beta.XC <- runif(n.feat, beta.XC.range[1], beta.XC.range[2])
    beta.XY <- runif(n.feat, beta.XY.range[1], beta.XY.range[2])
    if (generate.non.linear.data) {
      dat.train <- GenerateDataNL(n = n,
                                  n.feat = n.feat,
                                  rho.YC = rho.YC.train,
                                  beta.XC = beta.XC,
                                  beta.XY = beta.XY,
                                  rho.X = rho.X)
    }
    else {
      dat.train <- GenerateDataL(n = n,
                                 n.feat = n.feat,
                                 rho.YC = rho.YC.train,
                                 beta.XC = beta.XC,
                                 beta.XY = beta.XY,
                                 rho.X = rho.X)
    }
    
    lm.model <- TrainLm(dat.train, 
                        response.name,
                        feature.names)
    
    gam.model <- TrainGam(dat.train, 
                          response.name,
                          feature.names)
    
    for (j in seq(n.test)) {
      if (generate.non.linear.data) {
        dat.test <- GenerateDataNL(n = n,
                                   n.feat = n.feat,
                                   rho.YC = rho.YC.test[j],
                                   beta.XC = beta.XC,
                                   beta.XY = beta.XY,
                                   rho.X = rho.X)
      }
      else {
        dat.test <- GenerateDataL(n = n,
                                  n.feat = n.feat,
                                  rho.YC = rho.YC.test[j],
                                  beta.XC = beta.XC,
                                  beta.XY = beta.XY,
                                  rho.X = rho.X)
      }
      
      y.hat.lm <- predict(lm.model, newdata = dat.test)
      MSEs.lm[i, j] <- mean((dat.test[, response.name] - y.hat.lm)^2)
      
      y.hat.gam <- predict(gam.model, newdata = dat.test)
      MSEs.gam[i, j] <- mean((dat.test[, response.name] - y.hat.gam)^2)
      
      dat <- rbind(dat.train, dat.test)
      ca.lm <- CausalityAwareDataLm(dat = dat, 
                                    idx.train,
                                    idx.test,
                                    response.name = "Y", 
                                    confounder.name = "C", 
                                    feature.names = feature.names)
      ca.gam <- CausalityAwareDataGam(dat = dat, 
                                      idx.train,
                                      idx.test,
                                      response.name = "Y", 
                                      confounder.name = "C", 
                                      feature.names = feature.names)   
      
      lm.fit.ca <- TrainLm(ca.lm[idx.train,],
                           response.name = "Y",
                           feature.names = feature.names)
      y.hat.lm.ca <- predict(lm.fit.ca, newdata = ca.lm[idx.test,])
      MSEs.lm.ca[i, j] <- mean((ca.lm[idx.test, "Y"] - y.hat.lm.ca)^2)
      
      gam.fit.ca <- TrainGam(ca.gam[idx.train,],
                             response.name = "Y",
                             feature.names = feature.names)
      y.hat.gam.ca <- predict(gam.fit.ca, newdata = ca.gam[idx.test,])
      MSEs.gam.ca[i, j] <- mean((ca.gam[idx.test, "Y"] - y.hat.gam.ca)^2) 
    }
  }
  
  list(MSEs.lm = MSEs.lm,
       MSEs.lm.ca = MSEs.lm.ca,
       MSEs.gam = MSEs.gam,
       MSEs.gam.ca = MSEs.gam.ca)
}


####################################################
## additional functions for the real data analyses
####################################################

GetGlmAuc <- function(dat,
                      idx.train, 
                      idx.test.1, 
                      idx.test.2,
                      label.name, 
                      feature.names,
                      neg.class.name, 
                      pos.class.name) {
  dat <- dat[, c(label.name, feature.names)]
  dat[, label.name] <- factor(as.character(dat[, label.name]), 
                              levels = c(neg.class.name, pos.class.name)) 
  my.formula <- as.formula(paste(label.name, " ~ ", paste(feature.names, collapse = " + ")))
  fit <- glm(my.formula, data = dat[idx.train,], family = "binomial")
  pred.probs.test.1 <- predict(fit, dat[idx.test.1, -1, drop = FALSE], type = "response")
  pred.probs.test.2 <- predict(fit, dat[idx.test.2, -1, drop = FALSE], type = "response")
  pred.probs.train <- predict(fit, dat[idx.train, -1, drop = FALSE], type = "response")
  roc.obj.1 <- roc(dat[idx.test.1, 1], pred.probs.test.1, direction = "<", 
                   levels = c(neg.class.name, pos.class.name))    
  auc.test.1 <- pROC::auc(roc.obj.1)[1]
  roc.obj.2 <- roc(dat[idx.test.2, 1], pred.probs.test.2, direction = "<", 
                   levels = c(neg.class.name, pos.class.name))    
  auc.test.2 <- pROC::auc(roc.obj.2)[1]
  
  list(auc.test.1 = auc.test.1, 
       auc.test.2 = auc.test.2, 
       pred.probs.test.1 = pred.probs.test.1, 
       pred.probs.test.2 = pred.probs.test.2,
       pred.probs.train = pred.probs.train,
       roc.obj.1 = roc.obj.1, 
       roc.obj.2 = roc.obj.2, 
       fit = fit)
}


GetDiscretizedAge <- function(dat, n.levels, level.names = NULL, breaks = NULL) {
  if (is.null(breaks)) {
    if (!is.null(level.names)) {
      out <- cut(dat$age, breaks = n.levels, labels = level.names)
    }
    else {
      levelNames <- paste("age", seq(n.levels), sep = "")
      out <- cut(dat$age, breaks = n.levels, labels = level.names)
    }
  }
  else {
    if (!is.null(level.names)) {
      out <- cut(dat$age, breaks = breaks, labels = level.names)
    }
    else {
      nlevels <- length(breaks) + 1
      levelNames <- paste("age", seq(n.levels), sep = "")
      out <- cut(dat$age, breaks = breaks, labels = level.names)
    }
  }
  
  out
}


## prop has to be larger than 0.5
TrainTest1Test2DataSplit <- function(dat, 
                                     label.name,
                                     confounder.name,
                                     neg.class.name = "FALSE",
                                     pos.class.name = "TRUE",
                                     prop = 0.8) {
  confounder.levels <- sort(unique(dat[, confounder.name]))
  n.confounder.levels <- length(confounder.levels)
  idx.neg.bin <- vector(mode = "list", length = n.confounder.levels)
  idx.pos.bin <- vector(mode = "list", length = n.confounder.levels)
  for (i in seq(n.confounder.levels)) {
    idx.neg.bin[[i]] <- which(dat[, label.name] == neg.class.name & dat[, confounder.name] == confounder.levels[i])
    idx.pos.bin[[i]] <- which(dat[, label.name] == pos.class.name & dat[, confounder.name] == confounder.levels[i])
  }
  n.neg.bin <- unlist(lapply(idx.neg.bin, length))
  n.pos.bin <- unlist(lapply(idx.pos.bin, length))
  
  ## get indexes for test set 2 (where the association with the labels is flipped)
  idx.neg.bin.flipped <- vector(mode = "list", length = n.confounder.levels)
  idx.pos.bin.flipped <- vector(mode = "list", length = n.confounder.levels)
  for (i in seq(n.confounder.levels)) {
    if (n.neg.bin[i] >= n.pos.bin[i]) {
      idx.pos.bin.flipped[[i]] <- sample(idx.pos.bin[[i]], round(n.pos.bin[i] * prop), replace = FALSE)
      idx.neg.bin.flipped[[i]] <- sample(idx.neg.bin[[i]], round(n.pos.bin[i] * (1 - prop)), replace = FALSE)
    }
    if (n.neg.bin[i] < n.pos.bin[i]) {
      idx.neg.bin.flipped[[i]] <- sample(idx.neg.bin[[i]], round(n.neg.bin[i] * prop), replace = FALSE)
      idx.pos.bin.flipped[[i]] <- sample(idx.pos.bin[[i]], round(n.neg.bin[i] * (1 - prop)), replace = FALSE)
    }
  }
  idx.test.2 <- sort(c(unlist(idx.neg.bin.flipped), unlist(idx.pos.bin.flipped)))
  n.test.2 <- length(idx.test.2)
  
  ## remove the indexes from test.set.2, and then split the remaining 
  ## indexes into the training set and test.set.1 
  for (i in seq(n.confounder.levels)) {
    idx.neg.bin[[i]] <- setdiff(idx.neg.bin[[i]], idx.neg.bin.flipped[[i]])
    idx.pos.bin[[i]] <- setdiff(idx.pos.bin[[i]], idx.pos.bin.flipped[[i]])
  }
  n.train.test.1 <- nrow(dat) - n.test.2
  cell.proportions <- table(dat[-idx.test.2, label.name], dat[-idx.test.2, confounder.name])/nrow(dat[-idx.test.2,])
  ## get test.set.1 indexes (this test set has the same distribution as the training set)
  idx.neg.bin.test.1 <- vector(mode = "list", length = n.confounder.levels)
  idx.pos.bin.test.1 <- vector(mode = "list", length = n.confounder.levels)
  for (i in seq(n.confounder.levels)) {
    cell.proportion <- length(idx.neg.bin[[i]])/n.train.test.1
    idx.neg.bin.test.1[[i]] <- sample(idx.neg.bin[[i]], round(cell.proportion * n.test.2), replace = FALSE)
    cell.proportion <- length(idx.pos.bin[[i]])/n.train.test.1
    idx.pos.bin.test.1[[i]] <- sample(idx.pos.bin[[i]], round(cell.proportion * n.test.2), replace = FALSE)
  }
  idx.test.1 <- sort(c(unlist(idx.neg.bin.test.1), unlist(idx.pos.bin.test.1)))
  ## remove the indexes from test.set.1, to get the training set indexes 
  idx.train <- sort(setdiff(c(unlist(idx.neg.bin), unlist(idx.pos.bin)), idx.test.1))
  
  list(idx.train = idx.train,
       idx.test.1 = idx.test.1,
       idx.test.2 = idx.test.2)
}


RunAnalyses <- function(dat, 
                        label.name, 
                        feature.names, 
                        confounder.name,
                        neg.class.name, 
                        pos.class.name,
                        run.seeds,
                        adjustment.method,
                        prop = 0.8) {
  n.runs <- length(run.seeds)
  aucs.1 <- matrix(NA, n.runs, 1, dimnames = list(NULL, "auc"))
  aucs.2 <- aucs.1
  data.set.sizes <- matrix(NA, n.runs, 3)
  colnames(data.set.sizes) <- c("train", "test.1", "test.2")
  
  for (i in seq(n.runs)) {
    cat(i, "\n")
    
    set.seed(run.seeds[i])
    split.dat <- TrainTest1Test2DataSplit(dat, 
                                          label.name,
                                          confounder.name,
                                          neg.class.name,
                                          pos.class.name,
                                          prop)
    
    idx.train <- split.dat$idx.train
    idx.test.1 <- split.dat$idx.test.1
    idx.test.2 <- split.dat$idx.test.2
    
    aux <- ConfounderAdjustment(adjustment.method,
                                dat,
                                idx.train, 
                                idx.test.1,
                                idx.test.2, 
                                label.name, 
                                feature.names, 
                                confounder.name,
                                additional.adjustment.arguments)
    
    adat <- aux$adat
    idx.train <- aux$idx.train
    idx.test.1 <- aux$idx.test.1
    idx.test.2 <- aux$idx.test.2
    
    data.set.sizes[i, 1] <- length(idx.train)
    data.set.sizes[i, 2] <- length(idx.test.1)
    data.set.sizes[i, 3] <- length(idx.test.2)
    
    fit <- GetGlmAuc(adat,
                     idx.train,
                     idx.test.1,
                     idx.test.2,
                     label.name, 
                     feature.names, 
                     neg.class.name, 
                     pos.class.name)
    aucs.1[i, 1] <- fit$auc.test.1
    aucs.2[i, 1] <- fit$auc.test.2
  }
  
  list(aucs.1 = aucs.1,
       aucs.2 = aucs.2,
       data.set.sizes = data.set.sizes)
}


ConfounderAdjustment <- function(adjustment.method,
                                 dat,
                                 idx.train, 
                                 idx.test.1, 
                                 idx.test.2,
                                 label.name, 
                                 feature.names, 
                                 confounder.name,
                                 additional.arguments) {
  aux <- switch(adjustment.method,
                none = NoAdjustment(dat,
                                    idx.train,
                                    idx.test.1,
                                    idx.test.2),
                matching = MatchingAdjustment(dat, 
                                              idx.train, 
                                              idx.test.1, 
                                              idx.test.2,
                                              label.name, 
                                              feature.names, 
                                              confounder.name),
                causality.aware.lm = CausalityAwareAdjustmentLm.2(dat, 
                                                                  idx.train, 
                                                                  idx.test.1, 
                                                                  idx.test.2,
                                                                  response.name = label.name, 
                                                                  feature.names, 
                                                                  confounder.name),
                causality.aware.gam = CausalityAwareAdjustmentGam.2(dat, 
                                                                    idx.train, 
                                                                    idx.test.1, 
                                                                    idx.test.2,
                                                                    response.name = label.name, 
                                                                    feature.names, 
                                                                    confounder.name),
                causality.aware.lm.train.only = CausalityAwareAdjustmentLmTrainOnly.2(dat, 
                                                                                      idx.train, 
                                                                                      idx.test.1, 
                                                                                      idx.test.2,
                                                                                      response.name = label.name, 
                                                                                      feature.names, 
                                                                                      confounder.name),
                causality.aware.gam.train.only = CausalityAwareAdjustmentGamTrainOnly.2(dat, 
                                                                                        idx.train, 
                                                                                        idx.test.1, 
                                                                                        idx.test.2,
                                                                                        response.name = label.name, 
                                                                                        feature.names, 
                                                                                        confounder.name),
                approximate.IPW = ApproximateIpwAdjustment(dat, 
                                                           idx.train,
                                                           idx.test.1,
                                                           idx.test.2,
                                                           label.name,
                                                           confounder.name))
  adat <- aux$adat
  idx.train <- aux$idx.train
  idx.test.1 <- aux$idx.test.1
  idx.test.2 <- aux$idx.test.2
  
  list(adat = adat,
       idx.train = idx.train,
       idx.test.1 = idx.test.1,
       idx.test.2 = idx.test.2)
  
}


NoAdjustment <- function(dat,
                         idx.train,
                         idx.test.1,
                         idx.test.2) {
  adat <- dat
  
  list(adat = adat,
       idx.train = idx.train,
       idx.test.1 = idx.test.1,
       idx.test.2 = idx.test.2)
}


MatchingAdjustment <- function(dat,
                               idx.train,
                               idx.test.1,
                               idx.test.2,
                               label.name,
                               feature.names,
                               confounder.name) {
  OneToOneMatching <- function(dat,
                               label.name,
                               confounder.name) {
    GetMatchedPairs <- function(mdat) {
      aux <- table(mdat$subclass)
      subclasses <- names(aux)
      nclasses <- length(subclasses)
      keep <- c()
      for (i in seq(nclasses)) {
        idx0 <- which(mdat$subclass == subclasses[i] & mdat$numeric.label == 0)
        idx1 <- which(mdat$subclass == subclasses[i] & mdat$numeric.label == 1)
        mdat0 <- mdat[idx0,]
        mdat1 <- mdat[idx1,]
        n0 <- length(idx0) ## number of controls
        n1 <- length(idx1) ## number of cases
        
        ## when we have more controls (n0) than cases (n1),
        ## we select all cases, then randomly select n1 controls 
        if (n0 > n1) {
          keep <- c(keep, idx1, sample(idx0, n1, replace = FALSE))
        }
        
        ## when we have more cases (n1) than controls (n0),
        ## we select all controls, then randomly select n0 cases
        if (n1 > n0) {
          keep <- c(keep, idx0, sample(idx1, n0, replace = FALSE))
        }
        ## when the number of cases is the same as controls
        ## we just select all cases and all controls
        if (n0 == n1) { 
          keep <- c(keep, idx0, idx1)
        }
      }
      mdat[keep,]
    }
    dat$numeric.label <- as.numeric(as.factor(dat[, label.name])) - 1
    my.formula <- as.formula(paste("numeric.label ~ ", confounder.name, sep = ""))
    mm <- matchit(my.formula, data = dat, method = "exact")
    aux <- match.data(mm)
    mdat <- GetMatchedPairs(aux)
    mdat
  }
  ####
  dat.train <- dat[idx.train,]
  dat.test.1 <- dat[idx.test.1,]
  dat.test.2 <- dat[idx.test.2,]
  ## apply adjustment to the training data
  adat.train <- OneToOneMatching(dat.train,
                                 label.name,
                                 confounder.name)
  ## catenate the training and test data
  adat <- rbind(adat.train[, c(label.name, confounder.name, feature.names)], 
                dat.test.1[, c(label.name, confounder.name, feature.names)],
                dat.test.2[, c(label.name, confounder.name, feature.names)])
  
  ## re-assign the training and test set indexes
  idx.train <- seq(nrow(adat.train))
  idx.test.1 <- seq(nrow(adat.train) + 1, nrow(adat.train) + nrow(dat.test.1), by = 1)
  idx.test.2 <- seq(nrow(adat.train) + nrow(dat.test.1) + 1, 
                    nrow(adat.train) + nrow(dat.test.1) + nrow(dat.test.2), by = 1)
  
  list(adat = adat,
       idx.train = idx.train,
       idx.test.1 = idx.test.1,
       idx.test.2 = idx.test.2)
}


CausalityAwareAdjustmentLm.2 <- function(dat, 
                                         idx.train,
                                         idx.test.1,
                                         idx.test.2,
                                         response.name, 
                                         feature.names,
                                         confounder.name) {
  ## confounder is a numeric variable
  ## response can be a numeric or factor variable
  CausalityAwareAdjustmentLm1 <- function(dat, 
                                          idx.train,
                                          idx.test.1,
                                          idx.test.2,
                                          response.name, 
                                          feature.names,
                                          confounder.name) {
    n.feat <- length(feature.names)
    rhs.formula <- paste(" ~ ", response.name, " + ", confounder.name, sep = "")
    confounder.mean.train <- mean(dat[idx.train, confounder.name], na.rm = TRUE)
    confounder.mean.test.1 <- mean(dat[idx.test.1, confounder.name])
    confounder.mean.test.2 <- mean(dat[idx.test.2, confounder.name])
    for (i in seq(n.feat)) {
      my.formula <- as.formula(paste(feature.names[i], rhs.formula, sep = ""))
      fit <- lm(my.formula, data = dat[idx.train,])
      conf.coef.train <- fit$coefficients[confounder.name]
      aux.train <- predict(fit, newdata = dat[idx.train, c(response.name, confounder.name)], type = "terms")
      dat[idx.train, feature.names[i]] <- dat[idx.train, feature.names[i]] - aux.train[, -1] - conf.coef.train * confounder.mean.train
      ####
      aux.test.1 <- predict(fit, newdata = dat[idx.test.1, c(response.name, confounder.name)], type = "terms") 
      dat[idx.test.1, feature.names[i]] <- dat[idx.test.1, feature.names[i]] - aux.test.1[, -1] - conf.coef.train * confounder.mean.test.1
      ####
      aux.test.2 <- predict(fit, newdata = dat[idx.test.2, c(response.name, confounder.name)], type = "terms") 
      dat[idx.test.2, feature.names[i]] <- dat[idx.test.2, feature.names[i]] - aux.test.2[, -1] - conf.coef.train * confounder.mean.test.2    
    }
    list(adat = dat,
         idx.train = idx.train,
         idx.test.1 = idx.test.1,
         idx.test.2 = idx.test.2)
  }
  ## confounder is a factor variable
  ## response is a factor variable
  CausalityAwareAdjustmentLm2 <- function(dat, 
                                          idx.train,
                                          idx.test.1,
                                          idx.test.2,
                                          label.name, 
                                          feature.names,
                                          confounder.name) {
    dat.train <- dat[idx.train,]
    dat.test.1 <- dat[idx.test.1,]
    dat.test.2 <- dat[idx.test.2,]
    train.model.matrix <- model.matrix(as.formula(paste("~", confounder.name, sep = "")), 
                                       data = dat.train)[, -1, drop = FALSE] ## remove the intercept column
    test.model.matrix.1 <- model.matrix(as.formula(paste("~", confounder.name, sep = "")), 
                                        data = dat.test.1)[, -1, drop = FALSE] ## remove the intercept column
    test.model.matrix.2 <- model.matrix(as.formula(paste("~", confounder.name, sep = "")), 
                                        data = dat.test.2)[, -1, drop = FALSE] ## remove the intercept column
    n.feat <- length(feature.names)
    rhs.formula <- paste(" ~ ", label.name, " + ", confounder.name, sep = "")
    n.levels.label <- nlevels(dat[, label.name])
    for (i in seq(n.feat)) {
      my.formula <- as.formula(paste(feature.names[i], rhs.formula, sep = ""))
      fit <- lm(my.formula, data = dat.train)
      conf.coef.train <- fit$coefficients[-c(1:n.levels.label)] ## get only the confounder coefficient
      dat.train[, feature.names[i]] <- dat.train[, feature.names[i]] - as.numeric(train.model.matrix %*% conf.coef.train)
      dat.test.1[, feature.names[i]] <- dat.test.1[, feature.names[i]] - as.numeric(test.model.matrix.1 %*% conf.coef.train)
      dat.test.2[, feature.names[i]] <- dat.test.2[, feature.names[i]] - as.numeric(test.model.matrix.2 %*% conf.coef.train)
    }
    adat <- rbind(dat.train, dat.test.1, dat.test.2)
    idx.train <- seq(nrow(dat.train))
    idx.test.1 <- seq(nrow(dat.train) + 1, nrow(dat.train) + nrow(dat.test.1), by = 1)
    idx.test.2 <- seq(nrow(dat.train) + nrow(dat.test.1) + 1, 
                      nrow(dat.train) + nrow(dat.test.1) + nrow(dat.test.2), by = 1)
    list(adat = adat,
         idx.train = idx.train,
         idx.test.1 = idx.test.1,
         idx.test.2 = idx.test.2)
  }
  if (is.factor(dat[, confounder.name])) {
    out <- CausalityAwareAdjustmentLm2(dat, 
                                       idx.train,
                                       idx.test.1,
                                       idx.test.2,
                                       label.name = response.name, 
                                       feature.names,
                                       confounder.name)
  }
  else {
    out <- CausalityAwareAdjustmentLm1(dat, 
                                       idx.train,
                                       idx.test.1,
                                       idx.test.2,
                                       response.name, 
                                       feature.names,
                                       confounder.name)
  }
  
  out
}


CausalityAwareAdjustmentLmTrainOnly.2 <- function(dat, 
                                                  idx.train,
                                                  idx.test.1,
                                                  idx.test.2,
                                                  response.name, 
                                                  feature.names,
                                                  confounder.name) {
  ## confounder is a numeric variable
  ## response can be a numeric or factor variable
  CausalityAwareAdjustmentLmTrainOnly1 <- function(dat, 
                                                   idx.train,
                                                   idx.test.1,
                                                   idx.test.2,
                                                   response.name, 
                                                   feature.names,
                                                   confounder.name) {
    n.feat <- length(feature.names)
    rhs.formula <- paste(" ~ ", response.name, " + ", confounder.name, sep = "")
    confounder.mean.train <- mean(dat[idx.train, confounder.name], na.rm = TRUE)
    for (i in seq(n.feat)) {
      my.formula <- as.formula(paste(feature.names[i], rhs.formula, sep = ""))
      fit <- lm(my.formula, data = dat[idx.train,])
      conf.coef.train <- fit$coefficients[confounder.name]
      aux.train <- predict(fit, newdata = dat[idx.train, c(response.name, confounder.name)], type = "terms")
      dat[idx.train, feature.names[i]] <- dat[idx.train, feature.names[i]] - aux.train[, -1] - conf.coef.train * confounder.mean.train
    }
    list(adat = dat,
         idx.train = idx.train,
         idx.test.1 = idx.test.1,
         idx.test.2 = idx.test.2)
  }
  ## confounder is a factor variable
  ## response is a factor variable
  CausalityAwareAdjustmentLmTrainOnly2 <- function(dat, 
                                                   idx.train,
                                                   idx.test.1,
                                                   idx.test.2,
                                                   label.name, 
                                                   feature.names,
                                                   confounder.name) {
    dat.train <- dat[idx.train,]
    dat.test.1 <- dat[idx.test.1,]
    dat.test.2 <- dat[idx.test.2,]
    train.model.matrix <- model.matrix(as.formula(paste("~", confounder.name, sep = "")), 
                                       data = dat.train)[, -1, drop = FALSE] ## remove the intercept column
    n.feat <- length(feature.names)
    rhs.formula <- paste(" ~ ", label.name, " + ", confounder.name, sep = "")
    n.levels.label <- nlevels(dat[, label.name])
    for (i in seq(n.feat)) {
      my.formula <- as.formula(paste(feature.names[i], rhs.formula, sep = ""))
      fit <- lm(my.formula, data = dat.train)
      conf.coef.train <- fit$coefficients[-c(1:n.levels.label)] ## get only the confounder coefficient
      dat.train[, feature.names[i]] <- dat.train[, feature.names[i]] - as.numeric(train.model.matrix %*% conf.coef.train)
    }
    adat <- rbind(dat.train, dat.test.1, dat.test.2)
    idx.train <- seq(nrow(dat.train))
    idx.test.1 <- seq(nrow(dat.train) + 1, nrow(dat.train) + nrow(dat.test.1), by = 1)
    idx.test.2 <- seq(nrow(dat.train) + nrow(dat.test.1) + 1, 
                      nrow(dat.train) + nrow(dat.test.1) + nrow(dat.test.2), by = 1)
    list(adat = adat,
         idx.train = idx.train,
         idx.test.1 = idx.test.1,
         idx.test.2 = idx.test.2)
  }
  if (is.factor(dat[, confounder.name])) {
    out <- CausalityAwareAdjustmentLmTrainOnly2(dat, 
                                                idx.train,
                                                idx.test.1,
                                                idx.test.2,
                                                label.name = response.name, 
                                                feature.names,
                                                confounder.name)
  }
  else {
    out <- CausalityAwareAdjustmentLmTrainOnly1(dat, 
                                                idx.train,
                                                idx.test.1,
                                                idx.test.2,
                                                response.name, 
                                                feature.names,
                                                confounder.name)
  }
  
  out
}


CausalityAwareAdjustmentGam.2 <- function(dat, 
                                          idx.train,
                                          idx.test.1,
                                          idx.test.2,
                                          response.name, 
                                          feature.names,
                                          confounder.name) {
  n.feat <- length(feature.names)
  confName <- paste("s(", confounder.name, ")", sep = "")
  rhs.formula <- paste(" ~ ", response.name, " + ", confName, sep = "")
  confounder.mean.train <- mean(dat[idx.train, confounder.name])
  confounder.mean.test.1 <- mean(dat[idx.test.1, confounder.name])
  confounder.mean.test.2 <- mean(dat[idx.test.2, confounder.name])
  for (i in seq(n.feat)) {
    my.formula <- as.formula(paste(feature.names[i], rhs.formula, sep = ""))
    fit <- gam(my.formula, data = dat[idx.train,])
    conf.coef.train <- fit$coefficients[confName]
    aux.train <- predict(fit, newdata = dat[idx.train, c(response.name, confounder.name)], type = "terms")
    dat[idx.train, feature.names[i]] <- dat[idx.train, feature.names[i]] - aux.train[, confName] - conf.coef.train * confounder.mean.train
    ####
    aux.test.1 <- predict(fit, newdata = dat[idx.test.1, c(response.name, confounder.name)], type = "terms") 
    dat[idx.test.1, feature.names[i]] <- dat[idx.test.1, feature.names[i]] - aux.test.1[, confName] - conf.coef.train * confounder.mean.test.1
    ####
    aux.test.2 <- predict(fit, newdata = dat[idx.test.2, c(response.name, confounder.name)], type = "terms") 
    dat[idx.test.2, feature.names[i]] <- dat[idx.test.2, feature.names[i]] - aux.test.2[, confName] - conf.coef.train * confounder.mean.test.2
  }
  
  list(adat = dat,
       idx.train = idx.train, 
       idx.test.1 = idx.test.1,
       idx.test.2 = idx.test.2)
}


CausalityAwareAdjustmentGamTrainOnly.2 <- function(dat, 
                                                   idx.train,
                                                   idx.test.1,
                                                   idx.test.2,
                                                   response.name, 
                                                   feature.names,
                                                   confounder.name) {
  n.feat <- length(feature.names)
  confName <- paste("s(", confounder.name, ")", sep = "")
  rhs.formula <- paste(" ~ ", response.name, " + ", confName, sep = "")
  confounder.mean.train <- mean(dat[idx.train, confounder.name])
  for (i in seq(n.feat)) {
    my.formula <- as.formula(paste(feature.names[i], rhs.formula, sep = ""))
    fit <- gam(my.formula, data = dat[idx.train,])
    conf.coef.train <- fit$coefficients[confName]
    aux.train <- predict(fit, newdata = dat[idx.train, c(response.name, confounder.name)], type = "terms")
    dat[idx.train, feature.names[i]] <- dat[idx.train, feature.names[i]] - aux.train[, confName] - conf.coef.train * confounder.mean.train
  }
  
  list(adat = dat,
       idx.train = idx.train, 
       idx.test.1 = idx.test.1,
       idx.test.2 = idx.test.2)
}


ApproximateIpwAdjustment <- function(dat, 
                                     idx.train,
                                     idx.test.1,
                                     idx.test.2,
                                     label.name,
                                     confounder.name) { 
  PseudoPopulation <- function(dat, weights) {
    inv.weights <- round(weights)
    inv.weights[inv.weights == 0] <- 1
    idx <- rep(seq(nrow(dat)), times = inv.weights)
    dat[idx,]
  }
  WeightsPs <- function(dat,
                        label.name,
                        confounder.name) {
    my.formula <- as.formula(paste(label.name, " ~ ", paste(confounder.name, collapse = " + ")))
    aux <- WeightIt::weightit(my.formula, data = dat, estimand = "ATE", method = "ps")
    aux$weights
  }
  dat.train <- dat[idx.train,]
  dat.test.1 <- dat[idx.test.1,]
  dat.test.2 <- dat[idx.test.2,]
  
  w.train <- WeightsPs(dat.train,
                       label.name,
                       confounder.name)
  adat.train <- PseudoPopulation(dat.train, w.train)
  adat <- rbind(adat.train[, c(label.name, confounder.name, feature.names)], 
                dat.test.1[, c(label.name, confounder.name, feature.names)],
                dat.test.2[, c(label.name, confounder.name, feature.names)])
  
  idx.train <- seq(nrow(adat.train))
  idx.test.1 <- seq(nrow(adat.train) + 1, nrow(adat.train) + nrow(dat.test.1), by = 1)
  idx.test.2 <- seq(nrow(adat.train) + nrow(dat.test.1) + 1, 
                    nrow(adat.train) + nrow(dat.test.1) + nrow(dat.test.2), by = 1)
  
  list(adat = adat, 
       idx.train = idx.train, 
       idx.test.1 = idx.test.1,
       idx.test.2 = idx.test.2)
}


CatenateAucs1 <- function(file.names, n_runs = 100, method.names) {
  n_files <- length(file.names)
  aucs <- matrix(NA, n_runs, n_files)
  colnames(aucs) <- method.names
  for (i in seq(n_files)) {
    load(file.names[i])
    aucs[, i] <- out$aucs.1
  }
  
  aucs
}


CatenateAucs2 <- function(file.names, n_runs = 100, method.names) {
  n_files <- length(file.names)
  aucs <- matrix(NA, n_runs, n_files)
  colnames(aucs) <- method.names
  for (i in seq(n_files)) {
    load(file.names[i])
    aucs[, i] <- out$aucs.2
  }
  
  aucs
}


####################################################
## plotting functions for generating paper figures
####################################################

MyPlot1 <- function(e, my.ylim.1 = NULL, my.ylim.2 = NULL) {
  BestMSE <- function(x) {
    n.sim <- dim(x)[1]
    n.test <- dim(x)[2]
    best.mse <- matrix(NA, n.sim, n.test)
    colnames(best.mse) <- colnames(x[,, 1])
    for (i in seq(n.sim)) {
      aux <- x[i,,]
      for (j in seq(n.test)) {
        idx <- which.min(aux[j,])
        best.mse[i, j] <- aux[j, idx]
      }
    }
    best.mse
  }
  jit <- 0.15
  my.lwd <- 2
  my.lwd2 <- 1
  head.length <- 0
  
  col.na <- "darkorange"
  col.b1 <- "darkred"
  col.ca <- "darkblue"
  col.ar <- "cyan"
  
  lty.na <- 4
  lty.b1 <- 3
  lty.ca <- 1
  lty.ar <- 2
  lty <- 1
  ltype <- "l"
  my.cex <- 1
  my.pch <- 20
  
  mse.ar <- BestMSE(x = e$MSE.ar)
  
  err.ca <- apply(e$MSE.ca, 1, sd)
  err.b1 <- apply(e$MSE.b1, 1, sd)
  err.na <- apply(e$MSE.na, 1, sd)
  err.ar <- apply(mse.ar, 1, sd)
  errs <- data.frame(err.ca, err.b1, err.ar, err.na)
  names(errs) <- c("causal. aw.", "baseline", "anchor regr.", "no adjust.")
  
  m.ca <- apply(e$MSE.ca, 2, mean)
  m.b1 <- apply(e$MSE.b1, 2, mean)
  m.na <- apply(e$MSE.na, 2, mean)
  m.ar <- apply(mse.ar, 2, mean)
  
  sd.ca <- apply(e$MSE.ca, 2, sd)
  sd.b1 <- apply(e$MSE.b1, 2, sd)
  sd.na <- apply(e$MSE.na, 2, sd)
  sd.ar <- apply(mse.ar, 2, sd)
  
  if (is.null(my.ylim.1)) {
    my.ylim.1 <- c(0, max(unlist(e)))
  }
  if (is.null(my.ylim.2)) {
    my.ylim.2 <- c(0, max(errs))
  }
  cex1 <- 1
  n.test <- ncol(e$MSE.ca)
  xaxis <- seq(n.test)
  par(mfrow = c(1, 2), mar = c(6, 3, 1, 0.5), mgp = c(2, 0.75, 0))
  ##################
  plot(xaxis, m.na, ylim = my.ylim.1, type = "n", xaxt = "n", xlab = "", ylab = "MSE",
       main = "predictive perf.")
  axis(side = 1, at = seq(n.test), labels = paste0("test set ", seq(n.test)), las = 2)
  #####
  points(xaxis - jit, m.na, cex = my.cex, col = col.na, pch = my.pch)
  lines(xaxis - jit, m.na, type = ltype, col = col.na, lwd = my.lwd, lty = lty.na)
  arrows(xaxis - jit, m.na - sd.na, xaxis - jit, m.na + sd.na, col = col.na, 
         code = 3, angle = 90, length = head.length, lwd = my.lwd2, lty = lty)
  #####
  points(xaxis - jit/2, m.b1, cex = my.cex, col = col.b1, pch = my.pch)
  lines(xaxis - jit/2, m.b1, type = ltype, col = col.b1, lwd = my.lwd, lty = lty.b1)
  arrows(xaxis - jit/2, m.b1 - sd.b1, xaxis - jit/2, m.b1 + sd.b1, col = col.b1, 
         code = 3, angle = 90, length = head.length, lwd = my.lwd2, lty = lty)
  #####
  points(xaxis + jit/2, m.ar, cex = my.cex, col = col.ar, pch = my.pch)
  lines(xaxis + jit/2, m.ar, type = ltype, col = col.ar, lwd = my.lwd, lty = lty.ar)
  arrows(xaxis + jit/2, m.ar - sd.ar, xaxis + jit/2, m.ar + sd.ar, col = col.ar, 
         code = 3, angle = 90, length = head.length, lwd = my.lwd2, lty = lty)
  #####
  points(xaxis + jit, m.ca, cex = my.cex, col = col.ca, pch = my.pch)
  lines(xaxis + jit, m.ca, type = ltype, col = col.ca, lwd = my.lwd, lty = lty.ca)
  arrows(xaxis + jit, m.ca - sd.ca, xaxis + jit, m.ca + sd.ca, col = col.ca, 
         code = 3, angle = 90, length = head.length, lwd = my.lwd2, lty = lty)
  mtext(side = 3, at = 0, text = "(a)", cex = cex1)
  legend("topleft", legend = c("causality-aware", "baseline", "anchor regression", "no adjustment"),
         text.col = c(col.ca, col.b1, col.ar, col.na), bty = "n", 
         lty = c(lty.ca, lty.b1, lty.ar, lty.na),
         col = c(col.ca, col.b1, col.ar, col.na))
  ######
  boxplot(errs, border = c(col.ca, col.b1, col.ar, col.na),
          las = 2, ylab = "error", cex = 0.1,
          main = "stability error", ylim = my.ylim.2, col = "white")
  legend("topleft", legend = c("causality-aware", "baseline", "anchor regression", "no adjustment"),
         text.col = c(col.ca, col.b1, col.ar, col.na), bty = "n")
  mtext(side = 3, at = 0, text = "(b)", cex = cex1)
}


MyPlot2 <- function(MSEs.lm,
                    MSEs.lm.ca,
                    MSEs.gam,
                    MSEs.gam.ca,
                    ylim1 = c(0, 1.2),
                    ylim2 = c(0, 0.6),
                    jit = 0.1) {
  col.1 <- "darkorange"
  col.2 <- "darkblue"
  col.3 <- "purple"
  col.4 <- "darkgreen"
  n.test <- ncol(MSEs.lm)
  s.lm <- apply(MSEs.lm, 1, sd)
  s.lm.ca <- apply(MSEs.lm.ca, 1, sd)
  s.gam <- apply(MSEs.gam, 1, sd)
  s.gam.ca <- apply(MSEs.gam.ca, 1, sd)
  ss <- data.frame(s.lm, s.lm.ca, s.gam, s.gam.ca)
  names(ss) <- c("lin. model", "c. a. lin. m.", "add. model", "c. a. add. m.")
  
  par(mfrow = c(1, 3), mar = c(6, 3.5, 2, 1), mgp = c(2.3, 0.75, 0))
  boxplot(MSEs.lm, border = col.1, col = "white", ylim = ylim1, xaxt = "n", at = seq(n.test) - jit,
          main = "linear model", ylab = "MSE")
  boxplot(MSEs.lm.ca, border = col.2, col = rgb(0, 0, 1, 0.25), xaxt = "n", at = seq(n.test) + jit, add = TRUE)
  axis(side = 1, labels = colnames(MSEs.lm), at = seq(n.test), las = 2)
  legend("topleft", legend = c("linear model (no adjustment)", "causality-aware lin. model"), 
         text.col = c(col.1, col.2), bty = "n", 
         fill = c("white", rgb(0, 0, 1, 0.25)),
         border = c(col.1, col.2))
  mtext(side = 3, "(a)", at = 0)
  ####
  boxplot(MSEs.gam, border = col.3, col = "white", ylim = ylim1, xaxt = "n", at = seq(n.test) - jit,
          main = "additive model", ylab = "MSE")
  boxplot(MSEs.gam.ca, border = col.4, col = rgb(0, 1, 0, 0.25), xaxt = "n", at = seq(n.test) + jit, add = TRUE)
  axis(side = 1, labels = colnames(MSEs.lm), at = seq(n.test), las = 2)
  legend("topleft", legend = c("additive model (no adjust.)", "causality-aware add. model"), 
         text.col = c(col.3, col.4), bty = "n",
         fill = c("white", rgb(0, 1, 0, 0.25)),
         border = c(col.3, col.4))
  mtext(side = 3, "(b)", at = 0)
  ####
  boxplot(ss, ylab = "error", main = "stability error", ylim = ylim2,
          border = c(col.1, col.2, col.3, col.4),
          col = c("white", rgb(0, 0, 1, 0.25), "white", rgb(0, 1, 0, 0.25)),
          las = 2)
  mtext(side = 3, "(c)", at = 0)
  legend("topleft", 
         legend = c("linear model (no adjustment)", "causality-aware lin. model", 
                    "additive model (no adjust.)", "causality-aware add. model"), 
         text.col = c(col.1, col.2, col.3, col.4), 
         fill = c("white", rgb(0, 0, 1, 0.25), "white", rgb(0, 1, 0, 0.25)),
         border = c(col.1, col.2, col.3, col.4),
         bty = "n")
  par(mfrow = c(1, 3), mar = c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))
}

