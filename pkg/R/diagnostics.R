library(emulator)

## Example used in Section 4.1 of O'Hagan/Bastos diagnostics paper
simulator <- function(input){
  
  x1 = input[,1]
  x2 = input[,2]
  
  f = (1 - exp(-1/2*x2)) *
    (2300*x1^3 + 1900*x1^2 + 2092*x1 + 60)/(100*x1^3 + 500*x1^2 + 4*x1 + 20)

  return(f)
}


## Make predictions using the emulator for new input points
fit.emulator <- function(train.data, predict.data, hyperparams){
  
  scale = hyperparams[-length(hyperparams)] # same number of scale parameters as dimension of input
  sigmasq = hyperparams[length(hyperparams)] # one variance parameter
  
  ## Calc covariance matrices
  K.train = sigmasq * corr.matrix(train.data$input, scales = scale) 
  K.diag = sigmasq * corr.matrix(train.data$input, predict.data$input, scales = scale)
  K.predict = sigmasq * corr.matrix(predict.data$input, predict.data$input, scales = scale)

  ## Calc fitted mean and variance
  fitted.mean = K.diag %*% solve(K.train) %*% train.data$output
  fitted.var = K.predict -  K.diag %*% solve(K.train) %*% t(K.diag) 
  
  emulator = list()
  emulator$mean = fitted.mean
  emulator$var = fitted.var
  
  return(emulator)
}


## Calculate individual prediction errors, pit values and Mahalanobis distance
diagnostics <- function(verify,emulator){
  
  pred.errs = (verify$output - emulator$mean) / sqrt(diag(emulator$var))

  pit = pnorm(pred.errs)
  
  mahalanobis = t(verify$output - emulator$mean) %*%
    solve(emulator$var) %*% (verify$output - emulator$mean)
  
  return(list(pred.errs = pred.errs, pit = pit, mahalanobis = mahalanobis))
}


## Training data to build emulator
train.data = list()
train.data$input = latin.hypercube(20,2)
train.data$output = simulator(train.data$input)

## Verification data
verify.data = list()
verify.data$input = latin.hypercube(25,2)
verify.data$output = simulator(verify.data$input)

## Calculate correlation matrices
## N.B. corr.matrix function uses different parameterisation to O'Hagan
hyperparams = c(1/0.2421^2, 1/0.4240^2, 3.3316)

## Fit emulator
emulator = fit.emulator(train.data, verify.data, hyperparams)

## Calculate diagnostics
d = diagnostics(verify.data, emulator)
