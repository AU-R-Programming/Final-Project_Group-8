#' @title "Linear Regression Model"
#' 
#' 
#' @description Contains a linear regression function that also evaluates parameters such as beta hat, confidence 
#' intervals, and plots
#' @name my_lm
#' @param response Also known as beta, an estimation of the coefficients in the model
#' @param predictors Also known as sigma2, an estimation of the error  
#' @importFrom tidyverse ggplot2
#' @export my_lm 
library(ggplot2)

# add a column with ones to the feature matrix
shapeX = function(X) {
  X = as.matrix(X)
  X = cbind(1, X)
  X
}

shapeY = function(Y) {
  Y = as.matrix(Y)
  Y
}

# least squares estimators that minimizer
L  = function(response, predictors, beta){  
  loss = sum(t(response - predictors%*%beta)%*%(response - predictors%*%beta))
  loss 
}

my_lm = function(response, predictors, alpha) {
  
  # Make sure data formats are appropriate
  response <- shapeY(response)
  predictors <- shapeX(predictors)
  
  
  # Define parameters
  n <- length(response)
  p <- dim(predictors)[2]
  df <- n - p
  
  # Define the initial beta
  beta0 = matrix(0, nrow  = p, ncol = 1)
  beta0[1,1] = mean(response[,1])
  for (i in 2:p) {
    beta0[i,1] = cov(response[,1], predictors[,i]) / var(predictors[,i])
  }
  fit = optim(par = as.vector(beta0), fn = L, predictors = predictors, response = response)
  beta.hat = as.matrix(fit$par)
  
  # Estimate of the residual variance (sigma2) from Eq. (6.3)
  # Compute residuals
  resids <- response - predictors%*%as.vector(beta.hat) 
  sigma2.hat <- (1/df)*(t(resids)%*%resids)
  
  # Estimate of the variance of the estimated beta from Eq. (6.2)
  var.beta_mat <- as.vector(sigma2.hat)*solve(t(predictors)%*%predictors)
  var.beta = diag(var.beta_mat)
  
  # Estimate of the confidence interval based on alpha
  quant <- 1 - alpha/2
  ci.beta <- c(beta.hat - qnorm(p = quant)*sqrt(var.beta), beta.hat + qnorm(p = quant)*sqrt(var.beta))
  
  # Get the prediction values of the responses
  ypred<-predictors%*%as.vector(beta.hat) 
  
  # Calculate the R_squared values through SSE and SST
  SSM<-sum((ypred-mean(response))^2)
  SSE<-sum((response-ypred)^2)
  SST<-sum((response-mean(response))^2)
  R_squared<-1-(SSE/SST)
  
  # Calculate the Cp value
  Cp<-SSE+2*p*as.numeric(sigma2.hat)
  
  # Calculate the F_statistics and P-value
  DFM<-p-1
  DFE<-n-p
  F_statistics<-(SSM/DFM)/(SSE/DFE)
  P_value<-pf(F_statistics,DFM,DFE)
  
  #Adding other measurements for the model accuracy
  RMSE<-sqrt(sum(response-ypred)^2/n)
  MAE<-sum(abs(response-ypred))/n
  MAPE<-sum(abs((response-ypred)/response))/n
  
  #Plotting the Residual vs. fitted plots
  r_vs_f<-plot(resids,ypred)
  
  #Plotting the q-q plots
  qqnorm(resids)
  
  #Plotting Histogram
  plot(density(resids))
  
  
  # Return all estimated values
  return(list(beta = beta.hat, sigma2 = sigma2.hat, 
              variance_beta = var.beta, ci = ci.beta, Rsquared = R_squared, Cp = Cp, F_statistics = F_statistics, PValue = P_value, RMSE = RMSE, MAE = MAE, MAPE = MAPE))
}

