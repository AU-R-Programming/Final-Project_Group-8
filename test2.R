library(tidyverse)

df = read_csv("crop.data.csv")

X = df[,-4]
Y = df[,4]
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

# Define parameters

Y = shapeY(Y)
X = shapeX(X)
n <- length(Y)
p <- dim(X)[2]
df <- n - p


beta0 = matrix(0, nrow  = p, ncol = 1)
beta0[1,1] = mean(Y[,1])
for (i in 2:p) {
    beta0[i,1] = cov(Y[,1], X[,i]) / var(X[,i])
}

L  = function(Y, X, beta){  
    loss = sum(t(Y - X%*%beta)%*%(Y - X%*%beta))
    loss 
}

fit = optim(par = beta0, fn = L, X = X, Y = Y)
beta.hat = fit$par

resid <- Y - X%*%beta.hat
sigma2.hat <- (1/(n -p))*t(resid)%*%resid

# Estimate of the variance of the estimated beta from Eq. (6.2)
var.beta <- as.vector(sigma2.hat)*solve(t(X)%*%X)

# Estimate of the confidence interval based on alpha
alpha = 0.5
quant <- 1 - alpha/2
ci.beta <- c(beta.hat - qnorm(p = quant)*sqrt(var.beta), beta.hat + qnorm(p = quant)*sqrt(var.beta))


