library(tidyverse)

df = read_csv("forestfires.csv")

X = df[,-1]
Y = df[,1]
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



beta0 = matrix(0, nrow  = p, ncol = 1)
beta0[1,1] = mean(Y[,1])
for (i in 2:p) {
        beta0[i,1] = cov(Y[,1], X[,i]) / var(X[,i])
}

L  = function(Y, X, beta){  
        loss = sum(t(Y - X%*%beta)%*%(Y - X%*%beta))
        loss 
}

fit <- optim(par = beta0, fn = L, X = X, Y = Y)





