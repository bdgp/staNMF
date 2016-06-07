#######################################################################
# The following R code perform the LASSO+OLS procedure to generates sPP 
# coefficients for each expression pattern. The final product, alpha, 
# should be the same as the matrix stored in sPPCoefPatterns.csv
# Last update: March 28, 2016 by Siqi Wu
#######################################################################


source('loadFunctions.R');

# Make sure that you have the packages glmnet and nnls for some 
library(glmnet);
library(nnls);


# Import the spatial gene expression data matrix
X = read.table(file = '../expressionPatterns.csv',header = TRUE,
    		    check.names = FALSE,sep=',');
X = as.matrix(X);

# Import the learned dictionary with 21 principal patterns (PP)
D = read.table(file = '../PP.csv',header = TRUE,check.names = FALSE,sep = ',');
D = as.matrix(D);

K = ncol(D);

# 1. Model selection using LASSO and cross-validation

beta = matrix(0,nrow = K, ncol = ncol(X));
set.seed(215)
for (i in 1:ncol(X)){
    print(i)
    y = X[,i];
    res = cv.glmnet(D,y,family = 'gaussian',nfold = 10, nlambda = 100,intercept = TRUE,
        type.measure = 'mse', standardize = FALSE, keep = FALSE, lower.limits =0);
    idxTemp = which(res$lambda == res$lambda.1se);
    beta[,i] = res$glmnet.fit$beta[,idxTemp];
}

# 2. Refit using nonnegative least squares (NLS)
result = olsLassoFit(X,D,as.matrix(beta),intercept = TRUE);
alpha = result$alphaOLS;