###########################################################################
# Tutorial for data import, sptial expression pattern visualization,
# image reconstruction using sPP coefficients and gene categorization. 
# Last update: March 28, 2016 by Siqi Wu
###########################################################################

# 1. Import data 

# Import the spatial gene expression data matrix
X = read.table(file = '../expressionPatterns.csv',header = TRUE,check.names = FALSE,sep=',');
geneNames = colnames(X);

# Import the sparse PP (sPP) coefficients for all 1640 images
sPPCoef = read.table(file = '../sPPcoefPatterns.csv',header = TRUE,check.names = FALSE,sep = ',');

# Import the learned dictionary with 21 principal patterns (PP)
D = read.table(file = '../PP.csv',header = TRUE,check.names = FALSE,sep = ',');

# Import the elliptical template for spatial gene expression visualization
template = read.csv(file = './utilities/embryoTemplate.csv',sep = ',',header = TRUE,row.names = 1);


# 2. Data visualization

source('./loadFunctions.R')

# Display the first 25 spatial gene expression patterns in the data matrix X
imageBatchDisplay(X,template = template,nrow = 5, ncol = 5, font = 1)

# Display the 21 learned PP
imageBatchDisplay(D,template = template, nrow = 3, ncol = 7, font = 1, noNumber = T)



# 3. Reconstruct expression patterns using the 21 PP and the sPP coefficients

Dconst = rep(1,nrow(D));
Xhat = cbind(as.matrix(D),Dconst)%*%as.matrix(sPPCoef);

# Display the first 25 reconstructed spatial gene expression patterns
imageBatchDisplay(Xhat,template = template, nrow = 5, ncol = 5, font = 1)

# Evaluate image reconstruction quality
corr = rep(0,length= ncol(X));
for (i in 1:ncol(X)){
    corr[i] = cor(X[,i],Xhat[,i]);
}

# 0 divided by 0 is acceptable; some genes have faint or no expression.

meanCorr = signif(mean(corr, na.rm = TRUE),2);
mainTitle = paste("Histogram of correlation. Mean = ",meanCorr,sep = '');
hist(corr, main = mainTitle );
abline(v = meanCorr, col = 'red');



# 4. Use sPP coefficients to categorize gene expression patterns

# Suppose we want to find gene images that are expressed in PP6
pp6Coef = as.numeric(sPPCoef[6,]);

# Sort the sPP coefficient for PP6 in descending order
sortRes = sort.int(pp6Coef,index.return = TRUE,decreasing = TRUE);
idxTemp = sortRes$ix;

# Display 25 images with strongest expression in PP6
imageBatchDisplay(X[,idxTemp],template = template, imgNames = geneNames[idxTemp],
				       nrow = 5, ncol = 5, font = 1)


 