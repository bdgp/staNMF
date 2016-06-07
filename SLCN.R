###############################################
# Spatially Local Correlation Networks (SLCN)
# Last update: March 28, 2016 by Siqi Wu
###############################################

source('loadFunctions.R')

# Import the spatial gene expression data matrix
X = read.table(file = '../expressionPatterns.csv',header = TRUE,check.names = FALSE,sep=',');
geneNames = colnames(X);

# Import the sparse PP coefficients for all 1640 images
sPPCoef = read.table(file = '../sPPcoefPatterns.csv',header = TRUE,check.names = FALSE,sep = ',');

# Import the sparse PP coefficients for all 701 genes
sPPCoefGenes = read.table(file = '../sPPcoefGenes.csv',header = TRUE,check.names = FALSE,sep = ',');
ppNames = paste('PP',1:21,sep ='')

rownames(sPPCoefGenes) = ppNames;

# Import the learned dictionary with 21 principal patterns (PP)
D = read.table(file = '../PP.csv',header = TRUE,check.names = FALSE,sep = ',');

# Which of the 701 genes encode transcription factors (TF)?
temp = read.table(file = '../TF.csv',header = TRUE,check.names = FALSE,sep=',');
tfInd = which(temp==1);
tfNamesUniq = unique(geneNames[tfInd]);


# Anterior-posterior PP for SLCN construction
strIdx = c(4,6,7,8,9,17,20,21);

# The names of the six gap genes
gapGenes = c('kni','Kr','hkb','tll','gt','hb');

thrs = .1;

# Linearly ordered PP representation for the six gap genes
linearPPGap = t(sPPCoefGenes[strIdx[2:7],gapGenes]>thrs);


# Which genes are expressed in the above anterior-posterior PP?
strImgIdx = NULL;
for (i in 1:length(strIdx)){
       sPPCoefTemp = sPPCoef[strIdx[i],];
       idxTemp0 = which(sPPCoefTemp>thrs);
       strImgIdx = union(strImgIdx,idxTemp0);
}
strImgIdx = intersect(tfInd,strImgIdx)
nameStrUniq = unique(geneNames[strImgIdx]);


# Construct the six SLCN

corMatTemp = matrix(0,nrow = length(nameStrUniq),ncol = length(nameStrUniq));
rownames(corMatTemp) = nameStrUniq;
colnames(corMatTemp) = nameStrUniq;


corVecList = nameList = edgeList = corrList = list();
Thrs = matrix(0,nrow = length(strIdx) - 2,ncol =2);

for (i in 2:(length(strIdx)-1)){

    idxTemp = NULL;
    if (i>1){
       sPPCoefTemp = sPPCoef[strIdx[i-1],];
       idxTemp0 = which(sPPCoefTemp>thrs);
       idxTemp = union(idxTemp,idxTemp0);
    }

    sPPCoefTemp = sPPCoef[strIdx[i],];
    idxTemp1 = which(sPPCoefTemp>thrs);
    idxTemp = union(idxTemp,idxTemp1);

    if (i < length(strIdx)){
       sPPCoefTemp = sPPCoef[strIdx[i+1],];
       idxTemp2 = which(sPPCoefTemp>thrs);
       idxTemp = union(idxTemp,idxTemp2);
    }

    idxTemp = intersect(tfInd,idxTemp)

    Dmask = D[,strIdx[i]];
    
    Xtemp = X[,idxTemp];
    sPPCoefTemp = sPPCoef[,idxTemp];
    corTemp = matrix(1,nrow = ncol(Xtemp),ncol = ncol(Xtemp));

    for (k in 1:(ncol(Xtemp)-1)){
    	for (j in (k+1):ncol(Xtemp)){
	    corTemp[k,j] = corrWeighted(Xtemp[,k],Xtemp[,j],Dmask,method = 'pearson')$corrXY;	    
	    corTemp[j,k] = corTemp[k,j];
	 }
    }

    nameTemp = geneNames[idxTemp];
    nameUniq = unique(nameTemp);
    nameList[[i]] = nameUniq;
    corUniq = matrix(1,nrow = length(nameUniq),ncol = length(nameUniq));

    for (k in 1:(length(nameUniq)-1)){
    	for (j in (k+1):length(nameUniq)){
	    idxTemp1 = which(nameTemp == nameUniq[k]);
	    idxTemp2 = which(nameTemp == nameUniq[j]);
	    C = corTemp[idxTemp1,idxTemp2];
	    idxTemp3 = min(which(abs(C) == max(abs(C))))
	    corUniq[k,j] = C[idxTemp3];
	    corUniq[j,k] = corUniq[k,j];
	 }
    }

    colnames(corUniq) = nameUniq;
    rownames(corUniq) = nameUniq;
    
    corTemp = corUniq;

    upperIdx = which(upper.tri(corTemp));
    corVec = corTemp[upperIdx];
    
    corVecList[[i]] = corVec;    
    quan = quantile(corVec,c(0.05,0.95));

    Thrs[i-1,] = quan;


    corrList[[i]] = corMatTemp;
    corrList[[i]][nameUniq,nameUniq] = corUniq;

    edgeMatTemp = corUniq;
    edgeMatTemp[corUniq > quan[1] & corUniq < quan[2]] = 0;
    edgeList[[i]] = corMatTemp;
    edgeList[[i]][nameUniq,nameUniq] = edgeMatTemp;

    
    cat('PP',strIdx[i],' gap gene SLCN:\n', sep='');
    
    for (k in 1:(length(gapGenes)-1)){
    	for (j in (k+1):length(gapGenes)){
    	    cTemp = corTemp[which(nameUniq == gapGenes[k]), which(nameUniq==gapGenes[j])];

	    if (length(cTemp)){
	       if (max(cTemp) > quan[2]){
	       	  upperTemp = signif(sum(corVec>=max(cTemp))/length(corVec),2);
	       	  cat(paste(gapGenes[k],'<->',gapGenes[j],'\n',sep = ''));
	       }
	       if (min(cTemp) < quan[1]){
	       	  lowerTemp = signif(sum(corVec<=min(cTemp))/length(corVec),2);
	       	  cat(paste(gapGenes[k],'|-|',gapGenes[j],'\n',sep = ''));
	       }
	    }
	 }
    }

    cat('\n');
}


# Historgrams for the local correlations
par(mfrow=c(2,3));
for (i in 1:6){
    hist(corVecList[[i+1]],breaks=20, main = ppNames[strIdx[i+1]],
         xlab = 'local correlation',col='blue');
    quan = quantile(corVecList[[i+1]],c(0.05,0.95));
    abline(v = quan,col = 'red');
}



# Plotting the six TF SLCN  
# Save as SLCN1-3.pdf

numStrTF = length(nameStrUniq);
theta = seq(from = 360, to = 0, length.out = numStrTF);
lowerIdx = floor(numStrTF/4);
upperIdx = ceiling(numStrTF/4*3);
theta[lowerIdx:upperIdx] = 180 + theta[lowerIdx:upperIdx];
coord = genCircle(numStrTF+1);
coord$x = coord$x[1:numStrTF];
coord$y = coord$y[1:numStrTF];
idxGap = NULL;
for (i in 1:length(gapGenes)){
    idxGap = c(idxGap,which(nameStrUniq==gapGenes[i]));
}
for (ii in 1:3){
    savePath = paste('./SLCN', ii,'.pdf',sep='');
    pdf(savePath,width = 20, height = 10);
    par(mfrow=c(1,2));
    idxTemp0 = 2:3+2*(ii-1);
    for (i in idxTemp0){
    
	par(mar = c(3,3,3,3));
    	edgeMatTemp = edgeList[[i]];
	distX = as.dist(1 - edgeMatTemp);
    	hr <- hclust(distX, method="average"); # h-clustering
    	idx0 = hr$order;
    	nameStrUniq1 = nameStrUniq[idx0];
    	edgeMatTemp = edgeMatTemp[idx0,idx0];
    	
	nb = colSums(abs(edgeMatTemp)>=0.001);
    	nb[nb == 1] = 0;
    	cex = 3*sqrt(nb)/sqrt(max(nb));
    	mainTitle = paste('PP',strIdx[i],sep='');
    	plot(coord$x[-idxGap],coord$y[-idxGap],pch = 16,cex = cex[-idxGap],main = mainTitle,
         xaxt = 'n', yaxt = 'n',xlim = c(-1.1,1.1),ylim = c(-1.1,1.1),
         xlab = '',ylab = '', col = '#FF000080',bty = 'n');
    	 points(coord$x[idxGap],coord$y[idxGap],pch = 16,cex = cex[idxGap], col = '#FF000080');
    	 for (k in 1:(numStrTF-1)){
             for (j in (k+1):numStrTF){
             	 if (edgeMatTemp[k,j]<=-0.001){
                    segments(coord$x[k],coord$y[k],coord$x[j],coord$y[j],col = '#FF000099',lwd = 1);
            	   }
            	 if (edgeMatTemp[k,j]>=0.001){
                    segments(coord$x[k],coord$y[k],coord$x[j],coord$y[j],col = '#0000FF99',lwd = 1);
            	  }
             }
    	 }
    	for (k in 1:length(nameStrUniq)){
            if (sum(gapGenes == nameStrUniq1[k])){
               text(coord$x[k]*1.09,coord$y[k]*1.09,nameStrUniq1[k],srt = theta[k],cex = 1,col ='red');
            }else{
		text(coord$x[k]*1.09,coord$y[k]*1.09,nameStrUniq1[k],srt = theta[k],cex = 1);
            }

        }
    }
    dev.off();
}


# Plotting SLCN sub-network for the six gap genes
gapGenes = c('gt','hb','hkb','kni','Kr','tll');
numGap = length(gapGenes)
coord = genCircle(numGap);


par(mfrow=c(2,3));
for (i in 2:7){

    edgeMatTemp = edgeList[[i]][gapGenes,gapGenes];
    par(mar = c(1,1,1,1));
    plot(coord$x,coord$y,pch = 16,cex = 7,bty = 'n',
    xaxt = 'n', yaxt = 'n',xlim = c(-1.4,1.4),ylim = c(-1.4,1.4),
    xlab = '', main = ppNames[strIdx[i]],ylab = '', col = 'black');

   for (k in 1:(numGap-1)){
       for (j in (k+1):numGap){
    	   if (edgeMatTemp[k,j]<=-0.001){
    	      segments(coord$x[k],coord$y[k],coord$x[j],coord$y[j],col = '#FF0000',lwd = 4);
              }
 	   if (edgeMatTemp[k,j]>0.001){
    	      segments(coord$x[k],coord$y[k],coord$x[j],coord$y[j],col = '#0000FF',lwd = 4);
            }
   	}
    }
    points(coord$x,coord$y,pch = 16,cex = 7,col = 'black');
    points(coord$x,coord$y,pch = 16,cex = 6.5,col = 'white');
    text(coord$x,coord$y,gapGenes,cex = 1.5,col = 'black',font = 3)
}



# Global correlation analysis: correlating genes using whole embryos
Xtemp = X[,strImgIdx];
corGlobal = cor(Xtemp,method = 'pearson');
corGlobalUniq = matrix(1,nrow = length(nameStrUniq),ncol = length(nameStrUniq));
nameTemp = geneNames[strImgIdx]
for (k in 1:(length(nameStrUniq)-1)){
    for (j in (k+1):length(nameStrUniq)){
	    idxTemp1 = which(nameTemp == nameStrUniq[k]);
	    idxTemp2 = which(nameTemp == nameStrUniq[j]);
            C = corGlobal[idxTemp1,idxTemp2];
	    idxTemp3 = which(abs(C) == max(abs(C)))
            corGlobalUniq[k,j] = C[idxTemp3];
	    corGlobalUniq[j,k] = corGlobalUniq[k,j];
        }
}

colnames(corGlobalUniq) = nameStrUniq;
rownames(corGlobalUniq) = nameStrUniq;
upperIdx = which(upper.tri(corGlobalUniq));
corGlobalVec = corGlobalUniq[upperIdx];

# Count the total number of positive/negative edges for the local networks
totalEdgeCoef = edgeList[[2]]*0;
for (i in 2:(length(strIdx)-1)){
    temp = edgeList[[i]];
    totalEdgeCoef = totalEdgeCoef + abs(temp);
}
upperIdx = which(upper.tri(temp));
 
totalNumEdges = sum(totalEdgeCoef[upperIdx]>0);
totalNumPosEdges = ceiling(totalNumEdges/2)
totalNumNegEdges = ceiling(totalNumEdges/2)

# Sort the global correlation vector
sortCorTemp = sort(corGlobalVec);
lower = sortCorTemp[totalNumPosEdges];
upper = sortCorTemp[length(corGlobalVec)-totalNumPosEdges+1];

edgeGlobal = corGlobalUniq;
edgeGlobal[corGlobalUniq > lower & corGlobalUniq < upper] = 0;


# Gap gene network: global correlation
gapGenes = c('gt','hb','hkb','kni','Kr','tll');
numGap = length(gapGenes)
coord = genCircle(numGap);
edgeMatTemp = edgeGlobal[gapGenes,gapGenes]

X11();
plot(coord$x,coord$y,pch = 16,cex = 7,bty = 'n',
     xaxt = 'n', yaxt = 'n',xlim = c(-1.4,1.4),ylim = c(-1.4,1.4),
     xlab = '', main = 'Global correlation',ylab = '', col = 'black');

for (k in 1:(numGap-1)){
    for (j in (k+1):numGap){
        if (edgeMatTemp[k,j]<=-0.001){
            segments(coord$x[k],coord$y[k],coord$x[j],coord$y[j],col = '#FF0000',lwd = 4);
        }
        if (edgeMatTemp[k,j]>0.001){
            segments(coord$x[k],coord$y[k],coord$x[j],coord$y[j],col = '#0000FF',lwd = 4);
        }
    }
}
points(coord$x,coord$y,pch = 16,cex = 7,col = 'black');
points(coord$x,coord$y,pch = 16,cex = 6.5,col = 'white');

text(coord$x,coord$y,gapGenes,cex = 1.5,col = 'black',font = 3)

