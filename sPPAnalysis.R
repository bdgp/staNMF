###########################################################################
# sPP coefficients: heatmap visualization and pairwise interaction between
# anterior-posterior PP.
# Last update: April 5, 2016 by Siqi Wu
###########################################################################

# Import sparse PP coefficients for 701 genes
sPPCoefGenes = read.csv(file = '../sPPcoefGenes.csv',header = TRUE,check.names = FALSE);
geneNames = colnames(sPPCoefGenes);
K = 21
ppNames = paste('PP',1:K,sep ='')

# Import the learned dictionary with 21 principal patterns (PP)
D = read.csv(file = '../PP.csv',header = TRUE,check.names = FALSE);

# Import the elliptical template for spatial gene expression visualization
template = read.csv(file = './utilities/embryoTemplate.csv',sep = ',',header = TRUE,row.names = 1);


source('loadFunctions.R');

# Which of the 701 genes are Computed Genes (CG)?
geneNamesUniq = unique(geneNames)
characterized = rep(1,length(geneNamesUniq));
for (i in 1:length(geneNamesUniq)){
    temp = grep('CG',geneNamesUniq[i]);
    if (length(temp) > 0){
       characterized[i] = 0;
    }
}
cgIdx = which(characterized == 0)


# Cluster genes according to their highest coefficient PP
thrs = 0.1
labels = rep(0,length = ncol(sPPCoefGenes));
for (i in 1:ncol(sPPCoefGenes)){
    sPPCoefTemp = sPPCoefGenes[,i];
    if (max(sPPCoefTemp)>thrs){
       labels[i] = min(which(sPPCoefTemp == max(sPPCoefTemp)));
    }
}
names(labels) = colnames(sPPCoefGenes);


# Reorder the genes in each cluster by hierarchical clustering
idxReorder = rep(0,length = length(which(labels>0)));
i = 1;
for (k in 1:K){
    idxTemp = which(labels == k);
    corrTemp = cor(sPPCoefGenes[,idxTemp]);
    distTemp = 1 - corrTemp;
    hcTemp = hclust(as.dist(distTemp));
    idxTemp = idxTemp[hcTemp$order];
    lenTemp = length(idxTemp);
    start = i;
    idxTemp2 = start:(start+lenTemp-1);
    idxReorder[idxTemp2] = idxTemp;
    i = i + lenTemp;
}


# Plotting the heatmap
imagesc(sPPCoefGenes[,idxReorder],colorScale='blue',margin = c(6,6,2,2),printLabels = F);
        
temp1 = table(labels[idxReorder]);
temp = cumsum(temp1);
abline(v=(temp[1:20]-0.5)/(length(idxReorder)-1),col = 'red',lty =2,lwd = 2);
axis(1, at=(temp - temp1/2)/(length(idxReorder)-1),labels=ppNames, las = 1);
temp2 = (0:(K-1))/(K-1);
axis(2, at=1-temp2, labels=ppNames, las = 1);



# Number of genes expressed in each PP
ppCGIdx = ppGenesIdx = list();
numCGPP = numGenePP = rep(0,K);
# numGenePP: total number of genes expressed in each PP
# numCGPP: number of computed genes (CG) expressed in each PP

# Threshold for the sPP coefficients
thrs = .1;
for (i in 1:K){
    ppGenesIdx[[i]] = which(sPPCoefGenes[i,]>thrs);
    ppCGIdx[[i]] = which(sPPCoefGenes[i,geneNamesUniq[cgIdx]]>thrs);
    numGenePP[i] = length(ppGenesIdx[[i]]);
    numCGPP[i] = length(ppCGIdx[[i]]);
}


# Number of genes expressed in a PP pair
countCommonGenes = matrix(0,nrow=K,ncol=K);
fracCommonGenes = countCommonGenes;
q = 1;
for (i in 1:(K-1)){
    for (j in (i+1):K){
        intTemp = intersect(ppGenesIdx[[i]],ppGenesIdx[[j]]);
        intTempUnion = union(ppGenesIdx[[i]],ppGenesIdx[[j]]);
        countCommonGenes[i,j] = length(intTemp);
        countCommonGenes[j,i] = countCommonGenes[i,j];
	# Jaccard distance of two sets
        fracCommonGenes[i,j] = length(intTemp)/length(intTempUnion);
        fracCommonGenes[j,i] = fracCommonGenes[i,j];
    }
}


# Fraction of shared PP vs PP distance
patNameMat = matrix('',nrow=K,ncol=K);
for (i in 1:(K-1)){
    for (j in (i+1):K){
        patNameMat[i,j] = paste(i,'-',j,sep='');
        patNameMat[j,i] = paste(i,'-',j,sep='');
    }
}

# Obtain the centroid for each PP
L = patternCentroid(as.matrix(D),template=template);

# Compute the pairwise distances between PP centroids
distMatPat = as.matrix(dist(L))

patIdxTemp = setdiff(1:K,10:16);

corrPP = t(fracCommonGenes);
corrTemp = t(corrPP[patIdxTemp,patIdxTemp]);
distMatTemp = distMatPat[patIdxTemp,patIdxTemp];
patNameTemp = patNameMat[patIdxTemp,patIdxTemp];
upperIdx = which(upper.tri(distMatTemp));
distMatVec = distMatTemp[upperIdx];
corrVec = corrTemp[upperIdx];
patNameVec = patNameTemp[upperIdx];

plot(distMatVec,corrVec,type='n',xlab = 'PP distance', ylab = 'Fraction of common genes');
text(cbind(distMatVec,corrVec),patNameVec)






