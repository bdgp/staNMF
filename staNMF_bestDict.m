%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the optimal rank (K = 21 in our case), the following code
% finds the best dictionary that has the lowest NMF objective
% function value among the 100 dictionaries generated through NMF 
% with different initializations.  
% Last update: April 1, 2016 by Siqi Wu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

addpath('./utilities/matlabCode/');

% replace the following two paths with your SPAMS installation paths
addpath('~/matlab/packages/spams-matlab/');
addpath('~/matlab/packages/spams-matlab/build/');


loadPath = '../expressionPatterns.csv';

X = csvread(loadPath,1,0);
fid = fopen(loadPath);
tline = fgetl(fid);

% obtain the gene names for the images
geneNames = textscan(tline,'%s','delimiter',',')
geneNames = geneNames{1};
for i = 1:size(X,2)
    geneNames{i} = geneNames{i}(2:(end-1));
end

width = 32;
height = 16;

[m,n] = size(X);

% weighted NMF to consider replicates for the same gene:
weighted = true;
if weighted
    gnUniq = unique(geneNames);
    geneNum = zeros(1,length(gnUniq));
    weight = zeros(1,length(geneNames));
    for i = 1:length(gnUniq)
        geneNum(i) = length(strmatch(gnUniq(i),geneNames, 'exact'));
    end
    
    for i = 1:length(geneNames)
        idxTemp = strmatch(geneNames(i),gnUniq,'exact');
        weight(i) = 1/geneNum(idxTemp);        
    end
    
    for i = 1:length(geneNames)
        X(:,i) = sqrt(weight(i))*X(:,i);
    end
end

    
K = 21; % K = 21 the NMF rank selected by staNMF
initPath = './staNMFDicts/';
numReplicates = 100;
path = [initPath, 'K=', num2str(K),'/'];                        
R = zeros(numReplicates,1);
lambda = 0;         % sparsity control for the coefficients alpha
gamma1 = 0; % sparsity control on the dictionary patterns
for L = 1:numReplicates
    loadPath = [path,'rep',num2str(L),'Dict.mat']; 
    load(loadPath);
        
    param.mode = 2;     
    param.lambda=lambda; 
    param.pos = 1; % positive coefficients
    
    % standardize each column of the dictionary to have maximum
    % intensity equal to one.
    for k = 1:K
        D(:,k) = D(:,k)/max(D(:,k));
    end
    
    % nonnegative least squares
    alpha = mexLasso(X,D,param);    
    % compute objective function value:
    R(L) = mean(sum((X-D*alpha).^2));   

end


% find the dictionary with the lowest objective function value
L = find(R == min(R));    
loadPath = [path,'rep',num2str(L),'Dict.mat']; 
load(loadPath);

for k = 1:K
    D(:,k) = D(:,k)/max(D(:,k));
end



    
