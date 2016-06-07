%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following matlab code produces NMF dictionaries of different
% sizes K. For each K, 100 dictionaries are generated with different
% initial values. Results are stored in ./staNMFDicts/K=* as .mat 
% files (the files, rep*Dict.mat, are already there to save your 
% computation time.
% Please modify the save path to a different directory for your own 
% experiments).
% Last update: March 30, 2016 by Siqi Wu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NOTE: the following code requires the package SPAMS to be installed. 
% Please go to http://spams-devel.gforge.inria.fr/downloads.html
% for download and documentation of this package.

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


loadPath = '../embryoTemplate.csv';
template = csvread(loadPath,1,1);
ind = find(template==1);

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


numPatterns = 15:30; % dictionary sizes
lambda = 0;         % sparsity control for the coefficients alpha
gamma1 = 0; % sparsity control on the dictionary patterns


for k = 1:length(numPatterns)
    
    K =numPatterns(k);
    D0 = dictLearnInit(X,K); % randomly select K columns
                             % from the data matrix X
    
    path = ['./staNMFDicts/K=',num2str(K),'/'];
    mkdir(path);        
    
    param.mode = 2; 
    param.K=K;  
    param.lambda=lambda; 
    param.numThreads=-1; 
    param.batchsize=min(1024,n);
    param.posD = true;   % positive dictionary
    param.iter = 500;  % number of iteration 
    param.modeD = 0;
    param.verbose = 0; % print out update information?
    param.posAlpha = 1; % positive coefficients
    param.gamma1 = gamma1; % penalizing parameter on the dictionary patterns    
    param.D = D0; % set initial values
    
        
    Dtemplate = mexTrainDL(X,param); % learn a dictionary first to
                                     % be used as template for alignment
    
    % For each fixed dictionary K, repeat NMF for 100 times, each with a different initial value
   
    for i = 1:100
        D0 = dictLearnInit(X,K);        
        param.D = D0; 
        tic
        D = mexTrainDL(X,param); 
        toc
        % save the dictionary for future use
        save([path,'rep',num2str(i),'Dict.mat'],'D');
               
    end
end
