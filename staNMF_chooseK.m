%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code analyzes the stability of dictionaries stored
% in ./staNMFDicts/K=*/ 
% Last update: March 30, 2016 by Siqi Wu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('./utilities/matlabCode/');

% replace the following two paths with your SPAMS installation paths
addpath('~/matlab/packages/spams-matlab/');
addpath('~/matlab/packages/spams-matlab/build/');

numReplicates = 100;
numPatterns = 15:30;

initPath = './staNMFDicts/';
    
computeDictCorrError = true;

% for each K, compute the pairwise distance between
% dictionaries using an Amary-type metric. Save the
% pairwise distance matrix as distMatrixDictCorr.mat
% in ./staNMFDicts/K=*/

if computeDictCorrError
    K = numPatterns(1);
    path = [initPath, 'K=', num2str(K),'/'];                            
    loadPath = [path,'rep1Dict.mat']; 
    load(loadPath);
    d = size(D,1);

    estStability = zeros(1,length(numPatterns));
    
    for k = 1:length(numPatterns)
        K = numPatterns(k)
        path = [initPath, 'K=', num2str(K),'/'];                        
        Dhat = zeros(d,K,numReplicates);
        for L = 1:numReplicates
            loadPath = [path,'rep',num2str(L),'Dict.mat']; 
            load(loadPath);
            Dhat(:,:,L) = D;
        end
    
        distMat = zeros(numReplicates,numReplicates);    
        
        for q = 1:numReplicates
            for p = q:numReplicates
                CORR = corr(Dhat(:,:,q),Dhat(:,:,p));                
                distMat(q,p) = amariMaxError(CORR);
                distMat(p,q) = distMat(q,p);
            end
        end
        save([path,'distMatrixDictCorr.mat'],'distMat');    
    end
end


% load the pairwise dictionary distance matrix and average
% the results.

estStability = zeros(length(numPatterns),1);    
for k = 1:length(numPatterns)
    K = numPatterns(k);
    path = [initPath, 'K=', num2str(K),'/'];      
    load([path,'distMatrixDictCorr.mat']);        
    estStability(k) = sum(distMat(:))/(numReplicates*(numReplicates-1));
end

figure; plot(numPatterns,estStability);
xlabel('K');
ylabel('Instability');

    


