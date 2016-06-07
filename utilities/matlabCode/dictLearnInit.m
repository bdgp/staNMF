function D = dictLearnInit(X,K)
% X: the data matrix
% K: the number of dictinary elements

ind = randsample(1:size(X,2),K);
D = X(:,ind);

    
        
        
