function Dstd = dictStd(D,method)
% D: a dictionary to be standardized
% method: which standardization method? 
% 1: divided by L-1 norm
% 2: divided by L-2 norm
% 0: divided by maximum intensity
[p,K] = size(D);
Dstd = D;
for i = 1:K
    if method == 0
        Dstd(:,i) = D(:,i)/max(D(:,i));
    elseif method == 1
        l1Norm = sum(abs(D(:,i)));
        Dstd(:,i) = D(:,i)/l1Norm;
    elseif method == 2
        l2Norm = norm(D(:,i));
        Dstd(:,i) = D(:,i)/l2Norm;
    end
end

