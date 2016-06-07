function [distM,rowTemp,colTemp] = amariMaxError(A)
% A version of Amari error
% Computed as
%   1/(2K)*(sum_j (1 - max_{ij} A_{ij}) + sum_i (1 - max_{ij} A_{ij}))

[n,m] = size(A);


if n == 1 & m == 1
    distM = double(A==0);
    rowTemp = 0;
    cowTemp = 0;
    return;
end

maxCol = max(abs(A),[],1);
colTemp0 = 1 - maxCol;
colTemp = mean(colTemp0);

maxRow = max(abs(A),[],2);
rowTemp0 = 1 - maxRow;
rowTemp = mean(rowTemp0);

distM = (rowTemp + colTemp)/2;
