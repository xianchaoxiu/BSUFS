function [cluster,index] = CLUSTER(W,X,clustersNum,h)
% This code aims at clustering the W matrix, m features are selected by sorting 
% \|w^i\|^2 in descending order.
%--------------------------------------------------------------------------
% Inputs:
%     ~
% Outputs:
%     cluster: clustering result, a vector   
%--------------------------------------------------------------------------

n = size(X,1);
cluster=zeros(n,1);



sqW = (W.^2);
sumW = sum(sqW,2);
[~,index] = sort(sumW,'descend');

I = index(1:h);
D3 = X(:,I);

id=kmeans(D3, clustersNum,'MaxIter',200);  
cluster(:,1)=id;



