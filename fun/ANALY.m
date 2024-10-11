function [result] = ANALY(cluster,y,clustersNum,featuresNum)
% This code aims at evaluating ACC and NMI of clustering results 
%--------------------------------------------------------------------------
% Inputs:
%     ~
% Outputs:
%     cluster:    result(1,1) is ACC,result(1,2) is NMI
%--------------------------------------------------------------------------
result = zeros(2,1);

matchedCluster = match(cluster,clustersNum,featuresNum,y);

ACC = sum(matchedCluster == y) / numel(y);  % 准确率
NMI = getNMI(matchedCluster, y);
result(1,1) = ACC;
result(2,1) = NMI;


function [matchedCluster] = match(cluster,clustersNum,featuresNum,y)
        confusionMatrix = zeros(clustersNum, clustersNum); % 初始化混淆矩阵
        for i = 1:featuresNum
            confusionMatrix(cluster(i), y(i)) = confusionMatrix(cluster(i), y(i)) + 1;  %行 为cluster的信息
        end
        [match, cost] = munkres(-confusionMatrix); 
        matchedCluster = zeros(featuresNum,1);
        for i = 1:featuresNum
            matchedCluster(i) = match(cluster(i));
        end
end

end