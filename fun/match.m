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