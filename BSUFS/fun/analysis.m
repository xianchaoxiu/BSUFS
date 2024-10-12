clc; close all; warning off; addpath(genpath(pwd));

load('twomoonData_SPCAFS.mat')
clustersNum = Y(end);
y = Y(:,1);
featuresNum = length(y);  

T=(X-mean(X))';  % 规范后的转置
beta=[100,1e4];   % > 1 PSNP不收敛
[dim,num] = size(T);

m= 8 ;%Y(size(Y,1))-1;   %最后一个元素-1 
avg_times = 20;
m_range = 10:10:100;
cluster_times = length(beta)*avg_times*length(m_range);  %聚类次数
class=zeros(num,cluster_times);

%setting
t = 1.618;      % 1.618
H = eye(num)-(1/num)*ones(num); 
St = T*H*T';
Lambda = zeros(dim,m);
W = ones(dim,m);
V = zeros(dim,m);
lambda = 1;


t0 = tic;
for j=1:length(beta)
   [W,index] = ADMM(W,V,Lambda,beta(1),lambda,St,t,dim,m);
   tend = toc(t0)
   for f = 1:length(m_range)
        I = index(1:5);
        D3 = X(:,I);
        for z=1:20
            [id, I]=kmeans(D3, clustersNum,'MaxIter',200);  
            class(:,length(m_range)*avg_times*(j-1)+avg_times*(f-1)+z)=id;
        end
        disp(['当前选取特征值为ff=',num2str(m_range(f))])
        disp(['当前正则系数gamma=',num2str(beta(j))])
    end
end

result = zeros(2,length(beta)*length(m_range));
index = 1;
sumNMI = 0;
sumACC = 0;
i =1;
for j = 1:cluster_times
    cluster = class(:,j);  %聚类后的数据
    matchedCluster = match(cluster,clustersNum,featuresNum,y);

    ACC = sum(matchedCluster == y) / numel(y);  % 准确率
    NMI = getNMI(matchedCluster, y);
    sumACC = sumACC + ACC;
    sumNMI = sumNMI + NMI;
    if mod(j,avg_times) == 0
        str = strcat('finished :gamma:', num2str(beta(index)), '; m=', num2str(m_range(i)));
        disp(str);
        avgACC = sumACC / avg_times;
        avgNMI = sumNMI / avg_times;
        sumACC = 0;
        sumNMI = 0;
        disp(avgACC);
        disp(avgNMI);
        result(1,j/avg_times ) = avgACC;
        result(2,j/avg_times ) = avgNMI;
        i = i + 1;
    end
    if mod(j, length(m_range)*avg_times) == 0
        index = j/(length(m_range)*avg_times) + 1 ;
        i = 1;
    end
end

t = toc(t0)
%保存数据
save('result.mat','result');
filename = 'admm.xlsx';
xlswrite(filename, result, 'Sheet1');

