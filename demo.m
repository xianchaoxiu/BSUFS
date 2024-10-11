%%  paras setting
load("COIL20.mat");

current_lambda1 = 1000000; 
current_lambda2 = 10000;
current_beta1 = 100;
current_beta2 = 1000000;
current_para1 =  2/3;
current_para2 = 2/3;
  
avg_times = 50;        
cureent_features_select = 100;

%%  initialize
X =  (X-mean(X));
T = X';
[num,dim] = size(X);
clustersNum = max(Y);
y = Y(:,1);
samplesNum = length(y); 
m = clustersNum;

[W,~]  = qr(ones(dim,m),0);

options.sub2 = 'PP';
options.para = current_para1;
options.Epara = current_para2;

    
%% feature seleciton algorithms, find the best V
[V,iter] = PAM_3_PP(W,T,m,current_lambda1,current_lambda2,current_beta1,current_beta2,options);

%% select top features, and analyze
result = zeros(2,1);
result_sum = zeros(2,1);  
for time = 1:avg_times
    cluster = CLUSTER(V, X, clustersNum,cureent_features_select);
    result = ANALY(cluster, y, clustersNum, samplesNum);
    result_sum = result_sum + result;
end
acc = result_sum(1,1)/avg_times
nmi = result_sum(2,1)/avg_times
result(1,1) = acc;
result(2,1) = nmi;