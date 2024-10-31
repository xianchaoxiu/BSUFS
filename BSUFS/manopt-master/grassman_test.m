clc; close all; warning off; addpath(genpath(pwd));

% load('result_SPCAFS.mat' ,'result_MINST');
% W =  cell2mat(result_MINST.Ws(1,1));
% V =  cell2mat(result_MINST.Ws(1,1));

load('best_para_SPCA.mat');

W = results.USPS_W;
V = results.USPS_W;





load('USPS.mat');
X = X';
num = size(X,2); 
dim = size(X,1);
lambda = 1e4 ;    % gamma=[1e-6,1e-4,1e-2,1e0,1e2,1e4,1e6]; 
beta = 0.01;
clustersNum = max(Y);
m = clustersNum;

t = 1.1;      % 1.618

t0 = tic;

H = eye(num)-(1/num)*ones(num); 
St = X*H*X';


feature('memstats')
vars = whos;
total_memory = sum([vars.bytes]);
% 显示内存信息
disp(['使用内存：' num2str(total_memory/1024/1024) ' MB']);

Lambda = ones(dim,m);

    
    
    Stiefel = stiefelfactory(dim,m);
    problem.M = Stiefel;
    
    problem.cost  = @(U)robustpca_cost(U,V,St,beta,Lambda);
    problem.egrad = @(U)robustpca_gradient(U,V,St,beta,Lambda);

    [W, xcost, info, options] = trustregions(problem,W);

    checkgradient(problem);



    figure;
    semilogy([info.iter], [info.gradnorm], '.-');
    xlabel('Iteration number');
    ylabel('Norm of the gradient of f');


%     % Prepare a Manopt problem structure for optimization of the given
%     % cost (defined below) over the Grassmann manifold.
%     [p, n] = size(X);
%     manifold = grassmannfactory(p, m);
%     problem.M = manifold;
%     problem.cost = @robustpca_cost;
%     problem.egrad = @robustpca_gradient;
% 	
% 	% Do classical PCA for the initial guess.
% 	% This is just one idea: it is not necessarily useful or ideal.
%     % Using a random initial guess, and starting over for a few different
%     % ones is probably much better. For this example, we keep it simple.
%     [U, ~, ~] = svds(X, d);
% 
%     
% 	% Iteratively reduce the smoothing constant epsilon and optimize
% 	% the cost function over Grassmann.
%     epsilon = 1;
% 	n_iterations = 6;
% 	reduction = .5;
% 	options.verbosity = 2; % Change this number for more or less output
%     warning('off', 'manopt:getHessian:approx');
%     
%     % An alternative way to compute the egrad is to use automatic
%     % differentiation provided in the deep learning toolbox (slower).
%     % Call manoptAD to automatically obtain the egrad 
%     % problem = manoptAD(problem,'egrad');
%     
%     for iter = 1 : n_iterations
%         U = trustregions(problem, U, options);
%         epsilon = epsilon * reduction;
%     end
%     warning('on', 'manopt:getHessian:approx');
%     
%     
% 	% Return the cost as the actual sum of distances, not smoothed.
% 	epsilon = 0;
% 	cost = robustpca_cost(U,V,St,beta,Lambda);
%     
%     
%     
%     % Smoothed cost
%     function value = robustpca_cost(U,V,St,beta,Lambda)
% 
% 
%         value = -trace(U'*St*U) + 0.5*beta*(norm((U - V - Lambda/beta),'fro')^2);
% 
%     end
% 
%     % Euclidean gradient of the smoothed cost (it will be transformed into
%     % the Riemannian gradient automatically by Manopt).
%     function G = robustpca_gradient(U,V,St,beta,Lambda)
%         G = -2*St*W  + beta*(U - V - Lambda/beta);
%     end
% 

    function value = robustpca_cost(U,V,St,beta,Lambda)


        value = -trace(U'*St*U) + 0.5*beta*(norm((U - V - Lambda/beta),'fro')^2);

    end
        function G = robustpca_gradient(U,V,St,beta,Lambda) 
        G = -2*St*U  + beta*(U - V - Lambda/beta);
    end
