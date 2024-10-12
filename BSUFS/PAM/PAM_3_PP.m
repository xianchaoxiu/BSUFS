function [V,iter] = PAM_3_PP(W,X,m,lambda1,lambda2,beta1,beta2,options)

rho1 = 0.001;
rho2 = 0.001;
rho3 = 0.001;

num = size(X,2);
dim = size(X,1); 
H = eye(num)-(1/num)*ones(num); 
St = X*H*X';

%manopt迭代迭代选择的是100   框架迭代是50
%W = ones(size(W));
oldW = zeros(size(W));
oldU = zeros(size(W));
oldV = zeros(size(W));
oldf = 0 ;
tol = 100;
iter = 1; 

maxiter = 100;
OBJ = zeros(maxiter,1);
CHECH =  zeros(maxiter,5);
tic; % 开始计时

while tol > 1e-4 && iter < maxiter
    W = Stiefel_PAM3(W,oldU,oldV,St,dim,m,beta1,beta2,rho1,oldW);

    %元素稀疏
    Y = (beta1/(beta1 + rho2)) * W + (rho2/(beta1 + rho2)) * oldU;  
    T =   lambda1/(beta1+rho2);    
    
    U = fun(Y,T,'PP',options.Epara); %Hard Logarithm
%     disp( max((max(Y))) );
%     disp(T);

    Z = (beta2/(beta2 + rho3)) * W + (rho3/(beta2 + rho3)) * oldV;  
    V =  solve(Z, lambda2/(beta2 + rho3) ,options.sub2,options.para); %sqrt( 2*lambda2/(beta2 + rho3) ) 

        % 定义函数 G1 和 G2 (假设为平方函数)
    G1 = @(x) fun(x,T,'PP',options.Epara);
    G2 = @(x) fun(x, lambda2/(beta2 + rho3) ,'PP',options.para);
    part1 = -trace(W' * St * W); % -Tr(W' * S * W)
    part2 = lambda1 * sum(sum(arrayfun(G1, U))); % lambda1 * sum(sum(G1(w_ij)))
    part3 = lambda2 * sum(arrayfun(G2, vecnorm(V, 2, 2))); % lambda2 * sum(G2(||w_i||_2))
    f = part1 + part2 + part3;
    OBJ(iter) = f;
    tol = abs(f - oldf)/max(1,abs(oldf))
    oldf = f;
    
    tolW = norm(W - oldW,'fro')/(max(norm(W,'fro'),1));
    tolU = norm(U - oldU,'fro')/(max(norm(U,'fro'),1));
    tolV = norm(V - oldV,'fro')/(max(norm(V,'fro'),1));
    A = max([tolW, tolV, tolU]);
%     tol = max([tolW, tolV, tolU])
%     tol = A;

    CHECH(iter,1) = tolW;
    CHECH(iter,2) = tolU;
    CHECH(iter,3) = tolV;
    CHECH(iter,4) = A;
    CHECH(iter,5) = f;

    oldW = W;
    oldV = V;
    oldU = U;
    iter = iter + 1;

end

disp("iter:"+num2str(iter));

end
