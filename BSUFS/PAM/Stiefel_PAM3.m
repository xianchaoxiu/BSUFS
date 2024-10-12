function [W] = Stiefel_PAM3(W,U,V,St,dim,m,beta1,beta2,rho1,Wk)

options.maxiter = 100;
options.verbosity = 2;
Stiefel = stiefelfactory(dim,m);
problem.M = Stiefel;

problem.cost  = @(W)pca_cost(W,U,V,St,beta1,beta2,rho1,Wk);
problem.egrad = @(W)pca_gradient(W,U,V,St,beta1,beta2,rho1,Wk);
% problem.ehess = @(W)pca_hess(dim,m,V,St,beta1,beta2,rho1,Wk);

W = trustregions(problem,W,options);

function value = pca_cost(W,U,V,St,beta1,beta2,rho1,Wk)

   value = -trace(W' * St * W) + (beta1/2) * norm(W - U, 'fro')^2 + (beta2/2) * norm(W - V, 'fro')^2 + (rho1/2) * norm(W - Wk, 'fro')^2;

end

function G = pca_gradient(W,U,V,St,beta1,beta2,rho1,Wk) 
    G = -2*St*W  + beta1*(W - U) + beta2*(W - V) + rho1*(W - Wk);
end

% function H = pca_hess(dim,m,V,St,beta1,beta2,rho1,Wk) 
%     H = -2*kron(eye(dim),St)  + (beta1 + beta2 + rho1)*eye(dim*m);
% end

end