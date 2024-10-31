function [W,fval] = MY_GRASSMAN(W,V,Lambda,St,dim,m,beta)


options =struct();
options.debug = 0;

options.verbosity = 1;  % 0为无输出，越大越详细

manifold = grassmannfactory(dim, m);
problem.M = manifold;

problem.cost  = @(U)robustpca_cost(U,V,St,beta,Lambda);
problem.egrad = @(U)robustpca_gradient(U,V,St,beta,Lambda);

[W, fval, info, options] = trustregions(problem,W,options);

function value = robustpca_cost(U,V,St,beta,Lambda)

   value = -trace(U'*St*U) + 0.5*beta*(norm((U - V + Lambda/beta),'fro')^2);

end
function G = robustpca_gradient(U,V,St,beta,Lambda) 
    G = -2*St*U  + beta*(U - V + Lambda/beta);
end

end