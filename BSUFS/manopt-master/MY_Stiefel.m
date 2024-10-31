function [W] = MY_Stiefel(W,V,Lambda,St,dim,m,beta)



Stiefel = stiefelfactory(dim,m);
problem.M = Stiefel;

problem.cost  = @(U)pca_cost(U,V,St,beta,Lambda);
problem.egrad = @(U)pca_gradient(U,V,St,beta,Lambda);

W = trustregions(problem,W);

function value = pca_cost(U,V,St,beta,Lambda)

   value = -trace(U'*St*U) + 0.5*beta*(norm((U - V + Lambda/beta),'fro')^2);

end

function G = pca_gradient(U,V,St,beta,Lambda) 
    G = -2*St*U  + beta*(U - V + Lambda/beta);
end


end