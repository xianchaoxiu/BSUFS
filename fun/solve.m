function A = solve(Y, lambda, label,para)
    [d, m] = size(Y);
    normY = vecnorm(Y, 2, 2);
    A = zeros(d, m);
   
    % Logarithm(Y,lambda/beta,1e-1);
    % Generalized_Soft_Thresholding(Y, lambda/beta, 0.5);
    % 添加其他标签和函数的映射
    % functionHandleMap('其他标签') = @其他函数;
    
    for i = 1:d
        if Y(i, :) == 0
            A(i, :) = 0;
        else
            B = fun(normY(i),lambda,label,para);    %这个地方有问题
            A(i, :) = (Y(i, :) / normY(i)) * B;
        end
    end
end
