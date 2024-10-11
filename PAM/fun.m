function B = fun(x,lambda,label,para)
    if nargin < 4
        switch label
            case 'Hard'
                para = 0; % Hard_Thresholding doesn't need a parameter
            case 'soft'
                para = 0.5;
            case 'Logarithm'
                para = 1e-1;
            case 'LSP'
                para = 1e-6;
            case 'ETP'
                para = 1e-4;
            case 'LP'
                para = 1e-6;
            case 'LOG'
                para = 1e-1;
            case 'GEM'
                para = 1e-6;
            case 'PP'
                para = 0.5;
            case 'l1'
                para = 1;
            otherwise
                error('未知的阈值函数标签');
        end
    end

    switch label
        case 'Hard'
            B = Hard_Thresholding(x, lambda);
        case 'soft'
            B = Generalized_Soft_Thresholding(x, lambda, para);
        case 'Logarithm'
            B = Logarithm(x, lambda, para);
        case 'LSP'
            B = GAI_LSP(x,lambda,para,100);
        case 'ETP'
            B = GAI_ETP(x,lambda,para, 100);
        case 'LP'
            if(para == 0)
                B = Hard_Thresholding(x, lambda);
            else
                B = Generalized_Soft_Thresholding(x, lambda, para);
            end
        case 'LOG'  
            B = Log(x,lambda,para);
        case 'GEM'
            B = pnn(x,lambda,para);
        case 'PP'
            B = ProxmLq(x,lambda,para);
        case 'l1'
            B = soft(x,lambda);
        otherwise
            error('未知的阈值函数标签');
    end
end