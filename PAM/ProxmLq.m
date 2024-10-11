function [px,T] = ProxmLq(a, alam, q)
    switch q
        case 0
            t = sqrt(2*alam);
            T = find(abs(a) > t);
            px = zeros(size(a)); % initialize px with zeros
            px(T) = a(T);
            
        case 1/2
            t = (3/2)*alam^(2/3);
            T = find(abs(a) > t);
            aT = a(T);
            phi = acos((alam/4)*(3./abs(aT)).^(3/2));
            px = zeros(size(a)); % initialize px with zeros
            px(T) = (4/3)*aT.*( cos((pi-phi)/3) ).^2;
            
        case 2/3
            t = (128/27)^(1/4)*(alam)^(3/4);
            T = find(abs(a) > t);
            aT = a(T);
            tmp1 = aT.^2/2;
            tmp2 = sqrt(tmp1.^2 - (8*alam/9)^3);
            phi = (tmp1+tmp2).^(1/3)+(tmp1-tmp2).^(1/3);
            px = zeros(size(a)); % initialize px with zeros
            px(T) = sign(aT)/8.*( sqrt(phi)+sqrt(2*abs(aT)./sqrt(phi)-phi) ).^3;
            
        otherwise
            fprintf(' Input ''q'' is incorrect!!! \n');
            fprintf(' Please re-enter ''q'' one value of {0,1/2,2/3}!!!\n')
    end
end