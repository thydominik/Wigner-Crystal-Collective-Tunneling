function [shift] = f_initshift(eta, a, pos, NoP)  
    Q = 0;
    for i = 1:NoP
        F(i) = 0.25 * (pos(i,1)^2 - a)^2;
        for j = (i + 1):NoP
            Q = Q + eta/(abs(pos(i,1) - pos(j,1)));
        end
    end
    shift =  Q + sum(F);

end

