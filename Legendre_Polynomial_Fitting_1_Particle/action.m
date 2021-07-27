function [S_int] = action(q0, Legendre, Legendre_dif, C, r, a, div, z, N, dQ)
    dz = z(2) - z(1);
    tra = q0;
    N = div;
    
    der = zeros(1,N);
    der(1) = (tra(2)-tra(1))/dz;

    for i = 2:N-1
        %der(i) = ((tra(i+1) - tra(i-1))/(2*dz));
        der(i) = ((tra(i+1) - tra(i))/(dz));
    end
    der(N) = ((tra(N)-tra(N-1))/(dz));
    
    func = zeros(1,N);
    for k=1:N
        pre = (1 - (z(k)^2))/r;
        func(k) = pre*(1/2) * der(k)^2 + (1/pre)*(1/4)*(tra(k)^2 - a)^2;
    end
    
    %trapezoid rule
    actionn = 0;
    for L = 2:N
        actionn = actionn + (dz * (func(L-1)+func(L))/2);   
    end
    S_int = actionn;
end