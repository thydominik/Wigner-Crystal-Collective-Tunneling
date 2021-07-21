function action = actcalc(tra,r,N,z,dz,a,E_p)
%calculate first the derivative that emerges in the action.
    %the second order diff. can couse some oscillation, but will smoothed
    %out by the end low temperatures, the end points of the curve are first
    %order approx.
    der = zeros(1,N);
    der(1) = ((tra(2)-tra(1))/dz);

    for i = 2:N-1
        %der(i) = ((tra(i+1) - tra(i-1))/(2*dz));
        der(i) = ((tra(i) - tra(i-1))/(dz));
    end
    der(N) = ((tra(N)-tra(N-1))/(dz));
    
    func = zeros(1,N);
    for k=1:N
        pre = (1 - (z(k)^2))/r;
        func(k) = pre*(1/2) * der(k)^2 + (1/pre)*(1/4)*(tra(k)^2 - a)^2;    %(((tra(k)^4)/4) - (tra(k)^2)/2);
    end
    
    %trapezoid rule
    action = 0;
    for L = 2:N
        action = action + (dz * (func(L-1)+func(L))/2);   
    end
end