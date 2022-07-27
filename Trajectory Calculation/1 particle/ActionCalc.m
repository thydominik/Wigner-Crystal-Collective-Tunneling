function action = actioncalc(Pos, r , NoP, z, dz, alpha)
%calculate first the derivative that emerges in the action.
    %the second order diff. can couse some oscillation, but will smoothed
    %out by the end low temperatures, the end points of the curve are first
    %order approx.
    der = zeros(1,NoP);
    der(1) = (Pos(2)-Pos(1))/dz;

    for i = 2:NoP-1
        %der(i) = ((tra(i+1) - tra(i-1))/(2*dz));
        der(i) = ((Pos(i+1) - Pos(i))/(dz));
    end
    der(NoP) = ((Pos(NoP)-Pos(NoP-1))/(dz));
    
    func = zeros(1,NoP);
    for k = 2:NoP-1
        pre = (1 - (z(k)^2))/r;
        func(k) = pre*(1/2) * der(k)^2 + (1/pre)*(1/4)*(Pos(k)^2 - alpha)^2;
    end
    
    %trapezoid rule
    action = 0;
    for L = 2:NoP-1
        action = action + (dz * (func(L - 1) + func(L))/2);   
    end
end