function tra = newstep(tra,N,sigma,a)
    R = randi((N/2)-1) + 1;
    D = (normrnd(tra(R),sigma));
    temp = tra(R);
    tra(R) = D;
    while tra(R) > 0 || tra(R) < -sqrt(a)
       D = (normrnd(temp,sigma));
       tra(R) = D;
    end
    tra(N-(R-1)) = -tra(R);

end