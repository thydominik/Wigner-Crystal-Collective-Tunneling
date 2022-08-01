function pos = f_newstep(pos, alpha, eta, NoP, Nt, sigma, p_in, p_out, z)

    %brute force newstep algo:
    %generating random integers in the range of the discretization
    R1 = randi(Nt/2 - 1) + 1;
    R2 = randi(Nt/2 - 1) + 1;
    R3 = randi(Nt/2 - 1) + 1;
    R4 = randi(Nt/2 - 1) + 1;
    R5 = randi(Nt/2 - 1) + 1;
  
    %generating random new positions for those indices:
    D1 = normrnd(pos(1,R1), sigma * 0.1);
    D2 = normrnd(pos(2,R2), sigma * 0.5);
    D3 = normrnd(pos(3,R3), sigma);
    D4 = normrnd(pos(4,R4), sigma * 0.5);
    D5 = normrnd(pos(5,R5), sigma * 0.1);
    
    %the numerical position for the random index chosen:
    t1 = pos(1,R1);
    t2 = pos(2,R2);
    t3 = pos(3,R3);
    t4 = pos(4,R4);
    t5 = pos(5,R5);
    
    %assinging the new values:
    pos(1,R1)           = D1;
    pos(5,Nt -(R1-1))   = -D1;

    pos(2,R2)           = D2;
    pos(4,Nt -(R2-1))   = -D2;

    pos(3,R3)           = D3;
    pos(3,Nt -(R3-1))   = -D3;

    pos(4,R4)           = D4;
    pos(2,Nt -(R4-1))   = -D4;
    
    pos(5,R5)           = D5;
    pos(1,Nt -(R5-1))   = -D5;
    
    while pos(1, R1) < p_in(1) || pos(1, R1) > p_out(1)            
        D                       = normrnd(t1, sigma * 0.1);
        pos(1, R1)              = D;
        pos(5, Nt -(R1 - 1))    = -D;
    end

    while pos(2, R2) < p_in(2) || pos(2, R2) > p_out(2)            
        D                       = normrnd(t2, sigma * 0.5);
        pos(2, R2)              = D;
        pos(4, Nt -(R2 - 1))    = -D;
    end
     
    while pos(3, R3) < p_in(3) || (pos(3, R3)>=0) || pos(3, R3) > (z(R3) * abs(pos(3, 1)))
        D                       = normrnd(t3, sigma);
        pos(3, R3)              = D;
        pos(3, Nt -(R3 - 1))    = -D;
    end
     
    while pos(4, R4) < p_in(4) || pos(4, R4) > p_out(4)            
        D                       = normrnd(t4, sigma * 0.5);
        pos(4, R4)              = D;
        pos(2, Nt -(R4 - 1))    = -D;
    end

    while pos(5, R5) < p_in(5) || pos(5, R5) > p_out(5)            
        D                       = normrnd(t5, sigma * 0.1);
        pos(5, R5)              = D;
        pos(1, Nt -(R5 - 1))    = -D;
    end

end



