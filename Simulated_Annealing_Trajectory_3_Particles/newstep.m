function pos = newstep(pos,N,sigma,in1,fin1,in2,fin2,in3,fin3)

    R1 = randi(N - 1) + 1;
    R2 = randi(N - 2) + 1;
    R3 = randi(N - 1) + 1;
    
    D1 = normrnd(pos(1,R1),sigma);
    D2 = normrnd(pos(2,R2),sigma);
    D3 = normrnd(pos(3,R3),sigma);
    
    t1 = pos(1,R1);
    t2 = pos(2,R2);
    t3 = pos(3,R3);
    
    pos(1,R1) = D1;
    %pos(1,N-(R1-1)) = -D3;
    pos(2,R2) = D2;
    %pos(2,N-(R2-1)) = -D2;
    pos(3,R3) = D3;
    %pos(3,N-(R3-1)) = -D1;
    
    while pos(1,R1) < in1 || pos(1,R1) > fin1            
        D = normrnd(t1,sigma);
        pos(1,R1) = D;
        %pos(3,N-(R1-1)) = -D;
    end
%     
    while pos(2,R2) < in2 || pos(2,R2) > 0 %|| (R2 <= N/2 && pos(2,R2)>=0) || (R2 >= N/2 +1 && pos(2,R2)<=0)
        D = normrnd(t2,sigma);
        pos(2,R2) = D;
        %pos(2,N-(R2-1)) = -D;
    end
%     
    while pos(3,R3) < in3 || pos(3,R3) > fin3
        D = normrnd(t3,sigma);
        pos(3,R3) = D;
        %pos(1,N-(R3-1)) = -D;
    end
%     

end