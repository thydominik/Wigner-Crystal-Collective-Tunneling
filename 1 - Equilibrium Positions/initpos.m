function pos = initpos(N)
    if N == 3
        pos(1) = -5;
        pos(2) = -3;
        pos(3) = 5;
    else 
        pos=[];
        a = [-5,5];
        for i=1:N
            pos(i)=(a(2)-a(1)).*rand(1,1) + a(1);       
        end
        pos=sort(pos);
    end
end