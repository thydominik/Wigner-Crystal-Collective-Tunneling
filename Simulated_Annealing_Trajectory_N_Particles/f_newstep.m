function pos = f_newstep(pos, N, Nt, sigma, p_in, p_out)

    %brute force newstep algo:
    %generating random integers in the range of the discretization
    for i = 1:N
        RI(i) = randi(Nt/2 - 1) + 1;
        DI(i) = normrnd(pos(i, RI(i)), sigma);
        TI(i) = pos(i, RI(i));
    end
    
    for i = 1:N
        pos(i, RI(i))                       = DI(i);
        pos(N - i + 1, Nt - (RI(i) - 1))    = - DI(i);
    end
    
    for i = 1:N
        while pos(i, RI(i)) < p_in(i) || pos(i, RI(i)) > p_out(i)
            D                                   = normrnd(TI(i), sigma);
            pos(i, RI(i))                       = D;
            pos(N - i + 1, Nt - (RI(i) - 1))    = -D;
        end
        if i == floor(N/2) + 1
            while pos(i, RI(i)) < p_in(i) || pos(i, RI(i)) > p_out(i) || pos(i, RI(i)) >= 0
            D                                   = normrnd(TI(i), sigma);
            pos(i, RI(i))                       = D;
            pos(N - i + 1, Nt - (RI(i) - 1))    = -D;
            end
        end
    end
end

