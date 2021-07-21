function khi = initpos(N,a)
    khi = zeros(N,1);
    for i = 1:N
        if i <N/2
            khi(i) = -sqrt(a);%-2*rand();
        else
            khi(i) = sqrt(a);%2*rand();
        end
    end
    khi(1) = -sqrt(a);
    khi(N) = sqrt(a);
end