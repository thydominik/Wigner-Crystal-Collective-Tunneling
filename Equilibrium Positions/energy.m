function E=energy(pos,N,a,c,el_stre,E_0)
    E = 0;
    for i = 1:N
        E = E +  ((a/2 * pos(i)^2) + (1/4) * pos(i)^4 );
        for j = 1:N
            if j ~= i 
                E = E + (1/2)*(el_stre / abs(pos(i)-pos(j)));
            end
        end
    end 
end