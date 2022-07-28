function force = EL(xpos,a,N,el_stre)

    for i=1:N
        force = 0;
        force = xpos(i)^3 + a*xpos(i);
        for j = 1:N
            if j ~= i
                force = force +  el_stre*((xpos(i) - xpos(j))/((xpos(i) - xpos(j))^3));
            end
        end
        force
    end

end