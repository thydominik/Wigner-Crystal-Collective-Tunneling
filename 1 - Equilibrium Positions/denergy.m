function dE=denergy(xpos,a,c,el_stre)
    dE=0;
    for i=1:length(xpos)
        dE = dE + 2*a*xpos(i) + 4*xpos(i)^3 + c;
        for j=i:length(xpos)
            if j~=i
                dE = dE + (el_stre*(xpos(i)-xpos(j))/((xpos(i)-xpos(j))^3));
            end
        end
    end
    dE=abs(dE);
end