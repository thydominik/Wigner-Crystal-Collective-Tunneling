function [pos, sigma] = new_pos_7(xpos,T,a,c,el_stre,sigma)
    %Choosing a random electron
    j = randi(length(xpos));
    x_chosen = xpos(j);

    x_new = normrnd(x_chosen,sigma);

    %nem engedem meg hogy a r�szecsk�k �tl�pjenek egym�son
    if j>1 && j<length(xpos)
        while xpos(j+1) <= x_new || xpos(j-1) >= x_new
            x_new=normrnd(x_chosen,sigma);
        end
    elseif j == length(xpos) && length(xpos)~=1
        while xpos(j-1) >= x_new
            x_new=normrnd(x_chosen,sigma);
        end
    elseif length(xpos) ~= 1
        while xpos(2) <= x_new 
            x_new=normrnd(x_chosen,sigma);
        end
    end

    xpos(j) = x_new;
    pos = xpos;

end
