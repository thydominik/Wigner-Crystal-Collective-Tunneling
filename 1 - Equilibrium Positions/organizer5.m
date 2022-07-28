function [eq_pos] = organizer5(eq_pos)
eqpos = eq_pos;
    for i = 1:length(eq_pos)
        if eq_pos(3,i) > 0
            eq_pos(1,i) = -eqpos(5,i);
            eq_pos(2,i) = -eqpos(4,i);
            eq_pos(3,i) = -eqpos(3,i);
            eq_pos(4,i) = -eqpos(2,i);
            eq_pos(5,i) = -eqpos(1,i);
        end
    end
    
end