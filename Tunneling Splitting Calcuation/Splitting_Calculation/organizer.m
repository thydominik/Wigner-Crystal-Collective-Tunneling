function [eq_pos] = organizer(eq_pos)
eqpos = eq_pos;
    for i = 1:length(eq_pos)
        if eq_pos(2,i) > 0
            eq_pos(1,i) = -eqpos(3,i);
            eq_pos(2,i) = -eqpos(2,i);
            eq_pos(3,i) = -eqpos(1,i);
        end
    end
    
end

