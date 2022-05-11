function [eq_pos] = organizer(eq_pos, particles)
%ORGANIZER: We have a Z2 symmetry in the code, this function aims to
%organize the data so that floor(N/2) particles are on the right side and
%the mayority of them are on the left.

%Inputs:
    % eq_pos - Calculated equlibrium positions
    % particles - # of particles

    eqpos = eq_pos;
    if particles == 1
        for i = 1:length(eq_pos)
            if eq_pos(i) > 0
                eq_pos(i) = -eqpos(i);
            end
        end
    elseif particles == 3
        for i = 1:length(eq_pos)
            if eq_pos(2,i) > 0
                eq_pos(1,i) = -eqpos(3,i);
                eq_pos(2,i) = -eqpos(2,i);
                eq_pos(3,i) = -eqpos(1,i);
            end
        end
    elseif particles == 5
        for i = 1:length(eq_pos)
            if eq_pos(3,i) > 0
                eq_pos(1,i) = -eqpos(5,i);
                eq_pos(2,i) = -eqpos(4,i);
                eq_pos(3,i) = -eqpos(3,i);
                eq_pos(4,i) = -eqpos(2,i);
                eq_pos(5,i) = -eqpos(1,i);
            end
        end
    elseif particles == 7
        for i = 1:length(eq_pos)
            if eq_pos(4,i) > 0
                eq_pos(1,i) = -eqpos(7,i);
                eq_pos(2,i) = -eqpos(6,i);
                eq_pos(3,i) = -eqpos(5,i);
                eq_pos(4,i) = -eqpos(4,i);
                eq_pos(5,i) = -eqpos(3,i);
                eq_pos(6,i) = -eqpos(2,i);
                eq_pos(7,i) = -eqpos(1,i);
            end
        end
    elseif particles > 7
        error('Uh oh, no organizing function for particle number larger than 7')
    else
        warning('Even number of particles: No organizing needed!')
    end

end

