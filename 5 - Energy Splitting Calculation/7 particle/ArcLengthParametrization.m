function [ArcLength, VS] = ArcLengthParametrization(Trajectories, NoP, Nt, Alpha, Eta)
    %ARCLENGTHPARAMETRIZATION: requires a set of two dimensional curves, and
    %calculates the arc length parametrized version of that: f(x) -> f(x(s))
        
    ArcLength   = 0;
    dS          = 0;
    for i = 2:Nt
        for part_ind = 1:NoP
            x(part_ind) = (Trajectories(part_ind, i) - Trajectories(part_ind, i - 1))^2;
        end
        
        dS(i)           = sqrt(sum(x));
        ArcLength(i)    = ArcLength(i - 1) + dS(i);
    end
    
    for i = 1:Nt
        VS(i) = 0;
        for part_i_ind = 1:NoP
            VS(i) = VS(i) + 0.25 * (Trajectories(part_i_ind, i)^2 - Alpha)^2;
            for part_j_ind = (part_i_ind + 1):NoP
                VS(i) = VS(i) + Eta/abs(Trajectories(part_i_ind, i) - Trajectories(part_j_ind, i));
            end
        end
    end
    
    % shifting the arc length to have z = 0 at S = 0
    ArcLength = ArcLength - max(ArcLength)/2;

end

