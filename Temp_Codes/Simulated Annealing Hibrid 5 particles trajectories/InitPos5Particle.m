function [pos, shift] = InitPos5Particle(NoP, Nt, p_in, p_out, eta, alpha)
    %INITPOS5PARTICLE: Initialize the starting state for 3 particles
    %NoP - [int] Number of particles; Nt - [int] # of points in the trajectory; p_in -
    %[array, double] gives the initial poistions of the particles at z = -1; p_out -
    %[array, doulbe] -||- final at z = 1; eta - [double] dimless Coulomb int.;
    %alpha - [double] current alpha param.
    
    pos = zeros(NoP, Nt);
    
    for pInd = 1 : NoP
        pos(pInd, :) = linspace(p_in(pInd), p_out(pInd), Nt);
    end

    % Calculating the energy shift needed
        % The F's are coming from the quartic potential:
        % The Q's are from the interaction:
    Q = 0;
    for i = 1:NoP
        F(i) = 0.25 * (pos(i,1)^2 - alpha)^2;
        for j = (i + 1):NoP
            Q = Q + eta/(abs(pos(i,1) - pos(j,1)));
        end
    end
    shift =  Q + sum(F);

end

