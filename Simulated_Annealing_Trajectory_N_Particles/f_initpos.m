function [pos, shift] = f_initpos(N, Nt, p_in, p_out, rs, alpha)
%whats what:
    %N - # of particles
    %Nt - # of points on the trajectory
    %p_in - The starting positions of the particles
    %p_out - final positions of the particles
    %rs - interaction strength
    %alpha - alpha :)

    pos = zeros(N, Nt);
    
    
    for i = 1:N 
        pos(i, 1:end/2)     = p_in(i);
        pos(i, end/2 : end) = p_out(i);
    end
    
    %Calculating the energy shift needed
    % The F's are coming from the quartic potential:
    % The Q's are from the interaction:
    Q = 0;
    for i = 1:N
        F(i) = 0.25 * (pos(i,1)^2 + alpha)^2;
        for j = (i + 1):N
            Q = Q + rs/(abs(pos(i,1) - pos(j,1)));
        end
    end
    shift =  Q + sum(F);
end

