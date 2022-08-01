function [action, Kin] = ActionCalc(pos, r, alpha, eta, NoP, Nt, z, dz, shift)
    %ACTIONCALC: Calculating the action of a given trajectory set: {pos}
    % pos - [matrix, doulbe] gives the trajectories for NoP # of particles; r -
    % [doulbe] time scaling varieble; alpha - [double] alpha parameter; eta -
    % [double] dimless coulomb int param.; NoP - [int] # of particles; Nt -
    % [int] # of points in the trajectories; z - [array, double] imag time
    % variable; dz - [double] difference between two z times; shift - [double]
    % energy shift due to interaction from initial poistions (calculated in
    % position intialization function)

    action = 0;

    %first calculating the interaction term (without the prefactor)
    Q = zeros(Nt, 1);
    for i = 1:Nt
        for j = 1:NoP
            for k = (j+1):NoP
                Q(i) = Q(i) + eta/abs(pos(j, i) - pos(k, i));
            end
        end
    end

    %second: Potential term (without the prefactor)
    V = zeros(NoP, Nt);
    for i = 1:Nt
        for j = 1:NoP
            V(j, i) = 0.25 * (pos(j, i)^2 - alpha)^2;
        end
    end

    %integrating the potential and the interaction terms:
    for i = 2:Nt-1
        if Q(i) + V(1, i) + V(2, i) + V(3, i) + V(4, i) + V(5, i) - shift == 0
            action = action;
        elseif i > 1 || i < Nt
            PrefactorTemp   = ((1 - z(i)^2)/r);
            func_i           = sum(V(:, i)) + Q(i);
            func_ip1         = sum(V(:, i - 1)) + Q(i - 1);
            action  = action + dz * 0.5 * 1/PrefactorTemp * (func_i + func_ip1 - (2 * shift));
        end
    end

    %the derivative part:
    der     = zeros(NoP, Nt-1);
    z_temp  = z + dz/2;
    for i = 2:Nt
        for j = 1:NoP
            der(j, i) = (pos(j, i) - pos(j, i - 1))/dz;
        end
    end

    %integrating the kinetic part:
    for i = 1:Nt-1
        func = sum(der(:, i).^2);
        action = action + dz * 0.5 * (1 - z_temp(i)^2)/r * (func);
    end

end

