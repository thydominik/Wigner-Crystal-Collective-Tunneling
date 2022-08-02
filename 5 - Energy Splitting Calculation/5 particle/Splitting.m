function [shift, Spectra, splitAL, Psi, VSInterp, SS] = Splitting(Xi, EigVals, z, Alpha, R, NoP, Nt, Eta, Trajectory, Action, S, VS)
    %SPLITTING:

    % Calculate the energy shift:
    Q = 0;
    for part_i_ind = 1:NoP
        V(part_i_ind) = 0.25 * (Trajectory(part_i_ind, 1)^2 - Alpha)^2;
        for part_j_ind = (part_i_ind+1):NoP
            Q = Q + Eta/abs(Trajectory(part_i_ind, 1) - Trajectory(part_j_ind, 1));
        end
    end

    shift = Q + sum(V);

    omega = sqrt(EigVals(1, 1));

    % Calculating the Arc length 1D Schrödinger problem
    NtS = 100;
    NewS = linspace(min(S), max(S), NtS);
    VSInterp = interp1(S, VS, NewS, 'Spline');
    VSInterp = VSInterp - min(VSInterp);
    NewS = NewS - NewS/2;
    dS = NewS(2) - NewS(1);
    % New Part:
    [gof, fc] = f_fitting_VS_2(NewS(1:15), VSInterp(1:15));

    NoNP = 160; %The Number of New Points in V(S)
    SS = NewS;
    for i = 1:NoNP
        VSInterp = [fc.b * ((-i*dS + min(NewS)) - min(NewS))^4 VSInterp];
        SS = [(-i*dS + min(NewS)) SS];
        VSInterp = [VSInterp fc.b * ((i*dS + max(NewS)) - max(NewS))^4];
        SS = [SS (i*dS + max(NewS))];
    end

    % --------------------------------
    NtS = NtS + 2*NoNP;
    Hami        = zeros(NtS, NtS);
    Kinetic     = zeros(NtS, NtS);
    Potential   = zeros(NtS, NtS);

    K = 1/(2 * dS^2);

    for i = 2:NtS
        Kinetic(i - 1, i) = -K;
        Kinetic(i, i - 1) = -K;
    end

    for i = 1:NtS
        Potential(i, i) = VSInterp(i);
    end

    figure(100)
    clf(figure(100))
    hold on
    plot(SS, VSInterp, 'k.')
    plot(SS, diag(Potential), 'o')
    %plot(NewS(1:10), plotdiag(1:10), 'b.-')
    %xx = linspace(-1, -0.2, 1000);
    %plot(xx, fc.b * (xx - min(NewS)).^2)
    hold off
   

    % ---------------------------------------------------------------------

    Hamilton = Potential + Kinetic;
    diag(Hamilton)
    [Psi, Spectra] = eig(Hamilton);
    Spectra = diag(Spectra);
    figure(120)
    clf(figure(120))
    hold on
    title(Spectra(2) - Spectra(1))
    Spectra(2) - Spectra(1)
    plot(SS, 10^-2 * VSInterp)
    plot(SS, Psi(:, 1), '.-')
    plot(SS, Psi(:, 2), 'r.-')
%     plot(SS, Psi(:, 3), 'r.-')
%     plot(Psi(:, 4), 'r.-')
    yline(0)
    hold off

    disp('Hami')
    splitAL = Spectra(2) - Spectra(1);

end

