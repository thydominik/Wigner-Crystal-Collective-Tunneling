function [Split_EasyMilnikov, Split_Landau, Split_Instanton, Split_Prop_a, Split_Prop_b] = Splitting(Xi, EigVals, z, Alpha, r, NoP, Nt, Eta, Trajectory, Action, S, VS)
    %SPLITTING:

    Omega = diag(sqrtm(EigVals));
    Omega_0 = min(Omega);

    % Easy Milnikov:
    Split_EasyMilnikov = 4 * Omega_0 * max(S) * sqrt(Omega_0/pi) * exp(-Action);

    % Landau Approx:
    Split_Landau = Omega_0/pi * exp(-Action);

    %Instanton Prefactor:
    Split_Instanton = 2 * sqrt((6 * Action)/pi) * Omega_0 * exp(-Action);

    %Proper Milnikov
    P = sqrt(2 * VS);
    
    %calculating the derivative
    dP_s(1) = (P(3) - P(1))/(S(3) - S(1));
    for i = 2:length(P)-1
        dP_s(i) = (P(i + 1) - P(i - 1))/(S(i + 1) - S(i - 1));
    end
    dP_s(end+1) = (P(end) - P(end - 2))/(S(end) - S(end - 2));

    % Integral:
    %Integral:
    int     = 0;
    int2    = 0;
    func    = (r./(1 - z.^2)) .* (Omega_0 - dP_s);
    func2   = (r./(1 - z.^2)) .* (dP_s(1) - dP_s);

    for i = 2:length(P)/2 - 1
        dz = (z(i + 1) - z(i));
        int = int +  dz * 0.5 * (func(i + 1) + func(i));
        int2 = int2 + dz * 0.5 * (func2(i + 1) + func2(i));
    end

    Split_Prop_a = sqrt(4 * Omega_0) / pi * max(P(1:end/2-1)) *  exp(int) * exp(-Action);
    Split_Prop_b = sqrt(4 * Omega_0) / pi * max(P) *  exp(int2) * exp(-Action);


end


%     % Calculate the energy shift:
%     Q = 0;
%     for part_i_ind = 1:NoP
%         V(part_i_ind) = 0.25 * (Trajectory(part_i_ind, 1)^2 - Alpha)^2;
%         for part_j_ind = (part_i_ind+1):NoP
%             Q = Q + Eta/abs(Trajectory(part_i_ind, 1) - Trajectory(part_j_ind, 1));
%         end
%     end
% 
%     shift = Q + sum(V);
% 
%     omega = sqrt(EigVals(1, 1));
% 
%     % Calculating the Arc length 1D Schrödinger problem
%     NtS = 100;
%     NewS = linspace(min(S), max(S), NtS);
%     VSInterp = interp1(S, VS, NewS, 'Spline');
%     VSInterp = VSInterp - min(VSInterp);
%     NewS = NewS - NewS/2;
%     dS = NewS(2) - NewS(1);
%     % New Part:
%     [gof, fc] = f_fitting_VS_2(NewS(1:15), VSInterp(1:15));
% 
%     NoNP = 160; %The Number of New Points in V(S)
%     SS = NewS;
%     for i = 1:NoNP
%         VSInterp = [fc.b * ((-i*dS + min(NewS)) - min(NewS))^2 VSInterp];
%         SS = [(-i*dS + min(NewS)) SS];
%         VSInterp = [VSInterp fc.b * ((i*dS + max(NewS)) - max(NewS))^2];
%         SS = [SS (i*dS + max(NewS))];
%     end
% 
%     % --------------------------------
%     NtS = NtS + 2*NoNP;
%     Hami        = zeros(NtS, NtS);
%     Kinetic     = zeros(NtS, NtS);
%     Potential   = zeros(NtS, NtS);
% 
%     K = 1/(2 * dS^2);
% 
%     for i = 2:NtS
%         Kinetic(i - 1, i) = -K;
%         Kinetic(i, i - 1) = -K;
%     end
% 
%     for i = 1:NtS
%         Potential(i, i) = VSInterp(i);
%     end
% 
%     figure(100)
%     clf(figure(100))
%     hold on
%     plot(SS, VSInterp, 'k.')
%     plot(SS, diag(Potential), 'o')
%     %plot(NewS(1:10), plotdiag(1:10), 'b.-')
%     %xx = linspace(-1, -0.2, 1000);
%     %plot(xx, fc.b * (xx - min(NewS)).^2)
%     hold off
%    
% 
%     % ---------------------------------------------------------------------
% 
%     Hamilton = Potential + Kinetic;
%     diag(Hamilton)
%     [Psi, Spectra] = eig(Hamilton);
%     Spectra = diag(Spectra);
%     figure(120)
%     clf(figure(120))
%     hold on
%     title(Spectra(2) - Spectra(1))
%     Spectra(2) - Spectra(1)
%     plot(SS, 10^-2 * VSInterp)
%     plot(SS, Psi(:, 1), '.-')
%     plot(SS, Psi(:, 2), 'r.-')
% %     plot(SS, Psi(:, 3), 'r.-')
% %     plot(Psi(:, 4), 'r.-')
%     yline(0)
%     hold off
% 
%     disp('Hami')
%     splitAL = Spectra(2) - Spectra(1);


