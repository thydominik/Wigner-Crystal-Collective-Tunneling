clc
clear all

Nx      = 200;
Alpha   = -1:0.2:8;
Kappa   = -0.5:0.01:0.5;

WaveFunction = zeros(length(Alpha), length(Kappa), Nx);

% 1 particle Equilibrium positions

Eq_Pos = [];

for i = 1:length(Alpha)
    a           = Alpha(i);
    Potential   = @(x) 0.25 * (x^2 - a)^2;
    options     = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'MaxFunEvals', 10^9, 'MaxIter', 10^9);
    x_start     = 0;
    [x0, fval0] = fminsearch(Potential, x_start, options);
    Eq_Pos(i)   = x0;
end

figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 1 particle')
xlabel('\alpha')
ylabel('\chi_0')
plot(Alpha, Eq_Pos, '.-')
plot(Alpha, sqrt(abs(Alpha)), 'o')
hold off

Polarization                            = struct();

% Creating the potnential
for alphaInd = 1 : length(Alpha)
    a = Alpha(alphaInd);

    for kappaInd = 1:length(Kappa)
        k = Kappa(kappaInd);
        S   = linspace(-10, 10, Nx);
        dx  = S(2) - S(1);
        V   = 0.25 * (S.^2 - a).^2 + k * S;

        PotentialMtx = sparse(diag(V));
        for i = 2:Nx
            K = 1/(2 * dx^2);
            KineticMtx(i - 1, i) = -K;
            KineticMtx(i, i - 1) = -K;
        end
        Hamiltonian = PotentialMtx + sparse(KineticMtx);
        
        [Psi, E]    = eig(full(Hamiltonian));
        E           = diag(E);
        dE          = E(2) - E(1);
        
        Polarization.Splitting(alphaInd, kappaInd) = dE;

        figure(4)
        clf(figure(4))
        hold on
        title([a k])
        ylim([0 1])
        plot(S, 5*V/max(V))
        plot(S, (Psi(:, 1)/norm(Psi(:, 1))))
        plot(S, (Psi(:, 2)/norm(Psi(:, 2))))
        yline(0)
        ylim([-0.2 0.2])
        %Fr = getframe(gcf);
        %[Im(:, :, 1, alphaInd), Map] = rgb2ind(Fr.cdata, 8);
    
        hold off

        %Polarization calculation
        Int = 0;
        for i = 1:1
            for n = 1:Nx
                Int = Int + Psi(n, i)' * S(n) * Psi(n, i) * dx;
            end
            Pol(alphaInd, kappaInd) = Int;
        end

        WaveFunction(alphaInd, kappaInd, :) = Psi(:, 1);

        % Calculating the classical equilibrium positions
        a           = Alpha(alphaInd);
        k           = Kappa(kappaInd);
        Potential   = @(x) 0.25 * (x^2 - a)^2 + k * x;
        options     = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'MaxFunEvals', 10^9, 'MaxIter', 10^9);
        x_start     = 0;
        [x0, fval0] = fminsearch(Potential, x_start, options);
        Eq_Pos(i)   = x0;
    
        
        ClassicalPol(alphaInd, kappaInd) = sum(x0);
    end
    
    

    Split(alphaInd, 1) = E(2) - E(1);
    Split(alphaInd, 2) = a;
end


figure(5)
clf(figure(5))
hold on
surf(Alpha, Kappa, Pol.','FaceAlpha',0.5)
xlabel('Alpha')
ylabel('Kappa')
zlabel('P')
hold off




Polarization.Info                       = '1 particle polarization and data structure.';
Polarization.Resolution                 = Nx;
Polarization.Alpha                      = Alpha;
Polarization.Kappa                      = Kappa;
Polarization.PolarizationMtx            = Pol;
Polarization.Splittings                 = Split(:, 1);
Polarization.GSWaveFunctionMtx          = WaveFunction;
Polarization.ClassicalPolarizationMtx   = ClassicalPol;

save('OneParticlePolarization', 'Polarization')
disp('Done with calculation')
