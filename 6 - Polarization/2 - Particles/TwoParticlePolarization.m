clc
clear all

Nx      = 200;
Alpha   = 3:0.2:10;
Kappa   = -0.3:0.1:0.3;
Eta     = 20;
Beta    = 10^-5;

WaveFunction    = zeros(length(Alpha), length(Kappa), Nx^2);


Polarization    = struct();
Eq_Pos          = [];

for alphaInd = 1:length(Alpha)
    a           = Alpha(alphaInd);
    Potential   = @(x) 0.25 * (x(1)^2 - a)^2  + 0.25 * (x(2)^2 - a)^2 + Eta/(abs(x(2) - x(1)) + Beta) ;
    options     = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'MaxFunEvals', 10^9, 'MaxIter', 10^9);
    x_start     = [-1 1];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    Eq_Pos(alphaInd, :)   = x0;
    Error(alphaInd)    = fval0;
end

disp('Done with 2 particle equilibrium positions.')

figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 2 particles')
xlabel('\alpha')
ylabel('\chi_0')
plot(Alpha, Eq_Pos(:, 1), '.-', 'DisplayName', '1st particle')
plot(Alpha, Eq_Pos(:, 2), '.-', 'DisplayName', '2nd particle')
plot(Alpha, sqrt(Alpha), 'k.')
plot(Alpha, -sqrt(Alpha), 'k.')
legend
hold off

figure(2)
clf(figure(2))
hold on
title('Value of the minima of the function')
xlabel('\alpha')
ylabel('Fvals')
plot(Alpha, Error,                  'DisplayName', 'Function value')
plot(Alpha, 20./(2 * sqrt(Alpha)),  'DisplayName', 'Energy between the 2 minima')
% plot(alpha, (Eq_pos - sqrt(alpha)))
set(gca, 'Yscale', 'log')
legend
hold off

XMin    = -5;
XMax    = 5;

x = linspace(XMin, XMax, Nx);
dx = x(2) - x(1);

X1 = kron(sparse(diag(x)), speye(Nx));
X2 = kron(speye(Nx), sparse(diag(x)));
XX = full(diag(X1 + X2));




K = -1/(2 * dx^2);

for alphaInd = 1:Nx
    if alphaInd == 1
        KineticMtx(alphaInd, alphaInd)        = -2 * K;
    else
        KineticMtx(alphaInd - 1, alphaInd)    = 1 * K;
        KineticMtx(alphaInd, alphaInd - 1)    = 1 * K;
    end
end
KineticMtx = sparse(KineticMtx);

K1      = kron(KineticMtx, speye(Nx));
K2      = kron(speye(Nx), KineticMtx);
KineticMtx = K1 + K2;

for alphaInd = 1:length(Alpha)
    for kappaInd = 1:length(Kappa)
        clear Hamiltonian
        PotentialMtx = sparse(zeros(Nx));

        U = 0.25 * (x.^2 - Alpha(alphaInd)).^2 + Kappa(kappaInd) * x;
        PotentialMtx = sparse(diag(U));

        U1 = kron(PotentialMtx, speye(Nx));
        U2 = kron(speye(Nx), PotentialMtx);
        PotentialMtx = U1 + U2;

        Hamiltonian = KineticMtx + PotentialMtx;
    
        X_matrix    = sparse(diag(x));
        X1          = kron(X_matrix, speye(Nx));
        X2          = kron(speye(Nx), X_matrix);
        
        X12 = diag(X2) - diag(X1);
        
        InteractionMtx = Eta ./ sqrt(X12.^2 + Beta^2);
        InteractionMtx = sparse(1:Nx^2, 1:Nx^2, InteractionMtx);
        
        Hamiltonian = Hamiltonian + InteractionMtx;
        
        Hamiltonian = (Hamiltonian + Hamiltonian') * 0.5;
        
        Temp_Proj = zeros(Nx^2, 3);

        k = 1;
        q = 1;
        xtemp = [];
        XXX = zeros(Nx^2, 1);
        for j1 = 1:Nx
            jsize(j1) = 0;
            for j2 = 1:Nx
                j = (j1 - 1) * Nx + j2;
        
                if X12(j) < -eps
                    XXX(j)          = 1;
                    Temp_Proj(k, :) = [k, j, 1];
                    k               = k + 1;
                    jsize(j1)       = jsize(j1) + 1;
                else

                end
            end
        end
% 
%         k = 1;
%         for i = 1:length(XXX)
%             if XXX(i) == 1
%                 sel(k) = XX(i);
%                 k = k + 1;
%             end
%         end

        
        Temp_Proj = Temp_Proj(1:(k-1), :);
        Projector = sparse(Temp_Proj(:, 1), Temp_Proj(:, 2), Temp_Proj(:, 3), k - 1, Nx^2);
    
        Hamiltonian = Projector * Hamiltonian * Projector';
        X1          = Projector * X1 * Projector';
        X2          = Projector * X2 * Projector';

        [Psi, E]    = eigs(Hamiltonian, 2, 'sa');
        Spectrum    = diag(E);
        dE          = Spectrum(2) - Spectrum(1);


        Int1 = Psi(:, 1)' * (X1 * Psi(:, 1));
        Int2 = Psi(:, 1)' * (X2 * Psi(:, 1));

        P = Int1 + Int2;

        WaveFunction = Projector' * Psi(:, 1);
        WFT         = reshape(WaveFunction, [Nx Nx]);

        
        for Ind1 = 2:Nx
            PsiX(Ind1) = 0;
            PsiY(Ind1) = 0;
            for Ind2 = 1:Nx
                PsiX(Ind1) = PsiX(Ind1) + WFT(Ind1, Ind2);
                PsiY(Ind1) = PsiY(Ind1) + WFT(Ind2, Ind1);
            end
        end

        PsiX = PsiX ./norm(PsiX);
        PsiY = PsiY ./norm(PsiY);
        
        figure(5)
        clf(figure(5))
        hold on
        plot(x, 10^-3 * (0.25 * (x.^2 - Alpha(alphaInd)).^2 + Kappa(kappaInd) * x))
        plot(x, abs(PsiX))
        plot(x, abs(PsiY))

        hold off

        figure(6)
        clf(figure(6))
        hold on
        surf(x, x, abs(WFT), 'EdgeColor', 'none', 'FaceColor', 'interp')
        axis square
        hold off
    pause


        Polarization.Splittings(alphaInd, kappaInd)         = dE;
        Polarization.WaveFunction(alphaInd, kappaInd, :)   = WaveFunction;
        
        % QuantumPolarization matrix element:

        Pol(alphaInd, kappaInd) = P;


        a                       = Alpha(alphaInd);
        k                       = Kappa(kappaInd);
        Potential               = @(x) 0.25 * (x(1)^2 - a)^2 + k*x(1) + 0.25 * (x(2)^2 - a)^2 + k*x(2) + Eta/(abs(x(2) - x(1))) ;
        options                 = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'MaxFunEvals', 10^9, 'MaxIter', 10^9);
        x_start                 = [-1 1];
        [x0, fval0]             = fminsearch(Potential, x_start, options);
        Eq_Pos(alphaInd, :)     = x0;
        
        Polarization.ClassicalPolarizationMtx(alphaInd, kappaInd) = sum(x0);
        k
        a
    end
end

Polarization.Alpha              = Alpha;
Polarization.Kappa              = Kappa;
Polarization.PolarizationMtx    = Pol;
Polarization.Info               = '2 particle polarization and data structure.';


%%

figure(3)
clf(figure(3))
hold on
surf(Alpha, Kappa, Pol')

hold off