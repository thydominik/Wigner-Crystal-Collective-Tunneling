clc
clear all

disp('2 particle tunneling splitting calculation.')

% 2 particle Equilibrium positions

Nx = 200;
alpha = linspace(0, 15, Nx);
eta = 20;
Eq_Pos = [];
Beta = 10^-5;

for i = 1:Nx
    a           = alpha(i);
    Potential   = @(x) 0.25 * (x(1)^2 - a)^2 + 0.25 * (x(2)^2 - a)^2 + eta/(abs(x(2) - x(1)) + Beta) ;
    options     = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'MaxFunEvals', 10^9, 'MaxIter', 10^9);
    x_start     = [-1 1];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    Eq_Pos(i, :)   = x0;
    Error(i)    = fval0;
end

disp('Done with 2 particle equilibrium positions.')

figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 2 particles')
xlabel('\alpha')
ylabel('\chi_0')
plot(alpha, Eq_Pos(:, 1), '.-', 'DisplayName', '1st particle')
plot(alpha, Eq_Pos(:, 2), '.-', 'DisplayName', '2nd particle')
plot(alpha, sqrt(alpha), 'k.')
plot(alpha, -sqrt(alpha), 'k.')
legend
hold off

figure(2)
clf(figure(2))
hold on
title('Value of the minima of the function')
xlabel('\alpha')
ylabel('Fvals')
plot(alpha, Error,                  'DisplayName', 'Function value')
plot(alpha, 20./(2 * sqrt(alpha)),  'DisplayName', 'Energy between the 2 minima')
% plot(alpha, (Eq_pos - sqrt(alpha)))
set(gca, 'Yscale', 'log')
legend
hold off

%%

XMin    = -10;
XMax    = 10;

x = linspace(XMin, XMax, Nx);
dx = x(2) - x(1);

K = -1/(2 * dx^2);

for i = 1:Nx
    if i == 1
        KineticMtx(i, i)        = -2 * K;
    else
        KineticMtx(i - 1, i)    = 1 * K;
        KineticMtx(i, i - 1)    = 1 * K;
    end


end
KineticMtx = sparse(KineticMtx);

K1      = kron(KineticMtx, speye(Nx));
K2      = kron(speye(Nx), KineticMtx);
KineticMtx = K1 + K2;

for i = 1:length(alpha)
    PotentialMtx = sparse(zeros(Nx));
    
    U = 0.25 * (x.^2 - alpha(i)).^2;
    PotentialMtx = sparse(diag(U));
    
    U1 = kron(PotentialMtx, speye(Nx));
    U2 = kron(speye(Nx), PotentialMtx);
    PotentialMtx = U1 + U2;
    
    Hamiltonian = KineticMtx + PotentialMtx;
    
    X_matrix    = sparse(diag(x));
    X1          = kron(X_matrix, speye(Nx));
    X2          = kron(speye(Nx), X_matrix);
    
    X12 = diag(X2) - diag(X1);
    
    InteractionMtx = eta ./ sqrt(X12.^2 + Beta^2);
    InteractionMtx = sparse(1:Nx^2, 1:Nx^2, InteractionMtx);
    
    Hamiltonian = Hamiltonian + InteractionMtx;
    
    Hamiltonian = (Hamiltonian + Hamiltonian') * 0.5;
    
    Temp_Proj = zeros(Nx^2, 3);
    
    k = 1;
    for j = 1:Nx^2
        if X12(j) < -eps
            Temp_Proj(k, :) = [k, j, 1];
            k = k + 1;
        end
    end

    Temp_Proj = Temp_Proj(1:(k-1), :);
    Projector = sparse(Temp_Proj(:, 1), Temp_Proj(:, 2), Temp_Proj(:, 3), k - 1, Nx^2);

    Hamiltonian = Projector * Hamiltonian * Projector';

    [Psi, E]    = eigs(Hamiltonian, 2, 'sa');
    Spectrum    = diag(E);
    dE(i)       = Spectrum(2) - Spectrum(1);
end
%%
figure(3)
clf(figure(3))
hold on
plot(alpha, dE)
hold off


data(:, 1) = alpha;
data(:, 2) = dE;

name = ['EDSplitting_2_particles_Nx_' num2str(Nx)];
save(name, 'data');


