clc
clear all

disp('3 particle tunneling splitting calculation.')
tic
%% 3 particle Equilibrium positions

Nx      = 80;
Nx^3
alpha   = linspace(2, 15, Nx);
eta     = 20;
Eq_Pos  = [];
Beta    = 10^-5;

Eq_pos_3 = [];

for i = 1:Nx
    a = alpha(i);
    Potential = @(x) 0.25 * (x(1)^2 - a)^2 + 0.25 * (x(2)^2 - a)^2 + 0.25 * (x(3)^2 - a)^2 + eta/abs(x(1) - x(2)) + eta/abs(x(1) - x(3)) + eta/abs(x(2) - x(3));
    % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
   options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10^6, 'MaxIter', 10^6);
    x_start = [-sqrt(a)-1 -sqrt(a)+1 sqrt(a)];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    Eq_pos_3(i, :) = x0;
    FuncVal(i) = fval0;
end

disp('Done with 3 particles')

figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 3 particles')
xlabel('\alpha')
ylabel('\chi_0')
plot(alpha, Eq_pos_3(:, 1), '.-', 'DisplayName', '1st particle')
plot(alpha, Eq_pos_3(:, 2), '.-', 'DisplayName', '2nd particle')
plot(alpha, Eq_pos_3(:, 3), '.-', 'DisplayName', '3rd particle')
plot(alpha, sqrt(alpha), 'k')
plot(alpha, -sqrt(alpha), 'k')
yline(0)
legend
hold off


figure(2)
clf(figure(2))
hold on
title('Value of the minima of the function')
xlabel('\alpha')
ylabel('Fvals')
plot(alpha, FuncVal,                  'DisplayName', 'Function value')
% plot(alpha, (Eq_pos - sqrt(alpha)))
set(gca, 'Yscale', 'log')
yline(min(FuncVal))
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

K1      = kron(kron(KineticMtx, speye(Nx)), speye(Nx));
K2      = kron(speye(Nx), kron(KineticMtx, speye(Nx)));
K3      = kron(speye(Nx), kron(speye(Nx), KineticMtx));
KineticMtx = K1 + K2 + K3;

for i = 1:length(alpha)
    PotentialMtx = sparse(zeros(Nx));
    
    U = 0.25 * (x.^2 - alpha(i)).^2;
    PotentialMtx = sparse(diag(U));
    
    U1 = kron(kron(PotentialMtx, speye(Nx)), speye(Nx));
    U2 = kron(speye(Nx), kron(PotentialMtx, speye(Nx)));
    U3 = kron(kron(speye(Nx), speye(Nx)), PotentialMtx);
    PotentialMtx = U1 + U2 + U3;
    
    Hamiltonian = KineticMtx + PotentialMtx;
    
    X_matrix    = sparse(diag(x));
    X1          = kron(X_matrix, kron(speye(Nx), speye(Nx)));
    X2          = kron(speye(Nx), kron(X_matrix, speye(Nx)));
    X3          = kron(kron(speye(Nx), speye(Nx)), X_matrix);
    
    X12 = diag(X1) - diag(X2);
    X13 = diag(X1) - diag(X3);
    X23 = diag(X2) - diag(X3);

    InteractionMtx = eta ./ sqrt(X12.^2 + Beta^2) + eta ./ sqrt(X13.^2 + Beta^2) + eta ./ sqrt(X23.^2 + Beta^2);
    InteractionMtx = sparse(1:Nx^3, 1:Nx^3, InteractionMtx);
    
    Hamiltonian = Hamiltonian + InteractionMtx;
    
    Hamiltonian = (Hamiltonian + Hamiltonian') * 0.5;
    
    Temp_Proj = zeros(Nx^3, 3);
    
    k = 1;
    for j = 1:Nx^3
        if (X12(j) < -eps) && (X23(j) < -eps)
            Temp_Proj(k, :) = [k, j, 1];
            k = k + 1;
        end
    end

    Temp_Proj = Temp_Proj(1:(k-1), :);
    Projector = sparse(Temp_Proj(:, 1), Temp_Proj(:, 2), Temp_Proj(:, 3), k - 1, Nx^3);

    Hamiltonian = Projector * Hamiltonian * Projector';

    [Psi, E]    = eigs(Hamiltonian, 2, 'sa');
    Spectrum    = diag(E);
    dE(i)       = Spectrum(2) - Spectrum(1);
end
%%
toc
figure(3)
%clf(figure(3))
hold on
plot(alpha, dE, '.-', 'DisplayName', 'Nx = 80')
set(gca, 'Yscale', 'log')
legend
hold off

data(:, 1) = alpha;
data(:, 2) = dE;

name = ['EDSplitting_3_particles_restricted_Nx_' num2str(Nx) '_beta_' num2str(Beta)];
save(name, 'data');



