clc
clear all

disp('3 particle tunneling splitting calculation.')
tic
% 3 particle Equilibrium positions

Na      = 50;
alpha   = linspace(2, 15, Na);
eta     = 20;
Eq_Pos  = [];
Beta    = 0.1;

Eq_pos_3 = [];

for i = 1:Na
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
Nx1 = 40;
Nx2 = 90;
Nx3 = 40;
Nx1 * Nx2 * Nx3

XMin1    = -6;
XMax1    = 0;

XMin2    = -5;
XMax2    = 5;

XMin3    = 0;
XMax3    = 6;

x1 = linspace(XMin1, XMax1, Nx1);
x2 = linspace(XMin2, XMax2, Nx2);
x3 = linspace(XMin3, XMax3, Nx3);
dx1 = x1(2) - x1(1);
dx2 = x2(2) - x2(1);
dx3 = x3(2) - x3(1);

K1 = -1/(2 * dx1^2);
K2 = -1/(2 * dx2^2);
K3 = -1/(2 * dx3^2);

for i = 1:Nx1
    if i == 1
        KineticMtx1(i, i)        = -2 * K1;
    else
        KineticMtx1(i - 1, i)    = 1 * K1;
        KineticMtx1(i, i - 1)    = 1 * K1;
    end
end

for i = 1:Nx2
    if i == 1
        KineticMtx2(i, i)        = -2 * K2;
    else
        KineticMtx2(i - 1, i)    = 1 * K2;
        KineticMtx2(i, i - 1)    = 1 * K2;
    end
end

for i = 1:Nx3
    if i == 1
        KineticMtx3(i, i)        = -2 * K3;
    else
        KineticMtx3(i - 1, i)    = 1 * K3;
        KineticMtx3(i, i - 1)    = 1 * K3;
    end
end
KineticMtx1 = sparse(KineticMtx1);
KineticMtx2 = sparse(KineticMtx2);
KineticMtx3 = sparse(KineticMtx3);

K1      = kron(kron(KineticMtx1, speye(Nx2)), speye(Nx3));
K2      = kron(speye(Nx1), kron(KineticMtx2, speye(Nx3)));
K3      = kron(speye(Nx1), kron(speye(Nx2), KineticMtx3));
KineticMtx = K1 + K2 + K3;

for i = 1:length(alpha)
    PotentialMtx1 = sparse(zeros(Nx1));
    PotentialMtx2 = sparse(zeros(Nx2));
    PotentialMtx3 = sparse(zeros(Nx3));
    
    U1 = 0.25 * (x1.^2 - alpha(i)).^2;
    U2 = 0.25 * (x2.^2 - alpha(i)).^2;
    U3 = 0.25 * (x3.^2 - alpha(i)).^2;
    PotentialMtx1 = sparse(diag(U1));
    PotentialMtx2 = sparse(diag(U2));
    PotentialMtx3 = sparse(diag(U3));
    
    U1 = kron(kron(PotentialMtx1, speye(Nx2)), speye(Nx3));
    U2 = kron(speye(Nx1), kron(PotentialMtx2, speye(Nx3)));
    U3 = kron(kron(speye(Nx1), speye(Nx2)), PotentialMtx3);
    PotentialMtx = U1 + U2 + U3;
    
    Hamiltonian = KineticMtx + PotentialMtx;
    
    X_matrix1    = sparse(diag(x1));
    X_matrix2    = sparse(diag(x2));
    X_matrix3    = sparse(diag(x3));
    X1          = kron(X_matrix1, kron(speye(Nx2), speye(Nx3)));
    X2          = kron(speye(Nx1), kron(X_matrix2, speye(Nx3)));
    X3          = kron(kron(speye(Nx1), speye(Nx2)), X_matrix3);
    
    X12 = diag(X1) - diag(X2);
    X13 = diag(X1) - diag(X3);
    X23 = diag(X2) - diag(X3);

    InteractionMtx = eta ./ sqrt(X12.^2 + Beta^2) + eta ./ sqrt(X13.^2 + Beta^2) + eta ./ sqrt(X23.^2 + Beta^2);
    InteractionMtx = sparse(1:(Nx1 * Nx2 * Nx3), 1:(Nx1 * Nx2 * Nx3), InteractionMtx);
    
    Hamiltonian = Hamiltonian + InteractionMtx;
    
    Hamiltonian = (Hamiltonian + Hamiltonian') * 0.5;
    
    Temp_Proj = zeros((Nx1 * Nx2 * Nx3), 3);
    
    k = 1;
    for j = 1:(Nx1 * Nx2 * Nx3)
        if (X12(j) < -eps) && (X23(j) < -eps)
            Temp_Proj(k, :) = [k, j, 1];
            k = k + 1;
        end
    end

    Temp_Proj = Temp_Proj(1:(k-1), :);
    Projector = sparse(Temp_Proj(:, 1), Temp_Proj(:, 2), Temp_Proj(:, 3), k - 1, (Nx1 * Nx2 * Nx3));

    Hamiltonian = Projector * Hamiltonian * Projector';

    [Psi, E]    = eigs(Hamiltonian, 2, 'sa');
    Spectrum    = diag(E);
    dE(i)       = Spectrum(2) - Spectrum(1);
end
%%
toc
figure(3)
clf(figure(3))
hold on
plot(alpha, dE)
set(gca, 'Yscale', 'log')
hold off

data.alpha = alpha;
data.EnergySplitting = dE;
data.Particle1space = x1;
data.Particle2space = x2;
data.Particle3space = x3;

name = ['EDSplitting_3_particles_restricted_Nx1_' num2str(Nx1) '_Nx2_' num2str(Nx2) '_Nx3_' num2str(Nx3) '_beta_' '0_1'];
save(name, 'data');



