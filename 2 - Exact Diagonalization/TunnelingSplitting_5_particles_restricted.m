clc
clear all

disp('5 particle tunneling splitting calculation.')
tic
%% 5 particle Equilibrium positions

Na      = 20;
alpha   = linspace(7, 20, Na);
eta     = 20;
Eq_Pos  = [];
Beta    = 10^-5;

Eq_Pos = [];

for i = 1:Na
    a = alpha(i);
    Potential = @(x) 0.25 * ((x(1)^2 - a)^2 + (x(2)^2 - a)^2 + (x(3)^2 - a)^2 + (x(4)^2 - a)^2 + (x(5)^2 - a)^2) + eta * (1/abs(x(1) - x(2)) + 1/abs(x(1) - x(3)) + 1/abs(x(1) - x(4)) + 1/abs(x(1) - x(5)) + 1/abs(x(2) - x(3)) + 1/abs(x(2) - x(4)) + 1/abs(x(2) - x(5)) + 1/abs(x(3) - x(4)) + 1/abs(x(3) - x(5)) + 1/abs(x(4) - x(5)));
    % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
    options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10^6, 'MaxIter', 10^6);
    x_start = [-sqrt(a)-1 -sqrt(a)+0.5 -sqrt(a)+1 sqrt(a)-1 sqrt(a)+1];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    Eq_Pos(i, :) = sort(x0);
    FuncVal(i) = fval0;
end

disp('Done with 5 particle')

figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 5 particle')
xlabel('\alpha')
ylabel('\chi')
plot(alpha, Eq_Pos(:, 1), '.-')
plot(alpha, Eq_Pos(:, 2), '.-')
plot(alpha, Eq_Pos(:, 3), '.-')
plot(alpha, Eq_Pos(:, 4), '.-')
plot(alpha, Eq_Pos(:, 5), '.-')
plot(alpha, sqrt(alpha), 'k')
plot(alpha, -sqrt(alpha), 'k')
yline(0)
hold off


figure(2)
clf(figure(2))
hold on
title('Value of the minima of the function')
xlabel('\alpha')
ylabel('Fvals')
plot(alpha, FuncVal,                  'DisplayName', 'Function value')
% plot(alpha, (Eq_pos - sqrt(alpha)))
%set(gca, 'Yscale', 'log')
yline(min(FuncVal))
hold off

%%
for q = 1:length(alpha)
    alpha(q)
    Nx1 = 15;
    Nx2 = 20;
    Nx3 = 40;
    Nx4 = 20;
    Nx5 = 15;

    Nx1 * Nx2 * Nx3 * Nx4 * Nx5

    a = 1.5;
    XMin1   = Eq_Pos(q, 1) - a;
    XMax1   = -Eq_Pos(q, 5) + a;
    b = 3;
    XMin2   = Eq_Pos(q, 2) - b;
    XMax2   = -Eq_Pos(q, 4) + b;
    c = 3.5;
    XMin3   = Eq_Pos(q, 3) - c;
    XMax3   = -Eq_Pos(q, 3) + c;
    d = 3;
    XMin4   = Eq_Pos(q, 4) - d;
    XMax4   = -Eq_Pos(q, 2) + d;
    e = 1.5;
    XMin5   = Eq_Pos(q, 5) - e;
    XMax5   = -Eq_Pos(q, 1) + e;
    
    x1 = linspace(XMin1, XMax1, Nx1);
    x2 = linspace(XMin2, XMax2, Nx2);
    x3 = linspace(XMin3, XMax3, Nx3);
    x4 = linspace(XMin4, XMax4, Nx4);
    x5 = linspace(XMin5, XMax5, Nx5);
    
    dx1 = x1(2) - x1(1);
    dx2 = x2(2) - x2(1);
    dx3 = x3(2) - x3(1);
    dx4 = x4(2) - x4(1);
    dx5 = x5(2) - x5(1);
    
    K1 = -1/(2 * dx1^2);
    K2 = -1/(2 * dx2^2);
    K3 = -1/(2 * dx3^2);
    K4 = -1/(2 * dx4^2);
    K5 = -1/(2 * dx5^2);
    
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
    
    for i = 1:Nx4
        if i == 1
            KineticMtx4(i, i)        = -2 * K4;
        else
            KineticMtx4(i - 1, i)    = 1 * K4;
            KineticMtx4(i, i - 1)    = 1 * K4;
        end
    end
    
    for i = 1:Nx5
        if i == 1
            KineticMtx5(i, i)        = -2 * K5;
        else
            KineticMtx5(i - 1, i)   = 1 * K5;
            KineticMtx5(i, i - 1)   = 1 * K5;
        end
    end
    KineticMtx1 = sparse(KineticMtx1);
    KineticMtx2 = sparse(KineticMtx2);
    KineticMtx3 = sparse(KineticMtx3);
    KineticMtx4 = sparse(KineticMtx4);
    KineticMtx5 = sparse(KineticMtx5);
    
    K1          = kron(kron(KineticMtx1, speye(Nx2)), kron(speye(Nx3), kron(speye(Nx4), speye(Nx5))));
    K2          = kron(speye(Nx1), kron(KineticMtx2, kron(speye(Nx3), kron(speye(Nx4), speye(Nx5)))));
    K3          = kron(speye(Nx1), kron(speye(Nx2), kron(KineticMtx3, kron(speye(Nx4), speye(Nx5)))));
    K4          = kron(speye(Nx1), kron(speye(Nx2), kron(speye(Nx3), kron(KineticMtx4, speye(Nx5)))));
    K5          = kron(speye(Nx1), kron(speye(Nx2), kron(speye(Nx3), kron(speye(Nx4), KineticMtx5))));
    KineticMtx  = K1 + K2 + K3 + K4 + K5;

    PotentialMtx1 = sparse(zeros(Nx1));
    PotentialMtx2 = sparse(zeros(Nx2));
    PotentialMtx3 = sparse(zeros(Nx3));
    PotentialMtx4 = sparse(zeros(Nx4));
    PotentialMtx5 = sparse(zeros(Nx5));

    U1 = 0.25 * (x1.^2 - alpha(q)).^2;
    U2 = 0.25 * (x2.^2 - alpha(q)).^2;
    U3 = 0.25 * (x3.^2 - alpha(q)).^2;
    U4 = 0.25 * (x4.^2 - alpha(q)).^2;
    U5 = 0.25 * (x5.^2 - alpha(q)).^2;

    PotentialMtx1 = sparse(diag(U1));
    PotentialMtx2 = sparse(diag(U2));
    PotentialMtx3 = sparse(diag(U3));
    PotentialMtx4 = sparse(diag(U4));
    PotentialMtx5 = sparse(diag(U5));
    
    U1 = kron(kron(PotentialMtx1, speye(Nx2)), kron(speye(Nx3), kron(speye(Nx4), speye(Nx5))));
    U2 = kron(kron(speye(Nx1), PotentialMtx2), kron(speye(Nx3), kron(speye(Nx4), speye(Nx5))));
    U3 = kron(kron(speye(Nx1), speye(Nx2)), kron(PotentialMtx3, kron(speye(Nx4), speye(Nx5))));
    U4 = kron(kron(speye(Nx1), speye(Nx2)), kron(speye(Nx3), kron(PotentialMtx4, speye(Nx5))));
    U5 = kron(kron(speye(Nx1), speye(Nx2)), kron(speye(Nx3), kron(speye(Nx4), PotentialMtx5)));
    PotentialMtx = U1 + U2 + U3 + U4 + U5;
    
    Hamiltonian = KineticMtx + PotentialMtx;
    
    X_matrix1   = sparse(diag(x1));
    X_matrix2   = sparse(diag(x2));
    X_matrix3   = sparse(diag(x3));
    X_matrix4   = sparse(diag(x4));
    X_matrix5   = sparse(diag(x5));

    X1  = kron(X_matrix1, kron(speye(Nx2), kron(speye(Nx3), kron(speye(Nx4), speye(Nx5)))));
    X2  = kron(speye(Nx1), kron(X_matrix2, kron(speye(Nx3), kron(speye(Nx4), speye(Nx5)))));
    X3  = kron(speye(Nx1), kron(speye(Nx2), kron(X_matrix3, kron(speye(Nx4), speye(Nx5)))));
    X4  = kron(speye(Nx1), kron(speye(Nx2), kron(speye(Nx3), kron(X_matrix4, speye(Nx5)))));
    X5  = kron(speye(Nx1), kron(speye(Nx2), kron(speye(Nx3), kron(speye(Nx4), X_matrix5))));

    X12 = diag(X1) - diag(X2);
    X13 = diag(X1) - diag(X3);
    X14 = diag(X1) - diag(X4);
    X15 = diag(X1) - diag(X5);
    X23 = diag(X2) - diag(X3);
    X24 = diag(X2) - diag(X4);
    X25 = diag(X2) - diag(X5);
    X34 = diag(X3) - diag(X4);
    X35 = diag(X3) - diag(X5);
    X45 = diag(X4) - diag(X5);

    InteractionMtx = eta ./ sqrt(X12.^2 + Beta^2) + eta ./ sqrt(X13.^2 + Beta^2) + eta ./ sqrt(X14.^2 + Beta^2) + eta ./ sqrt(X15.^2 + Beta^2) + eta ./ sqrt(X23.^2 + Beta^2) + eta ./ sqrt(X24.^2 + Beta^2) + eta ./ sqrt(X25.^2 + Beta^2) + eta ./ sqrt(X34.^2 + Beta^2) + eta ./ sqrt(X35.^2 + Beta^2) + eta ./ sqrt(X45.^2 + Beta^2);
    InteractionMtx = sparse(1:(Nx1 * Nx2 * Nx3 * Nx4 * Nx5), 1:(Nx1 * Nx2 * Nx3 * Nx4 * Nx5), InteractionMtx);
    
    Hamiltonian = Hamiltonian + InteractionMtx;
    
    Hamiltonian = (Hamiltonian + Hamiltonian') * 0.5;
    
    Temp_Proj = zeros((Nx1 * Nx2 * Nx3 * Nx4 * Nx5), 3);
    
    k = 1;
    for j = 1:(Nx1 * Nx2 * Nx3 * Nx4 * Nx5)
        if (X12(j) < -eps) && (X23(j) < -eps) &&  (X34(j) < -eps) && (X45(j) < -eps)
            Temp_Proj(k, :) = [k, j, 1];
            k = k + 1;
        end
    end

    Temp_Proj = Temp_Proj(1:(k-1), :);
    Projector = sparse(Temp_Proj(:, 1), Temp_Proj(:, 2), Temp_Proj(:, 3), k - 1, (Nx1 * Nx2 * Nx3 * Nx4 * Nx5));

    Hamiltonian = Projector * Hamiltonian * Projector';

    [Psi, E]    = eigs(Hamiltonian, 2, 'sa');
    Spectrum    = diag(E);
    alpha(q)
    dE(q)       = Spectrum(2) - Spectrum(1)
end
%%
toc

%%
pa = [7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5];
ps = [1.77 1.57 1.33 1.068 0.76 0.48 0.23 0.0756 0.0275 0.01548 0.015 0.01803];
figure(3)
%clf(figure(3))
hold on
ylabel('Split')
xlabel('\alpha')
plot(alpha, dE, 'o')
plot(pa, ps, '.-')
set(gca, 'Yscale', 'log')
hold off

%%
figure(5)
clf(figure(5))
hold on
plot(x1, ones(length(x1), 1), '.-')
plot(x2, 0.01 + ones(length(x2), 1), '.-')
plot(x3, 0.02 + ones(length(x3), 1), '.-')
plot(x4, 0.01 + ones(length(x4), 1), '.-')
plot(x5, 0.00 + ones(length(x5), 1), '.-')
ylim([0.99 1.5])
hold off



data(:, 1) = alpha;
data(:, 2) = dE;

name = ['EDSplitting_5_particles_restricted_Nx1_' num2str(Nx1) '_Nx2_' num2str(Nx2) '_Nx3_' num2str(Nx3) '_Nx4_' num2str(Nx4) '_Nx5_' num2str(Nx5) '_beta_' num2str(Beta)];
save(name, 'data');

