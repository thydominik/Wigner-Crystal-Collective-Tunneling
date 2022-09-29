clc
clear all
FGS = 1;
disp('7 particle tunneling splitting calculation.')
tic
% 7 particle Equilibrium positions

Na      = 10;
alpha   = linspace(11, 17, Na);
Eta     = 20;
Eq_Pos  = [];
Beta    = 0.01;
% Calculating the Equilibrium positions

EqPos = [];

N = length(alpha);

EqPos = [];

for i = 1:N
    a = alpha(i);
    Potential = @(x) 0.25 * ((x(1)^2 - a)^2 + (x(2)^2 - a)^2 + (x(3)^2 - a)^2 + (x(4)^2 - a)^2 + (x(5)^2 - a)^2 + (x(6)^2 - a)^2 + (x(7)^2 - a)^2) + Eta * (1/abs(x(1) - x(2)) + 1/abs(x(1) - x(3)) + 1/abs(x(1) - x(4)) + 1/abs(x(1) - x(5)) + 1/abs(x(1) - x(6)) + 1/abs(x(1) - x(7)) + 1/abs(x(2) - x(3)) + 1/abs(x(2) - x(4)) + 1/abs(x(2) - x(5)) + 1/abs(x(2) - x(6)) + 1/abs(x(2) - x(7)) + 1/abs(x(3) - x(4)) + 1/abs(x(3) - x(5)) + 1/abs(x(3) - x(6)) + 1/abs(x(3) - x(7)) + 1/abs(x(4) - x(5)) + 1/abs(x(4) - x(6)) + 1/abs(x(4) - x(7)) + 1/abs(x(5) - x(6)) + 1/abs(x(5) - x(7)) + 1/abs(x(6) - x(7)));
    % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
    options = optimset('TolFun', 1e-300, 'TolX', 1e-300, 'MaxFunEvals', 10^16, 'MaxIter', 10^16);
    x_start = [(-sqrt(a)-2) (-sqrt(a)) (-sqrt(a)+1) (-sqrt(a)+2.1) (sqrt(a)-2) (sqrt(a)-1) (sqrt(a)+1)];
    %x_start = [(-sqrt(a)-2) (-sqrt(a)) (-sqrt(a)+1) 0 (sqrt(a)-2) (sqrt(a)-1) (sqrt(a)+1)];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    EqPos(i, :) = sort(x0);
    FuncVal = fval0;
end

disp('Done with 7 particle')
if FGS 
    figure(1)
    clf(figure(1))
    hold on
    title('Equilibrium positions for 7 particle')
    xlabel('\alpha')
    ylabel('\chi')
    plot(alpha, EqPos(:, 1), 'o-')
    plot(alpha, EqPos(:, 2), 'o-')
    plot(alpha, EqPos(:, 3), 'o-')
    plot(alpha, EqPos(:, 4), 'o-')
    plot(alpha, EqPos(:, 5), 'o-')
    plot(alpha, EqPos(:, 6), 'o-')
    plot(alpha, EqPos(:, 7), 'o-')
    plot(alpha, sqrt(alpha), 'k')
    plot(alpha, -sqrt(alpha), 'k')
    yline(0)
end
Eq_Pos = EqPos;
%
for q = 5:length(alpha)
    tic
    alpha(q)
    Nx1 = 15;
    Nx2 = 20;
    Nx3 = 25;
    Nx4 = 40;
    Nx5 = 25;
    Nx6 = 20;
    Nx7 = 15;

    Nx1 * Nx2 * Nx3 * Nx4 * Nx5 * Nx6 * Nx7

    a = 0.7;
    XMin1   = Eq_Pos(q, 1) - a;
    XMax1   = -Eq_Pos(q, 7) + a;
    b = 1;
    XMin2   = Eq_Pos(q, 2) - b;
    XMax2   = -Eq_Pos(q, 6) + b;
    c = 1.6;
    XMin3   = Eq_Pos(q, 3) - c;
    XMax3   = -Eq_Pos(q, 5) + c;
    d = 3;
    XMin4   = Eq_Pos(q, 4) - d;
    XMax4   = -Eq_Pos(q, 4) + d;
    e = c;
    XMin5   = Eq_Pos(q, 5) - e;
    XMax5   = -Eq_Pos(q, 3) + e;
    f = b;
    XMin6   = Eq_Pos(q, 6) - f;
    XMax6   = -Eq_Pos(q, 2) + f;
    g = a;
    XMin7   = Eq_Pos(q, 7) - g;
    XMax7   = -Eq_Pos(q, 1) + g;
    
    x1 = linspace(XMin1, XMax1, Nx1);
    x2 = linspace(XMin2, XMax2, Nx2);
    x3 = linspace(XMin3, XMax3, Nx3);
    x4 = linspace(XMin4, XMax4, Nx4);
    x5 = linspace(XMin5, XMax5, Nx5);
    x6 = linspace(XMin6, XMax6, Nx6);
    x7 = linspace(XMin7, XMax7, Nx7);
    
    dx1 = x1(2) - x1(1);
    dx2 = x2(2) - x2(1);
    dx3 = x3(2) - x3(1);
    dx4 = x4(2) - x4(1);
    dx5 = x5(2) - x5(1);
    dx6 = x6(2) - x6(1);
    dx7 = x7(2) - x7(1);
    
    K1 = -1/(2 * dx1^2);
    K2 = -1/(2 * dx2^2);
    K3 = -1/(2 * dx3^2);
    K4 = -1/(2 * dx4^2);
    K5 = -1/(2 * dx5^2);
    K6 = -1/(2 * dx6^2);
    K7 = -1/(2 * dx7^2);
    
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

    for i = 1:Nx6
        if i == 1
            KineticMtx6(i, i)        = -2 * K6;
        else
            KineticMtx6(i - 1, i)   = 1 * K6;
            KineticMtx6(i, i - 1)   = 1 * K6;
        end
    end

    for i = 1:Nx7
        if i == 1
            KineticMtx7(i, i)        = -2 * K7;
        else
            KineticMtx7(i - 1, i)   = 1 * K7;
            KineticMtx7(i, i - 1)   = 1 * K7;
        end
    end
    KineticMtx1 = sparse(KineticMtx1);
    KineticMtx2 = sparse(KineticMtx2);
    KineticMtx3 = sparse(KineticMtx3);
    KineticMtx4 = sparse(KineticMtx4);
    KineticMtx5 = sparse(KineticMtx5);
    KineticMtx6 = sparse(KineticMtx6);
    KineticMtx7 = sparse(KineticMtx7);

    K1          = kron(kron(kron(kron(kron(kron(KineticMtx1, speye(Nx2)), speye(Nx3)), speye(Nx4)), speye(Nx5)), speye(Nx6)), speye(Nx7));
    K2          = kron(kron(kron(kron(kron(kron(speye(Nx1), KineticMtx2), speye(Nx3)), speye(Nx4)), speye(Nx5)), speye(Nx6)), speye(Nx7));
    K3          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), KineticMtx3), speye(Nx4)), speye(Nx5)), speye(Nx6)), speye(Nx7));
    K4          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), KineticMtx4), speye(Nx5)), speye(Nx6)), speye(Nx7));
    K5          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), speye(Nx4)), KineticMtx5), speye(Nx6)), speye(Nx7));
    K6          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), speye(Nx4)), speye(Nx5)), KineticMtx6), speye(Nx7));
    K7          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), speye(Nx4)), speye(Nx5)), speye(Nx6)), KineticMtx7);
    KineticMtx  = K1 + K2 + K3 + K4 + K5 + K6 + K7;
    clear K1 K2 K3 K4 K5 K6 K7
    PotentialMtx1 = sparse(zeros(Nx1));
    PotentialMtx2 = sparse(zeros(Nx2));
    PotentialMtx3 = sparse(zeros(Nx3));
    PotentialMtx4 = sparse(zeros(Nx4));
    PotentialMtx5 = sparse(zeros(Nx5));
    PotentialMtx6 = sparse(zeros(Nx6));
    PotentialMtx7 = sparse(zeros(Nx7));

    U1 = 0.25 * (x1.^2 - alpha(q)).^2;
    U2 = 0.25 * (x2.^2 - alpha(q)).^2;
    U3 = 0.25 * (x3.^2 - alpha(q)).^2;
    U4 = 0.25 * (x4.^2 - alpha(q)).^2;
    U5 = 0.25 * (x5.^2 - alpha(q)).^2;
    U6 = 0.25 * (x6.^2 - alpha(q)).^2;
    U7 = 0.25 * (x7.^2 - alpha(q)).^2;

    PotentialMtx1 = sparse(diag(U1));
    PotentialMtx2 = sparse(diag(U2));
    PotentialMtx3 = sparse(diag(U3));
    PotentialMtx4 = sparse(diag(U4));
    PotentialMtx5 = sparse(diag(U5));
    PotentialMtx6 = sparse(diag(U6));
    PotentialMtx7 = sparse(diag(U7));
    
    U1          = kron(kron(kron(kron(kron(kron(PotentialMtx1, speye(Nx2)), speye(Nx3)), speye(Nx4)), speye(Nx5)), speye(Nx6)), speye(Nx7));
    U2          = kron(kron(kron(kron(kron(kron(speye(Nx1), PotentialMtx2), speye(Nx3)), speye(Nx4)), speye(Nx5)), speye(Nx6)), speye(Nx7));
    U3          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), PotentialMtx3), speye(Nx4)), speye(Nx5)), speye(Nx6)), speye(Nx7));
    U4          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), PotentialMtx4), speye(Nx5)), speye(Nx6)), speye(Nx7));
    U5          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), speye(Nx4)), PotentialMtx5), speye(Nx6)), speye(Nx7));
    U6          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), speye(Nx4)), speye(Nx5)), PotentialMtx6), speye(Nx7));
    U7          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), speye(Nx4)), speye(Nx5)), speye(Nx6)), PotentialMtx7);
    PotentialMtx = U1 + U2 + U3 + U4 + U5 + U6 + U7;
    clear U1 U2 U3 U4 U5 U6 U7
    Hamiltonian = KineticMtx + PotentialMtx;
    clear PotentialMtx KineticMtx
    X_matrix1   = sparse(diag(x1));
    X_matrix2   = sparse(diag(x2));
    X_matrix3   = sparse(diag(x3));
    X_matrix4   = sparse(diag(x4));
    X_matrix5   = sparse(diag(x5));
    X_matrix6   = sparse(diag(x6));
    X_matrix7   = sparse(diag(x7));

    X1          = kron(kron(kron(kron(kron(kron(X_matrix1, speye(Nx2)), speye(Nx3)), speye(Nx4)), speye(Nx5)), speye(Nx6)), speye(Nx7));
    X2          = kron(kron(kron(kron(kron(kron(speye(Nx1), X_matrix2), speye(Nx3)), speye(Nx4)), speye(Nx5)), speye(Nx6)), speye(Nx7));
    X3          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), X_matrix3), speye(Nx4)), speye(Nx5)), speye(Nx6)), speye(Nx7));
    X4          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), X_matrix4), speye(Nx5)), speye(Nx6)), speye(Nx7));
    X5          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), speye(Nx4)), X_matrix5), speye(Nx6)), speye(Nx7));
    X6          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), speye(Nx4)), speye(Nx5)), X_matrix6), speye(Nx7));
    X7          = kron(kron(kron(kron(kron(kron(speye(Nx1), speye(Nx2)), speye(Nx3)), speye(Nx4)), speye(Nx5)), speye(Nx6)), X_matrix7);

    clear X_matrix1 X_matrix2 X_matrix3 X_matrix4 X_matrix5 X_matrix6 X_matrix7
    X12 = diag(X1) - diag(X2);
    X13 = diag(X1) - diag(X3);
    X14 = diag(X1) - diag(X4);
    X15 = diag(X1) - diag(X5);
    X16 = diag(X1) - diag(X6);
    X17 = diag(X1) - diag(X7);
    X23 = diag(X2) - diag(X3);
    X24 = diag(X2) - diag(X4);
    X25 = diag(X2) - diag(X5);
    X26 = diag(X2) - diag(X6);
    X27 = diag(X2) - diag(X7);
    X34 = diag(X3) - diag(X4);
    X35 = diag(X3) - diag(X5);
    X36 = diag(X3) - diag(X6);
    X37 = diag(X3) - diag(X7);
    X45 = diag(X4) - diag(X5);
    X46 = diag(X4) - diag(X6);
    X47 = diag(X4) - diag(X7);
    X56 = diag(X5) - diag(X6);
    X57 = diag(X5) - diag(X7);
    X67 = diag(X6) - diag(X7);
    clear X1 X2 X3 X4 X5 X6 X7
    InteractionMtx = Eta ./ sqrt(X12.^2 + Beta^2) + Eta ./ sqrt(X13.^2 + Beta^2) + Eta ./ sqrt(X14.^2 + Beta^2) + Eta ./ sqrt(X15.^2 + Beta^2) + Eta ./ sqrt(X16.^2 + Beta^2) + Eta ./ sqrt(X17.^2 + Beta^2) + Eta ./ sqrt(X23.^2 + Beta^2) + Eta ./ sqrt(X24.^2 + Beta^2) + Eta ./ sqrt(X25.^2 + Beta^2) + Eta ./ sqrt(X26.^2 + Beta^2) + Eta ./ sqrt(X27.^2 + Beta^2) + Eta ./ sqrt(X34.^2 + Beta^2) + Eta ./ sqrt(X35.^2 + Beta^2) + Eta ./ sqrt(X36.^2 + Beta^2) + Eta ./ sqrt(X37.^2 + Beta^2) + Eta ./ sqrt(X45.^2 + Beta^2) + Eta ./ sqrt(X46.^2 + Beta^2) + Eta ./ sqrt(X47.^2 + Beta^2) + Eta ./ sqrt(X56.^2 + Beta^2) + Eta ./ sqrt(X57.^2 + Beta^2) + Eta ./ sqrt(X67.^2 + Beta^2);
    InteractionMtx = sparse(1:(Nx1 * Nx2 * Nx3 * Nx4 * Nx5 * Nx6 * Nx7), 1:(Nx1 * Nx2 * Nx3 * Nx4 * Nx5 * Nx6 * Nx7), InteractionMtx);
    
    clear X13 X14 X15 X16 X17 X21 X22 X24 X25 X26 X27 X35 X36 X37 X46 X57
    Hamiltonian = Hamiltonian + InteractionMtx;
    clear InteractionMtx
    Hamiltonian = (Hamiltonian + Hamiltonian') * 0.5;
    
    Temp_Proj = zeros((Nx1 * Nx2 * Nx3 * Nx4 * Nx5 * Nx6 * Nx7), 3);
    
    k = 1;
    for j = 1:(Nx1 * Nx2 * Nx3 * Nx4 * Nx5 * Nx6 * Nx7)
        if (X12(j) < -eps) && (X23(j) < -eps) &&  (X34(j) < -eps) && (X45(j) < -eps) && (X56(j) < -eps) && (X67(j) < -eps)
            Temp_Proj(k, :) = [k, j, 1];
            k = k + 1;
        end
    end
    clear X12 X23 X34 X45 X56 X67
    Temp_Proj = Temp_Proj(1:(k-1), :);
    Projector = sparse(Temp_Proj(:, 1), Temp_Proj(:, 2), Temp_Proj(:, 3), k - 1, (Nx1 * Nx2 * Nx3 * Nx4 * Nx5 * Nx6 * Nx7));
    toc

    tic
    Hamiltonian = Projector * Hamiltonian * Projector';
    clear Projector
    [Psi, E]    = eigs(Hamiltonian, 2, 'sa');
    Spectrum    = diag(E);
    alpha(q)
    dE(q)       = Spectrum(2) - Spectrum(1)
    toc
end
%%
toc

%%
pa = [7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5];
ps = [1.77 1.57 1.33 1.068 0.76 0.48 0.23 0.0756 0.0275 0.01548 0.015 0.01803];
figure(3)
clf(figure(3))
hold on
ylabel('Split')
xlabel('\alpha')
plot(alpha(1:5), dE, 'o')
%plot(pa, ps, '.-')
set(gca, 'Yscale', 'log')
hold off

%%
figure(5)
clf(figure(5))
hold on
plot(x1, ones(length(x1), 1), '.-')
plot(x2, 0.01 + ones(length(x2), 1), '.-')
plot(x3, 0.02 + ones(length(x3), 1), '.-')
plot(x4, 0.03 + ones(length(x4), 1), '.-')
plot(x5, 0.02 + ones(length(x5), 1), '.-')
plot(x6, 0.01 + ones(length(x6), 1), '.-')
plot(x7, 0.00 + ones(length(x7), 1), '.-')
xline([EqPos(7, :)])
xline([-EqPos(7, :)])
ylim([0.99 1.5])
hold off



data(:, 1) = alpha;
data(:, 2) = dE;

name = ['secondEDSplitting_7_particles_restricted_Nx1_' num2str(Nx1) '_Nx2_' num2str(Nx2) '_Nx3_' num2str(Nx3) '_Nx4_' num2str(Nx4) '_Nx5_' num2str(Nx5) '_Nx6_' num2str(Nx6) '_Nx7_' num2str(Nx7) '_beta_' num2str(Beta)];
save(name, 'data');

