clc
clear all

disp('2 particle tunneling splitting calculation.')
N = 100;
Nx1     = 70;
Nx2     = 70;
alpha   = linspace(0, 15, N);
eta     = 20;
Eq_Pos  = [];
Beta    = 0.01;

XMin1    = -1;
XMax1    = 5;

XMin2   = -5;
XMax2   = 1;

x1 = linspace(XMin1, XMax1, Nx1);
x2 = linspace(XMin2, XMax2, Nx2);
dx1 = x1(2) - x1(1);
dx2 = x2(2) - x2(1);

K1 = -1/(2 * dx1^2);
K2 = -1/(2 * dx2^2);

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

KineticMtx1 = sparse(KineticMtx1);
KineticMtx2 = sparse(KineticMtx2);

K1      = kron(KineticMtx1, speye(Nx2));
K2      = kron(speye(Nx1), KineticMtx2);
KineticMtx = K1 + K2;

for i = 1:length(alpha)
    PotentialMtx1 = sparse(zeros(Nx1));
    PotentialMtx2 = sparse(zeros(Nx2));
    
    U1 = 0.25 * (x1.^2 - alpha(i)).^2;
    U2 = 0.25 * (x2.^2 - alpha(i)).^2;

    PotentialMtx1 = sparse(diag(U1));
    PotentialMtx2 = sparse(diag(U2));

    U1 = kron(PotentialMtx1, speye(Nx2));
    U2 = kron(speye(Nx1), PotentialMtx2);
    PotentialMtx = U1 + U2;
    
    Hamiltonian = KineticMtx + PotentialMtx;
    
    X_matrix1    = sparse(diag(x1));
    X_matrix2    = sparse(diag(x2));

    X1          = kron(X_matrix1, speye(Nx2));
    X2          = kron(speye(Nx1), X_matrix2);
    
    X12 = diag(X2) - diag(X1);
    
    InteractionMtx = eta ./ sqrt(X12.^2 + Beta^2);
    InteractionMtx = sparse(1:(Nx1 * Nx2), 1:(Nx1 * Nx2), InteractionMtx);
    
    Hamiltonian = Hamiltonian + InteractionMtx;
    
    Hamiltonian = (Hamiltonian + Hamiltonian') * 0.5;
    
    Temp_Proj = zeros((Nx1 * Nx2), 3);
    
    k = 1;
    for j = 1:(Nx1 * Nx2)
        if X12(j) < -eps
            Temp_Proj(k, :) = [k, j, 1];
            k = k + 1;
        end
    end

    Temp_Proj = Temp_Proj(1:(k-1), :);
    Projector = sparse(Temp_Proj(:, 1), Temp_Proj(:, 2), Temp_Proj(:, 3), k - 1, (Nx1 * Nx2));

    Hamiltonian = Projector * Hamiltonian * Projector';

    [Psi, E]    = eigs(Hamiltonian, 2, 'sa');
    Spectrum    = diag(E);
    dE(i)       = Spectrum(2) - Spectrum(1);
end
%%
figure(3)

hold on
plot(alpha, dE)
hold off

