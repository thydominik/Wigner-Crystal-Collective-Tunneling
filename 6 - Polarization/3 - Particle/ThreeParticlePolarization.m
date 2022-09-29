clc
clear all


Nx      = 100;
Alpha   = 3:0.2:12;
Kappa   = -0.2:0.01:0.2;
Eta     = 20;
Eq_Pos  = [];
Beta    = 10^-5;
Polarization = struct();
Eq_pos_3 = [];

XMin    = -7;
XMax    = 7;


x = linspace(XMin, XMax, Nx);
dx = x(2) - x(1);

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

K1      = kron(kron(KineticMtx, speye(Nx)), speye(Nx));
K2      = kron(speye(Nx), kron(KineticMtx, speye(Nx)));
K3      = kron(speye(Nx), kron(speye(Nx), KineticMtx));
KineticMtx = K1 + K2 + K3;

for alphaInd = 1:length(Alpha)
    for kappaInd = 1:length(Kappa)
        PotentialMtx = sparse(zeros(Nx));

        U = 0.25 * (x.^2 - Alpha(alphaInd)).^2 + Kappa(kappaInd) * x;
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

        InteractionMtx = Eta ./ sqrt(X12.^2 + Beta^2) + Eta ./ sqrt(X13.^2 + Beta^2) + Eta ./ sqrt(X23.^2 + Beta^2);
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
        X1          = Projector * X1 * Projector';
        X2          = Projector * X2 * Projector';
        X3          = Projector * X3 * Projector';

        [Psi, E]    = eigs(Hamiltonian, 2, 'sa');
        Spectrum    = diag(E);
        dE       = Spectrum(2) - Spectrum(1);
        
        Polarization.Splittings(alphaInd, kappaInd) = dE;

        Int1 = Psi(:, 1)' * (X1 * Psi(:, 1));
        Int2 = Psi(:, 1)' * (X2 * Psi(:, 1));
        Int3 = Psi(:, 1)' * (X3 * Psi(:, 1));

        P = Int1 + Int2 + Int3;

        WaveFunction    = Projector' * Psi(:, 1);
        WFT             = reshape(WaveFunction, [Nx Nx Nx]);
        
        Polarization.Polarization(alphaInd, kappaInd) = P;

        for Ind1 = 2:Nx
            PsiX(Ind1) = 0;
            PsiY(Ind1) = 0;
            PsiZ(Ind1) = 0;
            for Ind2 = 1:Nx
                for Ind3 = 1: Nx
                    PsiX(Ind1) = PsiX(Ind1) + WFT(Ind1, Ind2, Ind3);
                    PsiY(Ind1) = PsiY(Ind1) + WFT(Ind2, Ind1, Ind3);
                    PsiZ(Ind1) = PsiZ(Ind1) + WFT(Ind2, Ind3, Ind1);
                end
            end
        end

        PsiX = PsiX ./norm(PsiX);
        PsiY = PsiY ./norm(PsiY);
        PsiZ = PsiZ ./norm(PsiZ);
        
        NameX = ['X' '_' num2str(alphaInd) '_' num2str(kappaInd) ];
        NameY = ['Y' '_' num2str(alphaInd) '_' num2str(kappaInd) ];
        NameZ = ['Z' '_' num2str(alphaInd) '_' num2str(kappaInd) ];

        save(NameX, 'PsiX')
        save(NameY, 'PsiY')
        save(NameZ, 'PsiZ')

        Polarization.PsiX = PsiX;
        Polarization.PsiY = PsiY;
        Polarization.PsiZ = PsiZ;

        figure(5)
        clf(figure(5))
        hold on
        plot(x, 10^-3 * (0.25 * (x.^2 - Alpha(alphaInd)).^2 + Kappa(kappaInd) * x))
        plot(x, abs(PsiX))
        plot(x, abs(PsiY))
        plot(x, abs(PsiZ))
        hold off

        figure(6)
        clf(figure(6))
        hold on
        surf(x, x, abs(reshape(WFT(:, :, 10), [Nx Nx])), 'EdgeColor', 'none', 'FaceColor', 'interp')
        axis square
        hold off
        Alpha(alphaInd)
        Kappa(kappaInd)

    end
end
Polarization.Alpha              = Alpha;
Polarization.Kappa              = Kappa;
%Polarization.PolarizationMtx    = Pol;
Polarization.Info               = '3 particle polarization and data structure.';

save('Polarization', 'Polarization')







