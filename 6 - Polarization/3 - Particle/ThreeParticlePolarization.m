clc
clear all

Nx      = 200;
Alpha   = 10;
Kappa   = 0.1;
Eta     = 20;
Beta    = 10^-5;
Polarization = struct();

XMin    = -6;
XMax    = 6;


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
%%
        WaveFunction1    = sqrt(3) * Projector' * Psi(:, 1);
        WFT1             = reshape(WaveFunction1, [Nx Nx Nx]);

        Polarization.Polarization(alphaInd, kappaInd) = P;

        for Ind1 = 2:Nx
            rho_even(Ind1) = 0;
            PsiX1(Ind1) = 0;
            PsiY1(Ind1) = 0;
            PsiZ1(Ind1) = 0;
            for Ind2 = 1:Nx
                for Ind3 = 1: Nx
                    PsiX1(Ind1) = PsiX1(Ind1) + WFT1(Ind1, Ind2, Ind3)^1;
                    PsiY1(Ind1) = PsiY1(Ind1) + WFT1(Ind2, Ind1, Ind3)^1;
                    PsiZ1(Ind1) = PsiZ1(Ind1) + WFT1(Ind2, Ind3, Ind1)^1;

                end
            end
            rho_even(Ind1) = PsiX1(Ind1) + PsiY1(Ind1) + PsiZ1(Ind1);
        end

        PsiX1 = sqrt(PsiX1 ./norm(PsiX1));
        PsiY1 = sqrt(PsiY1 ./norm(PsiY1));
        PsiZ1 = sqrt(PsiZ1 ./norm(PsiZ1));

        WaveFunction2    = sqrt(3) * Projector' * Psi(:, 2);
        WFT2             = reshape(WaveFunction2, [Nx Nx Nx]);

        Polarization.Polarization(alphaInd, kappaInd) = P;

        for Ind1 = 2:Nx
            PsiX2(Ind1) = 0;
            PsiY2(Ind1) = 0;
            PsiZ2(Ind1) = 0;
            for Ind2 = 2:Nx
                for Ind3 = 2:Nx
                    PsiX2(Ind1) = PsiX2(Ind1) + WFT2(Ind1, Ind2, Ind3)^2;
                    PsiY2(Ind1) = PsiY2(Ind1) + WFT2(Ind2, Ind1, Ind3)^2;
                    PsiZ2(Ind1) = PsiZ2(Ind1) + WFT2(Ind2, Ind3, Ind1)^2;
                end
            end
            rho_odd(Ind1) = PsiX2(Ind1) + PsiY2(Ind1) + PsiZ2(Ind1);
        end

        PsiX2 = sqrt(PsiX2 ./norm(PsiX2));
        PsiY2 = sqrt(PsiY2 ./norm(PsiY2));
        PsiZ2 = sqrt(PsiZ2 ./norm(PsiZ2));

        NameX = ['X' '_' num2str(alphaInd) '_' num2str(kappaInd) ];
        NameY = ['Y' '_' num2str(alphaInd) '_' num2str(kappaInd) ];
        NameZ = ['Z' '_' num2str(alphaInd) '_' num2str(kappaInd) ];

        save(NameX, 'PsiX1')
        save(NameY, 'PsiY1')
        save(NameZ, 'PsiZ1')

        Polarization.PsiX = PsiX1;
        Polarization.PsiY = PsiY1;
        Polarization.PsiZ = PsiZ1;

        figure(5)
        clf(figure(5))
        hold on
        plot(x, 10^-3 * (0.25 * (x.^2 - Alpha(alphaInd)).^2 + Kappa(kappaInd) * x))
        plot(x, (PsiX2 + PsiY2 + PsiZ2), 'o')
        plot(x, (PsiX2), '.-')
        plot(x, (PsiY2), '.-')
        plot(x, (PsiZ2), '.-')
        hold off

        % figure(6)
        % clf(figure(6))
        % hold on
        % surf(x, x, abs(reshape(WFT1(:, :, 10), [Nx Nx])), 'EdgeColor', 'none', 'FaceColor', 'interp')
        % axis square
        % hold off
        % Alpha(alphaInd)
        % Kappa(kappaInd)

    end
end
Polarization.Alpha              = Alpha;
Polarization.Kappa              = Kappa;
%Polarization.PolarizationMtx    = Pol;
Polarization.Info               = '3 particle polarization and data structure.';

save('Polarization', 'Polarization')
%%
WF1 = (PsiX1 + PsiY1 + PsiZ1); %WF1 = interp1(x, WF1, linspace(min(x), max(x), 10^3), 'cubic');% WF1 = 3 * WF1 / norm(WF1);
WF2 = (PsiX2 + PsiY2 + PsiZ2);% WF2 = interp1(x, WF2, linspace(min(x), max(x), 10^3), 'cubic');% WF2 = 3*  WF2 / norm(WF2);
%xnew = linspace(min(x), max(x), 10^3);
figure(6)
clf
hold on
plot(x, WF1.^2)
plot(x, WF2.^2)

hold off

%%
figure(5)
clf(figure(5))
hold on
plot(x, 10^-3 * (0.25 * (x.^2 - Alpha(alphaInd)).^2 + Kappa(kappaInd) * x))
plot(xnew, WF2 , '-', 'LineWidth', 3)
plot(x, (PsiX2), '.-')
plot(x, (PsiY2), '.-')
plot(x, (PsiZ2), '.-')
hold off

figure(6)
clf(figure(6))
hold on
plot(x, 10^-3 * (0.25 * (x.^2 - Alpha(alphaInd)).^2 + Kappa(kappaInd) * x))
plot(xnew, WF1, '-', 'LineWidth', 3)
plot(x, (PsiX1), '.-')
plot(x, (PsiY1), '.-')
plot(x, (PsiZ1), '.-')
hold off
%%
figure(8)
clf(figure(8))
hold on
plot(xnew ,(WF1 + WF2).^2 - (WF1 - WF2).^2)
plot(xnew, WF1 - WF2)
hold off
%%

figure(10)
clf(figure(10))
hold on
plot(x, rho_even, 'DisplayName', '\rho_E')
plot(x, rho_odd, 'DisplayName', '\rho_O')
box
legend
hold off
%%
figure(11)
clf(figure(11))
hold on
plot(x, rho_even + rho_odd, 'DisplayName', '\rho_L')
plot(x, rho_even - rho_odd, 'DisplayName', '\rho_R')
legend
box
hold off
%%
figure(12)
clf(figure(12))
hold on
plot(x, rho_odd*2)
hold off

%%

for Ind1 = 2:Nx
    PsiX3(Ind1) = 0;
    PsiY3(Ind1) = 0;
    PsiZ3(Ind1) = 0;
    for Ind2 = 2:Nx
        for Ind3 = 2:Nx
            PsiX3(Ind1) = PsiX3(Ind1) + (WFT2(Ind1, Ind2, Ind3) * WFT1(Ind1, Ind2, Ind3));
            PsiY3(Ind1) = PsiY3(Ind1) + (WFT2(Ind2, Ind1, Ind3) * WFT1(Ind2, Ind1, Ind3));
            PsiZ3(Ind1) = PsiZ3(Ind1) + (WFT2(Ind2, Ind3, Ind1) * WFT1(Ind2, Ind3, Ind1));
        end
    end
    rho_mix(Ind1) = PsiX3(Ind1) + PsiY3(Ind1) + PsiZ3(Ind1);
end
rho_mix = rho_mix / norm(rho_mix);
rho_even = rho_even / norm(rho_even);
rho_odd = rho_odd / norm(rho_odd);
figure(13)
clf(figure(13))
hold on
plot(x, rho_mix, 'DisplayName', 'Mix')
plot(x, rho_even, '.', 'DisplayName', 'even')
plot(x, rho_odd, 'DisplayName', 'odd')
box
legend

hold off

%%
clc

figure(14)
clf(figure(14))
hold on

%plot(rho_even)
%plot(flip(rho_even))
plot(conv(rho_even - flip(rho_even),normpdf([-2:0.10:2],0,1.1)));
box
hold off


close all
FontSize = 30;
Position = [1 1 6 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);
hold on
C = conv(rho_even - flip(rho_even), normpdf([-2:0.10:2],0,1.1));
C = C / norm(C);

I = 0;
x = linspace(-6, 6 , length(C)) * 160;
for i = 2:length(C)
    dx = x(i) - x(i - 1);
    I = I + dx (C(i) + C(i-1))/2; 
end
I
norm(C)
%%
plot(linspace(-6, 6, length(C)) * 160,1/160 * C , 'r-', 'LineWidth', 6)

set(gca, 'FontSize', FontSize)
%xlim([-7*160 7*160])
%xticks([-500 0 500])
%ylim([-0.02 0.5])
%yticks([0 0.5])
xlabel('$z [{\rm{nm}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
%label = ylabel('$\rho(z) \left[ \frac{1}{nm}\right] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
x = linspace(-7, 7, 1000);
mult = 100;
yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
trshld = 8*10^-4;
shift1 = -9.5;
%text( 10, 0.42, 0, '$\epsilon = 10^{-4}$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_difference.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);