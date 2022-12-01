clc
clear all

AlphaInd = 20;
FontSize = 25;
% Figure 1: Quartic potential with some finite linear factor, indicating the wavefunctions and polarization 
addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\6 - Polarization\3 - Particle\WaveFunction data')
Pol = load('Polarization3Particle.mat'); Pol = Pol.Polarization;
WF1 = load(['X_' num2str(AlphaInd) '_1']); WF1 = WF1.PsiX; WF1 = interp1(linspace(-7, 7, 100), WF1, linspace(-7, 7, 1000), 'spline');
WF2 = load(['Y_' num2str(AlphaInd) '_1']); WF2 = WF2.PsiY; WF2 = interp1(linspace(-7, 7, 100), WF2, linspace(-7, 7, 1000), 'spline');
WF3 = load(['Z_' num2str(AlphaInd) '_1']); WF3 = WF3.PsiZ; WF3 = interp1(linspace(-7, 7, 100), WF3, linspace(-7, 7, 1000), 'spline');

figure(1)
clf(figure(1))
hold on
xlabel('\chi', 'FontSize', 25)
yyaxis right
    
    %xlim([-6 6])
    ylim([-10 70])
    x = linspace(-7, 7, 1000);
    plot(x, ((Pol.Kappa(1) * x) + 0.012 * (x.^2 - Pol.Alpha(AlphaInd)).^2) / 0.03, 'k-', 'LineWidth', 4)
    ax = gca;
    ax.YColor = 'k';
    ax.FontSize = FontSize - 5;
    ylabel('V_{scaled}(\chi)', 'FontSize', FontSize)
    
yyaxis left
    ylabel('\Psi_j(\chi)', 'FontSize', FontSize)
    xlim([-5 5])
    ylim([-0.1 0.7])
    trshld = 10^-4;
    plot(x(abs(WF1) > trshld), abs(WF1(abs(WF1) > trshld)), '-', 'LineWidth', 2.0)
    plot(x(abs(WF2) > trshld), abs(WF2(abs(WF2) > trshld)), '-', 'LineWidth', 2.0)
    plot(x(abs(WF3) > trshld), abs(WF3(abs(WF3) > trshld)), '-', 'LineWidth', 2.0)
    ax = gca;
    %ax.YColor = 'b';
    ax.FontSize = FontSize - 5;
    ylabel('\Psi_j(\chi)', 'FontSize', FontSize)

yline(0, 'k-', 'LineWidth', 3)
axis square
hold off
%
for alphaInd = AlphaInd:AlphaInd
    for kappaInd = 1:201
        name1 = ['X_' num2str(alphaInd) '_' num2str(kappaInd)];
        name2 = ['Y_' num2str(alphaInd) '_' num2str(kappaInd)];
        name3 = ['Z_' num2str(alphaInd) '_' num2str(kappaInd)];
        PsiX = load(name1); PsiX = PsiX.PsiX;
        PsiY = load(name2); PsiY = PsiY.PsiY;
        PsiZ = load(name3); PsiZ = PsiZ.PsiZ;
        LD1 = PsiX .* PsiX;
        LD2 = PsiY .* PsiY;
        LD3 = PsiZ .* PsiZ;
        LDOS(alphaInd, kappaInd, :) = LD1 + LD2 + LD3;
    end
end

%
figure(2)
clf(figure(2))
hold on
%title(['Polarization at a specific \alpha = ' num2str(Pol.Alpha(alphaInd))])
xlabel('\chi', 'FontSize', FontSize)
ylabel('\epsilon', 'FontWeight', 'bold',  'FontSize', FontSize)
%colormap turbo
colorbar
ax = gca;
ax.FontSize = 20;
surf(linspace(-7, 7, 100), Pol.Kappa, reshape(LDOS(alphaInd, :, :), [201 100]), 'EdgeColor', 'none', 'FaceColor', 'interp')
plot3(Pol.Polarization(alphaInd, :), Pol.Kappa, 10 * ones(26,201), 'k--', 'LineWidth', 5)
%plot3(linspace(-7, 7, length(Pol.Polarization(alphaInd, :))), 0 * ones(1, length(Pol.Kappa)), 10 * ones(26,201), 'k-', 'LineWidth', 1.5)
xlim([-5.5 5.5])
xticks([-5 -2.5 0 2.5 5])
yticks([-0.1 -0.05 0 0.05 0.1])
axis square
hold off

    