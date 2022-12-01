clc
clear all

addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\6 - Polarization\3 - Particle\WaveFunction data')
Polarization = load('Polarization3Particle.mat'); Polarization = Polarization.Polarization
for alphaInd = 1:26
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
%%
alphaIndicies = [5 8 11 14 17 18 19 20];
k = 1;
figure(11)
clf(figure(11))
hold on
for i = alphaIndicies
    subplot(2, 4, k)
        alphaInd = i;
        hold on
        surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [201 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
        plot3(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(26,201),'k--', 'LineWidth', 3)
        xlim([-5 5])
        ax = gca;
    ax.FontSize = 14;
        %axis square
        title(['\alpha = ' num2str(Polarization.Alpha(alphaInd))], 'FontSize', 20)
        xlabel('\chi', 'FontSize', 20)
        ylabel('\epsilon', 'FontSize', 20)
        view([00 90])
        colormap turbo
        colorbar
        hold off
    k = k + 1;
    
end


