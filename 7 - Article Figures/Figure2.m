clc
clear all

FontSize = 25;

Polarization = load('Polarization3Particle.mat')
Polarization = Polarization.Polarization;

k = 1;
figure(3)
clf(figure(3))
hold on
title('Polarization: P(\alpha, \kappa)')
xlabel('')
ylabel('')
for alphaInd = 1:26
    %surf(Polarization.Alpha(alphaInd) * ones(26, 1), Polarization.Kappa, Polarization.Polarization(alphaInd, :, :), 'FaceColor','interp', 'EdgeAlpha', 0.5)
    plot3(Polarization.Alpha(alphaInd) * ones(201, 1), Polarization.Kappa, Polarization.Polarization(alphaInd, :, :), 'Color', [alphaInd/26 0.5 (1 - alphaInd/26)], 'LineWidth', 2)
end
%surf(Polarization.Alpha, Polarization.Kappa, Polarization.Polarization.', 'FaceColor','interp', 'EdgeAlpha', 0.5)
%colorbar
xlabel('\alpha', 'FontSize', FontSize)
xlim([6 9])
ylabel('\kappa', 'FontSize', FontSize)
zlabel('P(\alpha, \kappa)', 'FontSize', FontSize)
view([90 0])
grid on
axis square
ax = gca;
ax.FontSize = 20;
% Colorbar
hold off