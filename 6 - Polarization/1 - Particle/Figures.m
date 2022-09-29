clc
clear all

% Load Data

Polarization = load('OneParticlePolarization.mat');
Polarization = Polarization.Polarization;

disp(Polarization.Info)
k = 1;
figure(k)
clf(figure(k))
k = k + 1;
hold on
title('Polarization: P(\alpha, \kappa)')
xlabel('')
ylabel('')
surf(Polarization.Alpha, Polarization.Kappa, Polarization.PolarizationMtx.')
colorbar
xlabel('\alpha')
ylabel('\kappa')
view([140 50])
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 6;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\kappa')
ylabel('P(\alpha, \kappa)')
plot(Polarization.Kappa, Polarization.PolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'QM')
plot(Polarization.Kappa, Polarization.ClassicalPolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'Classical')
legend
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 9;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\kappa')
ylabel('P(\alpha, \kappa)')
plot(Polarization.Kappa, Polarization.PolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'QM')
plot(Polarization.Kappa, Polarization.ClassicalPolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'Classical')
legend
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 12;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\kappa')
ylabel('P(\alpha, \kappa)')
plot(Polarization.Kappa, Polarization.PolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'QM')
plot(Polarization.Kappa, Polarization.ClassicalPolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'Classical')
legend
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 15;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\kappa')
ylabel('P(\alpha, \kappa)')
plot(Polarization.Kappa, Polarization.PolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'QM')
plot(Polarization.Kappa, Polarization.ClassicalPolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'Classical')
legend
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 18;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\kappa')
ylabel('P(\alpha, \kappa)')
plot(Polarization.Kappa, Polarization.PolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'QM')
%plot(Polarization.Kappa, Polarization.ClassicalPolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'Classical')
legend
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 21;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\kappa')
ylabel('P(\alpha, \kappa)')
plot(Polarization.Kappa, Polarization.PolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'QM')
%plot(Polarization.Kappa, Polarization.ClassicalPolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'Classical')
legend
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
title(['Polarization at a specific \alpha s'])
xlabel('\kappa')
ylabel('P(\alpha, \kappa)')
plot(Polarization.Kappa, Polarization.PolarizationMtx(6, :), '.-', 'DisplayName', '\alpha = 0')
plot(Polarization.Kappa, Polarization.PolarizationMtx(9, :), '.-', 'DisplayName', '\alpha = 1')
plot(Polarization.Kappa, Polarization.PolarizationMtx(12, :), '.-', 'DisplayName', '\alpha = 2')
plot(Polarization.Kappa, Polarization.PolarizationMtx(15, :), '.-', 'DisplayName', '\alpha = 3')
plot(Polarization.Kappa, Polarization.PolarizationMtx(18, :), '.-', 'DisplayName', '\alpha = 4')
plot(Polarization.Kappa, Polarization.PolarizationMtx(21, :), '.-', 'DisplayName', '\alpha = 5')
%plot(Polarization.Kappa, Polarization.ClassicalPolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'Classical')
legend
hold off

% Calculating the LDOS

LDOS = zeros(46, 101, 200);
for alphaInd = 1:length(Polarization.Alpha)
    for kappaInd = 1:length(Polarization.Kappa)
        LDOS(alphaInd, kappaInd, :) = reshape(Polarization.GSWaveFunctionMtx(alphaInd, kappaInd, :), [1 200]) .* reshape(Polarization.GSWaveFunctionMtx(alphaInd, kappaInd, :), [1 200]);
    end
end

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 6;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\chi')
ylabel('\kappa')
colormap hot
colorbar
surf(linspace(-10, 10, 200), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [101 200]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.PolarizationMtx(alphaInd, :), Polarization.Kappa, 10 * ones(101,1), 'LineWidth', 3)
xlim([-3 3])
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 9;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\chi')
ylabel('\kappa')
colormap hot
colorbar
surf(linspace(-10, 10, 200), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [101 200]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.PolarizationMtx(alphaInd, :), Polarization.Kappa, 10 * ones(101,1), 'LineWidth', 3)
xlim([-3 3])
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 12;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\chi')
ylabel('\kappa')
colormap hot
colorbar
surf(linspace(-10, 10, 200), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [101 200]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.PolarizationMtx(alphaInd, :), Polarization.Kappa, 10 * ones(101,1), 'LineWidth', 3)
xlim([-3 3])
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 15;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\chi')
ylabel('\kappa')
colormap hot
colorbar
surf(linspace(-10, 10, 200), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [101 200]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.PolarizationMtx(alphaInd, :), Polarization.Kappa, 10 * ones(101,1), 'LineWidth', 3)
xlim([-3 3])
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 18;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\chi')
ylabel('\kappa')
colormap hot
colorbar
surf(linspace(-10, 10, 200), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [101 200]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.PolarizationMtx(alphaInd, :), Polarization.Kappa, 10 * ones(101,1), 'LineWidth', 3)
xlim([-3 3])
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 21;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\chi')
ylabel('\kappa')
colormap hot
colorbar
surf(linspace(-10, 10, 200), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [101 200]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.PolarizationMtx(alphaInd, :), Polarization.Kappa, 10 * ones(101,1), 'LineWidth', 3)
xlim([-3 3])
hold off

