clc
clear all

% Load Data

Polarization = load('Polarization.mat')
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
surf(Polarization.Alpha(1:end-5), Polarization.Kappa, Polarization.Polarization(1:end-1, :).', 'FaceColor','interp')
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
plot(Polarization.Kappa, Polarization.Polarization(alphaInd, :), '.-', 'DisplayName', 'QM')
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
plot(Polarization.Kappa, Polarization.Polarization(alphaInd, :), '.-', 'DisplayName', 'QM')
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
plot(Polarization.Kappa, Polarization.Polarization(alphaInd, :), '.-', 'DisplayName', 'QM')
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
plot(Polarization.Kappa, Polarization.Polarization(alphaInd, :), '.-', 'DisplayName', 'QM')
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
plot(Polarization.Kappa, Polarization.Polarization(alphaInd, :), '.-', 'DisplayName', 'QM')
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
plot(Polarization.Kappa, Polarization.Polarization(alphaInd, :), '.-', 'DisplayName', 'QM')
legend
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 24;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\kappa')
ylabel('P(\alpha, \kappa)')
plot(Polarization.Kappa, Polarization.Polarization(alphaInd, :), '.-', 'DisplayName', 'QM')
legend
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 27;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\kappa')
ylabel('P(\alpha, \kappa)')
plot(Polarization.Kappa, Polarization.Polarization(alphaInd, :), '.-', 'DisplayName', 'QM')
legend
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
title(['Polarization at a specific \alpha s'])
xlabel('\kappa')
ylabel('P(\alpha, \kappa)')
plot(Polarization.Kappa, Polarization.Polarization(12, :), '.-', 'DisplayName', '\alpha = 5.2')
plot(Polarization.Kappa, Polarization.Polarization(15, :), '.-', 'DisplayName', '\alpha = 5.8')
plot(Polarization.Kappa, Polarization.Polarization(18, :), '.-', 'DisplayName', '\alpha = 6.4')
plot(Polarization.Kappa, Polarization.Polarization(21, :), '.-', 'DisplayName', '\alpha = 7')
plot(Polarization.Kappa, Polarization.Polarization(24, :), '.-', 'DisplayName', '\alpha = 7.6')
plot(Polarization.Kappa, Polarization.Polarization(27, :), '.-', 'DisplayName', '\alpha = 8.2')
%plot(Polarization.Kappa, Polarization.ClassicalPolarizationMtx(alphaInd, :), '.-', 'DisplayName', 'Classical')
legend
hold off

% Calculating the LDOS

LDOS = zeros(41, 41, 100);
for alphaInd = 1:41
    for kappaInd = 1:41
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
surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [41 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(41,1), 'LineWidth', 3)
xlim([-5 5])
axis square
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
surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [41 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(41,1), 'LineWidth', 3)
xlim([-5 5])
axis square
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
surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [41 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(41,1), 'LineWidth', 3)
xlim([-5 5])
axis square
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
surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [41 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(41,1), 'LineWidth', 3)
xlim([-5 5])
axis square
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
surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [41 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(41,1), 'LineWidth', 3)
xlim([-5 5])
axis square
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
surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [41 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(41,1), 'LineWidth', 3)
xlim([-5 5])
axis square
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 22;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\chi')
ylabel('\kappa')
colormap hot
colorbar
surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [41 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(41,1), 'LineWidth', 3)
xlim([-5 5])
axis square
hold off

figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 23;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\chi')
ylabel('\kappa')
colormap hot
colorbar
surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [41 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(41,1), 'LineWidth', 3)
xlim([-5 5])
axis square
hold off
figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 24;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\chi')
ylabel('\kappa')
colormap hot
colorbar
surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [41 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(41,1), 'LineWidth', 3)
xlim([-5 5])
axis square
hold off
figure(k)
clf(figure(k))
k = k + 1;
hold on
alphaInd = 25;
title(['Polarization at a specific \alpha = ' num2str(Polarization.Alpha(alphaInd))])
xlabel('\chi')
ylabel('\kappa')
colormap hot
colorbar
surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [41 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
line(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(41,1), 'LineWidth', 3)
xlim([-5 5])
axis square
hold off

