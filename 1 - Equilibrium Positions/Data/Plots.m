clc
clear all

Data = load('ThreeParticleEquilibriumPositions.mat'); Data = Data.DataStruct;

figure(17)
clf(figure(18))
hold on
plot(Data.Alpha, Data.EquilibriumPositions(:, 1), 'k', 'LineWidth', 2)
plot(Data.Alpha, Data.EquilibriumPositions(:, 2), 'k', 'LineWidth', 2)
plot(Data.Alpha, Data.EquilibriumPositions(:, 3), 'k', 'LineWidth', 2)
plot(Data.Alpha, -Data.EquilibriumPositions(:, 1), 'b--', 'LineWidth', 2)
plot(Data.Alpha, -Data.EquilibriumPositions(:, 2), 'b--', 'LineWidth', 2)
plot(Data.Alpha, -Data.EquilibriumPositions(:, 3), 'b--', 'LineWidth', 2)
xlabel('\alpha', 'FontSize', 20)
ylabel('\chi_0', 'FontSize', 20)

hold off