clc
clear all

addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\1 - Equilibrium Positions\Data')

P1 = load('OneParticleEquilibriumPositions.mat');
P3 = load('ThreeParticleEquilibriumPositions.mat');
P5 = load('FiveParticleEquilibriumPositions.mat');
P7 = load('SevenParticleEquilibriumPositions.mat');

%% 1 particle
P1 = P1.DataStruct;
COM1(1, :) = P1.Alpha;
COM1(2, :) = P1.EquilibriumPositions;

figure(1)
clf(figure(1))
hold on
plot(COM1(1, :), -COM1(2, :))

hold off
%% 3 particle
P3 = P3.DataStruct;
COM3(1, :) = P3.Alpha;
for i = 1:length(P3.Alpha)
    COM3(2, i) = 1/P3.NumberOfParticles * sum(P3.EquilibriumPositions(i, :));
end
figure(2)
clf(figure(2))
hold on

plot(COM3(1, :), COM3(2, :), '.-', 'DisplayName', 'Centr of Mass 3P')
legend
grid
ylabel('Center of Mass')
xlabel('\alpha')
hold off
%% 5 particle

P5 = P5.DataStruct;
COM5(1, :) = P5.Alpha;
for i = 1:length(P5.Alpha)
    COM5(2, i) = 1/P5.NumberOfParticles * sum(P5.EquilibriumPositions(i, :));
end
figure(3)
clf(figure(3))
hold on
plot(COM5(1, :), COM5(2, :), '.-', 'DisplayName', 'Centr of Mass 5P')
legend
grid
ylabel('Center of Mass')
xlabel('\alpha')
hold off

%% 7 particle

P7 = P7.eqpos;
COM7(1, :) = P7(end, :);
for i = 1:length(P7(1, :))
    COM7(2, i) = 1/7 * sum(P7(1:end-1, i));
end
figure(4)
clf(figure(4))
hold on
plot(-COM7(1, :), COM7(2, :), '.-', 'DisplayName', 'Centr of Mass 7P')
legend
grid
ylabel('Center of Mass')
xlabel('\alpha')

hold off

%%

figure(5)
clf(figure(5))
hold on
x = linspace(0, 8, 40);
f1 = -sqrt(x);
f3 = -0.4197 * x.^(0.4571);
f5 = -0.2927 * x.^(0.415);
f7 = -0.2942 * x.^(0.3413) + 0.0679;
a1 = COM1(1, :);
a3 = COM3(1, :) - 4.45;
a5 = COM5(1, :) - 7.8;
a7 = -COM7(1, :) - 10.6;
c1 = COM1(2, :);
c3 = COM3(2, :);
c5 = COM5(2, :);
c7 = COM7(2, :);



plot(-COM7(1, :) - 10.6, COM7(2, :), 'DisplayName', '7 particle')
plot(COM5(1, :) - 7.8, COM5(2, :), 'DisplayName', '5 particle')
plot(COM3(1, :) - 4.45, COM3(2, :), 'DisplayName', '3 particle')
plot(COM1(1, :), -COM1(2, :), 'DisplayName', '1 particle')

plot(x, f1, 'o', 'DisplayName', '1 particle fit -sqrt(x)')
plot(x, f3, 'o', 'DisplayName', '3 particle fit -0.4197 * x.^{0.4571}')
plot(x, f5, 'o', 'DisplayName', '5 particle fit -0.2927 * x.^{0.415};')
plot(x, f7, 'o', 'DisplayName', '7 particle fit -0.2942 * x.^{0.3413} + 0.0679')
legend
xlim([0 10])
ylim([-3 0])
grid
hold off