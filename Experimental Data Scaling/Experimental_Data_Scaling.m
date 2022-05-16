clc
clear all

% THIS WILL CHANGE: Adding the folders with the data to path
addpath("D:\BME PhD\Wigner Crystal Collective Tunneling\Data\Shahal Experimental Data")
%addpath("D:\BME PhD\Temp\#Collective tunneling\Data\Pascu Experimental Data");
addpath("D:\BME PhD\Wigner Crystal Collective Tunneling\Data");
% Loading experimental data (from P1.tsv ect.)
P1 = load('P1.tsv');
P3 = load('P3.tsv');
P5 = load('P5.tsv');
P7 = load('P7.tsv');

% Loading the experimental elbow fit data
CLP3 = load('3p_classics.tsv');

% We have ED data for 1 and 3 particles, this is the basis of the scalings
ED3 = load('E_Schrodinger_3e_eta_20.00_beta_0.01_N_100.dat');
ED1 = load('splitting energy from Schrödinger.mat');
ED1 = ED1.dE;

%% Classical data
% The critical classical alpha values fro 1,3,5 and 7 particles:
alpha_c_1 = 0;
alpha_c_3 = 4.45;
alpha_c_5 = 7.81;
alpha_c_7 = 10.61;
a_c = [alpha_c_1 alpha_c_3 alpha_c_5 alpha_c_7];

% Arc length parametrized, effective potential approximated with a 4th
% order potential and fitted as such. alpha_eff = a * (alpha - alpha_c);
% eff_beta = b
%the a scaling values:
alpha_eff_a(1) = 1;
alpha_eff_a(2) = 1.37; %1.401;
alpha_eff_a(3) = 1.552;
%alpha_eff_a(4) = ???

%effective beta values:
beta_eff_b(1) = 1;
beta_eff_b(2) = 1.5; %1.454;
beta_eff_b(3) = 1.5;
%beta_eff_b(4) = ???

%the scaling of alphas -> r = a/b^(2/3)
r13 = alpha_eff_a(2)/beta_eff_b(2)^(2/3);

%% Scaling the 1 particle ED with the 3 particle ED

figure(1)
clf(figure(1))
hold on
plot(ED1(:, 1), ED1(:, 2), '-s', 'DisplayName', '1P ED')
E_d0 = 1.5;
p = 1.3;
%plot((abs(ED3(:, 1)) - a_c(2)), abs(ED3(:, 2) - ED3(:, 3)), 'd-', 'DisplayName', '3P ED \alpha = \alpha - \alpha_c')
%plot((abs(ED3(:, 1)) - a_c(2)) * r13, abs(ED3(:, 2) - ED3(:, 3)), 'o-', 'DisplayName', '3P ED \alpha = (\alpha - \alpha_c) * a * b^{-(2/3)}')
plot((abs(ED3(:, 1)) - a_c(2)) * r13, abs(ED3(:, 2) - ED3(:, 3)) / E_d0, '.-', 'DisplayName', '3P ED \alpha = (\alpha - \alpha_c) * a * b^{-(2/3)} & \Delta = \Delta / E_{d0}')
%plot((abs(ED3(:, 1)) - a_c(2)), abs(ED3(:, 2) - ED3(:, 3)) / E_d0, '*-', 'DisplayName', '3P ED \alpha = (\alpha - \alpha_c) & \Delta = \Delta / E_{d0}')
%plot((abs(ED3(:, 1)) - 5.21) * p, abs(ED3(:, 2) - ED3(:, 3)), '*-')
text(7, 0.4, 'Best fit from these: green stars')
title('Scalings of 1 and 3 particle EDs together')
xlabel('alpha')
ylabel('\Delta')
legend
grid on
xlim([0 5])
% set(gca, 'Yscale', 'log')
hold off

%% Scaling the 1 particle Ed and the 1 particle experimental data together
alpha = abs(ED1(:, 1));
Delta = ED1(:, 2);

Voltage_1 = P1(:, 1);
dE_1 = P1(:, 2);

%scaling variables:
% shift in the volage:
V0 = 69.0;
% strecth in the voltage:
Vsc = 18.0;
% Vertical stretch in dE:
E0 = 1/ 0.030;

figure(2)
clf(figure(2))
hold on
plot(ED1(:, 1), ED1(:, 2), '-s', 'DisplayName', '1P ED')
plot((Voltage_1 - V0)/Vsc, dE_1 / E0, 'o-', 'DisplayName', '1P exp. scaled')
title('Scalings of 1 particle ED and 1 partical experiment')
xlabel('alpha')
ylabel('\Delta')
legend
grid on
% xlim([0 5])
% set(gca, 'Yscale', 'log')
hold off


% Now if the experimental data is scalable, then we would be able to scale
% the 3 particle Ed with the 3 particle experiment using {V0, Vsc, E0} and
% {a_c, D, r}
%% Scaling the 3 particle ED with the 3 particle Experiment

alpha = abs(ED3(:, 1));
Delta = ED3(:, 3) - ED3(:, 2);

Voltage_3 = P3(:, 1);
dE_3 = P3(:, 2);

%scaling variables:
% shift in the volage:
V0 = 69.0;
% strecth in the voltage:
Vsc = 18.0;
% Vertical stretch in dE:
E0 = 1/0.030;

% VoltageScaled3particle
Vs3 = ((Voltage_3 - V0) / Vsc - a_c(2));
%Vs3 = Voltage_3 - V0 - a_c(2)*Vsc;
% dEScaled3particle:
dEs3 = dE_3 / E0;
dEs3 = dEs3 / E_d0;

figure(3)
clf(figure(3))
hold on
%plot(ED1(:, 1), ED1(:, 2), '-s', 'DisplayName', '1P ED')
plot((alpha - a_c(2)), Delta / E_d0, '*-', 'DisplayName', '3P ED \alpha = (\alpha - \alpha_c) & \Delta = \Delta / E_{d0}')

plot(Vs3, dEs3, 'o-', 'DisplayName', '3P exp. scaled')

% x = 220;
% y = 45;
% z = 15;
plot(Vs3  / 3.6, dEs3 * 1.65, 'o-', 'DisplayName', '3P exp scaled w/ new factors')
%plot((Voltage_1 - V0)/Vsc, dE_1 * E0, 'o-', 'DisplayName', '1P exp. scaled')

title('Scalings of 3 particle ED and 3 partical experiment')
xlabel('alpha')
ylabel('\Delta')
legend
grid on
% xlim([0 5])
% set(gca, 'Yscale', 'log')
hold off




%%
figure(3)
clf(figure(3))
hold on
plot(ED1(:, 1), ED1(:, 2))
%Vsc = 2.2247;
Vsc = 65;
%V0 = abs((2.95 * Vsc) - 2.81776 * 10^2);
V0 = .03;
V00 = 20;
%plot((CLP3(:, 1) - Vsc)/V00 , CLP3(:, 2) * V0, 'o-')

plot((V1 - Vsc)/V00, D1 * V0, 'o')
D_2 = 6;
plot(((V2 - Vsc)/V00-13) * r , D2 * V0 * D_2)
%alpha = (abs(alpha) - 4.45);% * 1.401 / (1.454^(2/3));

%plot(alpha, Delta, '.-')

%x = linspace(0, 4, 100);
%a = 1.46;
%b = 0.495;
%func = a - b * x;
%plot(x, func)
%yline([0 10^-1 10^-2])
%xline([3 a/b])
title('x_0 = 2.96')

hold off

disp('Done.')
%%
figure(2)
clf(figure(2))
hold on



plot(alpha, Delta, '.-', 'DisplayName', 'DMRG')

x = linspace(1.2, 3.2, 100);
a = 1.46;
b = 0.495;
func = a - b * x;
plot(x, func)
yline([0 10^-1 10^-2])
xline([3 a/b])
legend
title('x_0 = 2.96')
Vsc = 2.2247;
Vsc = -220;
V00 = abs((2.95 * Vsc) - 2.81776 * 10^2);
V0 = .0145;
V00 = 1;
slope = 0.309224;
V00 = a * (slope * V0);
Vsc = CLP3(2, 1) - 2.95/V00;
plot((CLP3(:, 1) - Vsc) * V00 , CLP3(:, 2) * V0, 'o-')

plot((V2 - Vsc) * V00, D2 * V0, 'o')


hold off