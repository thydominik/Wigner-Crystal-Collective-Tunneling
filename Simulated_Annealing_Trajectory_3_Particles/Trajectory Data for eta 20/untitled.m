clc
clear all

T = load('Pos_eta_20_alpha_8_r_1_1.mat');
T = T.position;

figure(1)
clf(figure(1))
hold on
z = linspace(-1, 1, 200);

plot(z, T(1, :), '-', 'LineWidth', 2)
plot(z, T(2, :), '-', 'LineWidth', 2)
plot(z, T(3, :), '-', 'LineWidth', 2)
xlabel('z', 'FontSize', 12)
ylabel('\chi_i', 'FontSize', 18)
% yline(0)
% xline(0)
hold off
%%
figure(2)
clf(figuer(2))
hold on

hold off