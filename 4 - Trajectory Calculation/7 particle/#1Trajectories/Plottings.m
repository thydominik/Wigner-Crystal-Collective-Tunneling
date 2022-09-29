clc
clear all

T = load('Traj_5p_10.mat');
T = T.Position;
z = linspace(-1, 1, 160);

figure(1)
clf(figure(1))
hold on
plot(z, T, 'LineWidth', 2)
xline(0)
yline(0)
xlabel('z @ \alpha = 13.5', 'FontSize', 20)
ylabel('\chi_i', 'FontSize', 20)
hold off