%% 7 particles

clc
clear all

N = 500;
alpha = linspace(10.5, 10.7, N)
%alpha = 8:0.1:15;
N = length(alpha);
eta = 20;

Eq_pos_7 = [];

for i = 1:N
    a = alpha(i);
    Potential = @(x) 0.25 * ((x(1)^2 - a)^2 + (x(2)^2 - a)^2 + (x(3)^2 - a)^2 + (x(4)^2 - a)^2 + (x(5)^2 - a)^2 + (x(6)^2 - a)^2 + (x(7)^2 - a)^2) + eta * (1/abs(x(1) - x(2)) + 1/abs(x(1) - x(3)) + 1/abs(x(1) - x(4)) + 1/abs(x(1) - x(5)) + 1/abs(x(1) - x(6)) + 1/abs(x(1) - x(7)) + 1/abs(x(2) - x(3)) + 1/abs(x(2) - x(4)) + 1/abs(x(2) - x(5)) + 1/abs(x(2) - x(6)) + 1/abs(x(2) - x(7)) + 1/abs(x(3) - x(4)) + 1/abs(x(3) - x(5)) + 1/abs(x(3) - x(6)) + 1/abs(x(3) - x(7)) + 1/abs(x(4) - x(5)) + 1/abs(x(4) - x(6)) + 1/abs(x(4) - x(7)) + 1/abs(x(5) - x(6)) + 1/abs(x(5) - x(7)) + 1/abs(x(6) - x(7)));
    % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
    options = optimset('TolFun', 1e-30, 'TolX', 1e-30, 'MaxFunEvals', 10^16, 'MaxIter', 10^16);
    x_start = [(-sqrt(a)-2) (-sqrt(a)) (-sqrt(a)+1) (-sqrt(a)+2.1) (sqrt(a)-2) (sqrt(a)-1) (sqrt(a)+1)];
    %x_start = [(-sqrt(a)-2) (-sqrt(a)) (-sqrt(a)+1) 0 (sqrt(a)-2) (sqrt(a)-1) (sqrt(a)+1)];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    Eq_pos_7(i, :) = sort(x0);
    FuncVal = fval0;
end

disp('Done with 7 particle')



data = load('SevenParticleEquilibriumPositions.mat');
figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 7 particle')
xlabel('\alpha')
ylabel('\chi')
% plot(-data.eqpos(8, :), data.eqpos(1, :), '.-')
% plot(-data.eqpos(8, :), data.eqpos(2, :), '.-')
% plot(-data.eqpos(8, :), data.eqpos(3, :), '.-')
% plot(-data.eqpos(8, :), data.eqpos(4, :), '.-')
% plot(-data.eqpos(8, :), data.eqpos(5, :), '.-')
% plot(-data.eqpos(8, :), data.eqpos(6, :), '.-')
% plot(-data.eqpos(8, :), data.eqpos(7, :), '.-')

plot(alpha, Eq_pos_7(:, 1), '.-')
plot(alpha, Eq_pos_7(:, 2), '.-')
plot(alpha, Eq_pos_7(:, 3), '.-')
plot(alpha, Eq_pos_7(:, 4), '.-')
plot(alpha, Eq_pos_7(:, 5), '.-')
plot(alpha, Eq_pos_7(:, 6), '.-')
plot(alpha, Eq_pos_7(:, 7), '.-')
plot(alpha, sqrt(alpha), 'k')
plot(alpha, -sqrt(alpha), 'k')
yline(0)
hold off