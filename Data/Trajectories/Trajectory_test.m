clc
clear all

addpath('D:\BME PhD\Wigner Crystal Collective Tunneling\Data\Trajectories\eta 20\3 Particles');

traj    = load("P_1.mat");
traj = traj.position;
eqpos   = load("EqPos_eta20_alpha_5_20.mat");
eqpos = eqpos.eqpos;
alpha = eqpos(4, :);
alpha = alpha(1);

z = linspace(-1, 1, 200);
eta = 20;

figure(1)
clf(figure(1))
hold on
plot(z, traj(1, :))
plot(z, traj(2, :))
plot(z, traj(3, :))
hold off

T = [];
for i = 1:200
    xm = traj(2, i);
    options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10^6, 'MaxIter', 10^6);
    x0 = [-3 3];
    pot = @(chi) 0.25 * (chi(1)^2 + alpha)^2 + 0.25 * (chi(2)^2 + alpha)^2 + eta * (1/abs(chi(1) - chi(2)) + 1/abs(chi(1) - xm) + 1/abs(chi(2) - xm));
    [xe, fval] = fminsearch(pot, x0, options);
    T(1, i) = xe(1);
    T(2, i) = xm;
    T(3, i) = xe(2);
    F(i) = fval;
end

%%
figure(1)
clf(figure(1))
hold on
plot(z, traj(1, :))
plot(z, traj(2, :))
plot(z, traj(3, :))
plot(z, T(2, :), 'o')
plot(z, T(1, :), 'o')
plot(z, T(3, :), 'o')
hold off
%%
f_actioncalc(traj, 1.1, alpha, eta, 3, 200, z, z(2) - z(1), 1) - f_actioncalc(T, 1.1, alpha, eta, 3, 200, z, z(2) - z(1), 1)
