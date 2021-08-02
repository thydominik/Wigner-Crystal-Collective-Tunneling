%Test: How well we can fit the 1 particle analytical trajectory with
%Legendre polynomials
clear all
clc

%the analytical solution with respect to the alpha parameter
alpha   = 20;
% z = tanh(tau/r) -> tau = atanh(z)*r; 
r = 1/(sqrt(alpha/2));
%number of points
N   = 1000;
tau = linspace(-50,50,N);
z   = linspace(-1,1,N);

Trajectory_tau  = sqrt(alpha) * tanh(sqrt(alpha/2)*tau);
Trajectory_z    = sqrt(alpha) * tanh(sqrt(alpha/2)*atanh(z) * r);

%(Max)Number of Legendre polynomials:
N_Lp = 19;

%Legendre polynomial matrix:
syms x
for i = 1:N_Lp      
    Legendre_symbolical(i) = legendreP(i, x);
end
Lp_mtx = eval(subs(Legendre_symbolical',x,z));

for j = 2:N_Lp
    disp(j)
    %integration for the coefficitents
    Coeffs = zeros(1, j);
    for i = 1:j
        Lp = Lp_mtx(i, :);
        Coeffs(i) = ((2 * i) + 1)/2 * f_integration(Trajectory_z, Lp, N, z);
    end
    
    SS_tot = 0;
    SS_res = 0;
    
    fit_fun = Lp_mtx(1:j, :) .* Coeffs';
    fit_fun = sum(fit_fun,1);
    
    traj_sum = sum(Trajectory_z);
    for i = 1:N
        SS_res = SS_res + (Trajectory_z(i) - fit_fun(i))^2;
        SS_tot = SS_tot + (Trajectory_z(i) - ((1/N) * traj_sum))^2;
    end
    
    figure(2)
    clf(figure(2))
    hold on
    plot(z, Trajectory_z, 'r')
    plot(z, fit_fun, 'b')
    hold off
    
    R_sq(j) = 1 - (SS_res/SS_tot);
    R_x(j) = j;
end

figure(1)
clf(figure(1))
hold on
plot(R_x, abs(R_sq))
set(gca, 'YScale', 'log')
hold off