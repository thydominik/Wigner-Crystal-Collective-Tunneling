%1 particle Legendre polynomial action integral calculation
clc
clear all
tic

N_Lp        = 15;       % number of polynomials
Iteration   = 10000;
div         = 400;      % number of points
r           = 0.5;        %we can get back a linear traj. with r = (alpha/2)^-(1/2)
syms x
eps     = 10^-10;    % only significant for the S_0 calculation, all first and second derivatives are non divergent integrals.
z   = linspace(-1 + 10^-10, 1 - 10^-10, div);
%z       = linspace(-1, 1, div);
alpha   = 5;

%from the analytical solution this is the exact solution for the 1 p case:
%this should give {c} = 0 at the end might be good for checking
traj = sqrt(alpha) * tanh(atanh(z)*sqrt((alpha/2)) * r );

%this will give a slightly sharper slope in the middle
traj_apr = sqrt(alpha) * tanh(atanh(z)*sqrt((alpha+0.1)/2) * r);

%It's more memory friendly if we symbolically give the expressions for the
%polynomials and store them in a matrix. (I do this with the derivatives as
%well)

for i = 1:N_Lp
    Leg_sym(i)        = (1 - x^2) * legendreP(i, x);
    Dif_Leg_sym(i)    = diff(Leg_sym(i), x);
end
Legendre        = eval(subs(Leg_sym',x,z));
Legendre_diff   = eval(subs(Dif_Leg_sym',x,z));


%guess curve (q0)
[C, q_0, dQ] = f_Guesscurve(alpha, N_Lp, div, r, z);
q_0 = traj;
%calculating the initial curves action integral.
S_0 = f_Action(q_0, r, alpha, div, z);
disp(['initial action: ' num2str(S_0)])

disp(['initial parameters: ' num2str(C)])

figure(1)
clf(figure(1))
hold on
title('initial guess curve')
plot(z, q_0, 'b')
plot(z, traj, 'r')
xlabel('z')
ylabel('q_0 (z)(blue) and exact solution(red)')
hold off

%although one could separate some part of the integrals in order to save
%some time calculating the integrals, but as of now its not that important

%iteration for loop for the C coefficients
for i = 1:Iteration
    %first derivatives -> S1
    S1 = f_int_S1(z, q_0, C, div, Legendre, Legendre_diff, N_Lp, alpha, r);
    %second derivatives -> S2
    S2 = f_int_S2(z, q_0, div, Legendre, Legendre_diff, N_Lp, alpha, r);
    
    %determining the coefficients from S1 and S2:
    C_mtx(i,:) = C;
    C = C + f_Coefficients(S1, S2, N_Lp);
    %calculating the new q_0 trajectory:
    q_0 = f_q_traj(z, C, Legendre, N_Lp, r, alpha);
    
    %calculating the action from the new q_0 traj
    Action(i) = f_Action(q_0, r, alpha, div, z);

    if mod(i,1000) == 0
        disp(Action(end) - S_0)
        
        figure(2)
        clf(figure(2))
        hold on
        plot(q_0)
        plot(traj)
        hold off
    end
end

disp('done')

figure(3)
clf(figure(3))
hold on
plot(Action)

hold off

figure(4)
clf(figure(4))
hold on
plot(C_mtx)
hold off