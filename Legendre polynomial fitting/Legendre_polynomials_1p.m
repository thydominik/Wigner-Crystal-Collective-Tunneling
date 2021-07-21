%1 particle Legendre polynomial action integral calculation

clc
clear all
tic

N           = 20;       % number of polynomials
Iteration   = 100000;
div         = 100;      % number of points
r           = 2;
syms x
eps     = 10^-10;    % only significant for the S_0 calculation, all first and second derivatives are non divergent integrals.
z_res   = linspace(-1+eps, 1-eps, div);
z       = linspace(-1,1,div);
alpha   = 10;

func = sqrt(alpha) * tanh(atanh(z)*(r + sqrt(r)) );

%It's more memory friendly if we symbolically give the expressions for the
%polynomials and store them in a matrix.

for i = 1:N
    Leg_sym(i)        = legendreP(i, x);
    Dif_Leg_sym(i)    = diff(Leg_sym(i), x);
end
Leg     = eval(subs(Leg_sym',x,z));
Dif_Leg = eval(subs(Dif_Leg_sym',x,z));

%guess curve (q0): let's give this in terms of legendre polynomials as well
[C_0, q_0] = guesscurve(alpha, N, div, Leg);
disp(['initial parameters: ' num2str(C_0)])

figure(1)
clf(figure(1))
hold on
title('initial guess curve')
plot(z, q_0)
xlabel('z')
ylabel('q_0 (z)')
hold off

[S_21, S_23] = const_integrals(Leg, Dif_Leg, N, div, z, r, alpha);

for i = 1:Iteration
    [S1, S2] = integrals(Leg, Dif_Leg, C_0, N, div, z, r, alpha, S_21, S_23);
    C_0 = (C_0 + Coeffs(S1, S2, N));
    q_0 = q_curve(C_0, Leg, N, div, alpha);

    Action(i) = action(q_0, Leg, Dif_Leg, C_0, r, alpha, div, z_res, N);

    if rem(i,1000) == 0
        disp(i)
        %disp(['Coefficients: ' num2str(C_0)])
        figure(2)
        clf(figure(2))
        hold on
        plot(z, q_0)
        plot(z, func)
        hold off
    end

end
