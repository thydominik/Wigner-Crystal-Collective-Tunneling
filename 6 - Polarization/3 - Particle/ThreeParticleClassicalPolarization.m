clc
clear all

disp('3 particle classical polarization calculatio')
Na = 100;
Nk = 100;
Alpha   = linspace(3, 10, Na);
Kappa   = linspace(-0.1, 0.1, Nk);

Eta     = 20; 
Beta    = 10^-5;
ClassicalPolarization = struct();

for AlphaInd = 1:length(Alpha)
    for KappaInd = 1:length(Kappa)
        a = Alpha(AlphaInd);
        k = Kappa(KappaInd);
        Potential = @(x) 0.25 * (x(1)^2 - a)^2 + (k * x(1)) + 0.25 * (x(2)^2 - a)^2 + (k * x(2)) + 0.25 * (x(3)^2 - a)^2 + (k * x(3)) + Eta/abs(x(1) - x(2) + Beta) + Eta/abs(x(1) - x(3) + Beta) + Eta/abs(x(2) - x(3) + Beta);
        options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10^6, 'MaxIter', 10^6);
        x_start = [-sqrt(a)-2 0 sqrt(a)];
        [x0, fval0] = fminsearch(Potential, x_start, options);
        Eq_pos_3(AlphaInd, KappaInd, :) = sort(x0);
        Polarization(AlphaInd, KappaInd) = sum(x0);
        FuncVal(AlphaInd, KappaInd) = fval0;
    end
end
%%
figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 3 particle')
xlabel('\alpha')
ylabel('\chi')
plot(Alpha, Eq_pos_3(:, 1, 1), '.-', 'DisplayName', '\chi_1')
plot(Alpha, Eq_pos_3(:, 1, 2), '.-', 'DisplayName', '\chi_2')
plot(Alpha, Eq_pos_3(:, 1, 3), '.-', 'DisplayName', '\chi_3')
grid
legend
hold off

%%
figure(2)
clf(figure(2))
hold on
ylabel('\kappa')
xlabel('\alpha')
title('Classical Polarization')
surf(Alpha, Kappa, Polarization.', 'FaceColor', 'interp')
hold off
