%% 1 particle

clc
clear all

N = 500;
alpha = linspace(0, 10, N);

Eq_Pos = [];

for i = 1:N
    a = alpha(i);
    Potential = @(x) 0.25 * (x^2 - a)^2;
    
    % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
    
    options = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'MaxFunEvals', 10^8, 'MaxIter', 10^8);
    x_start = 0;
    [x0, fval0] = fminsearch(Potential, x_start, options);
    Eq_Pos(i) = x0;
    Error(i)    = fval0;
end

disp('Done with 1 particle')

figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 1 particle')
xlabel('\alpha')
ylabel('\chi')
plot(alpha, Eq_Pos, '.-')
plot(alpha, sqrt(alpha), 'o')
hold off

figure(2)
clf(figure(2))
hold on
title('Error of minima')
xlabel('\alpha')
ylabel('Error')
plot(alpha, Error)
set(gca, 'Yscale', 'log')
hold off

DataStruct.NumberOfParticles = 1;
DataStruct.NumberOfAlphaPoints = N;
DataStruct.Alpha = alpha;
DataStruct.EquilibriumPositions = Eq_Pos;
DataStruct.FunctionValue = Error;


save('OneParticleEquilibriumPositions', 'DataStruct');
%% 3 particles

clc
clear all

N = 1000;
alpha = linspace(4, 14, N);
eta = 20;

Eq_pos_3 = [];

for i = 1:N
    a = alpha(i);
    Potential = @(x) 0.25 * (x(1)^2 - a)^2 + 0.25 * (x(2)^2 - a)^2 + 0.25 * (x(3)^2 - a)^2 + eta/abs(x(1) - x(2)) + eta/abs(x(1) - x(3)) + eta/abs(x(2) - x(3));
    % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
   options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10^6, 'MaxIter', 10^6);
    x_start = [-sqrt(a)-1 -sqrt(a)+1 sqrt(a)];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    Eq_pos_3(i, :) = x0;
    FuncVal(i) = fval0;
end

disp('Done with 3 particle')

figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 3 particle')
xlabel('\alpha')
ylabel('\chi')
plot(alpha, Eq_pos_3(:, 1), '.-', 'DisplayName', '\chi_1')
plot(alpha, Eq_pos_3(:, 2), '.-', 'DisplayName', '\chi_2')
plot(alpha, Eq_pos_3(:, 3), '.-', 'DisplayName', '\chi_3')
plot(alpha, sqrt(alpha), 'k', 'DisplayName', 'sqrt(\alpha)')
plot(alpha, -sqrt(alpha), 'k', 'DisplayName', '-sqrt(\alpha)')
grid
legend
hold off

for k = 1:length(alpha)
    for i = 1:3
        for j = (i + 1):3
            distances(k, i, j) = Eq_pos_3(k, i) - Eq_pos_3(k, j);
        end
    end
end

figure(3)
clf(figure(3))
hold on
plot(alpha, distances(:, 1, 2))
plot(alpha, distances(:, 2, 3))
hold off

DataStruct.NumberOfParticles = 3;
DataStruct.NumberOfAlphaPoints = N;
DataStruct.Alpha = alpha;
DataStruct.EquilibriumPositions = Eq_pos_3;
DataStruct.FunctionValue = FuncVal;


save('ThreeParticleEquilibriumPositions', 'DataStruct');
%% 5 Particles

clc
clear all

N = 1000;
alpha = linspace(7.8, 7.82, N);
eta = 20;

Eq_pos_5 = [];

for i = 1:N
    a = alpha(i);
    Potential = @(x) 0.25 * ((x(1)^2 - a)^2 + (x(2)^2 - a)^2 + (x(3)^2 - a)^2 + (x(4)^2 - a)^2 + (x(5)^2 - a)^2) + eta * (1/abs(x(1) - x(2)) + 1/abs(x(1) - x(3)) + 1/abs(x(1) - x(4)) + 1/abs(x(1) - x(5)) + 1/abs(x(2) - x(3)) + 1/abs(x(2) - x(4)) + 1/abs(x(2) - x(5)) + 1/abs(x(3) - x(4)) + 1/abs(x(3) - x(5)) + 1/abs(x(4) - x(5)));
    % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
    options = optimset('TolFun', 1e-25, 'TolX', 1e-25, 'MaxFunEvals', 10^14, 'MaxIter', 10^14);
    x_start = [-sqrt(a)-1 -sqrt(a)+0.5 -sqrt(a)+1 sqrt(a)-1 sqrt(a)+1];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    Eq_pos_5(i, :) = sort(x0);
    FuncVal(i) = fval0;
end

disp('Done with 5 particle')

figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 5 particle')
xlabel('\alpha')
ylabel('\chi')
plot(alpha, Eq_pos_5(:, 1), '.-', 'DisplayName', '\chi_1')
plot(alpha, Eq_pos_5(:, 2), '.-', 'DisplayName', '\chi_2')
plot(alpha, Eq_pos_5(:, 3), '.-', 'DisplayName', '\chi_3')
plot(alpha, Eq_pos_5(:, 4), '.-', 'DisplayName', '\chi_4')
plot(alpha, Eq_pos_5(:, 5), '.-', 'DisplayName', '\chi_5')
plot(alpha, sqrt(alpha), 'k', 'DisplayName', 'sqrt(\alpha)')
plot(alpha, -sqrt(alpha), 'k', 'DisplayName', '-sqrt(\alpha)')
grid
legend
hold off

for k = 1:length(alpha)
    for i = 1:5
        for j = (i + 1):5
            distances(k, i, j) = Eq_pos_5(k, i) - Eq_pos_5(k, j);
        end
    end
end

figure(3)
clf(figure(3))
hold on
plot(alpha, distances(:, 1, 2))
%plot(alpha, distances(:, 1, 3))
%plot(alpha, distances(:, 1, 4))
%plot(alpha, distances(:, 1, 5))
plot(alpha, distances(:, 2, 3))
%plot(alpha, distances(:, 2, 4))
%plot(alpha, distances(:, 2, 5))
plot(alpha, distances(:, 3, 4))
%plot(alpha, distances(:, 3, 5))
plot(alpha, distances(:, 4, 5))

hold off

DataStruct.NumberOfParticles = 5;
DataStruct.NumberOfAlphaPoints = N;
DataStruct.Alpha = alpha;
DataStruct.EquilibriumPositions = Eq_pos_5;
DataStruct.FunctionValue = FuncVal;


save('FiveParticleEquilibriumPositions', 'DataStruct');

%% 7 particles

clc
clear all

N = 50;
alpha = linspace(4, 28, N);
eta = 20;

Eq_pos_5 = [];

for i = 1:N
    a = alpha(i);
    Potential = @(x) 0.25 * ((x(1)^2 - a)^2 + (x(2)^2 - a)^2 + (x(3)^2 - a)^2 + (x(4)^2 - a)^2 + (x(5)^2 - a)^2 + (x(6)^2 - a)^2 + (x(7)^2 - a)^2) + eta * (1/abs(x(1) - x(2)) + 1/abs(x(1) - x(3)) + 1/abs(x(1) - x(4)) + 1/abs(x(1) - x(5)) + 1/abs(x(1) - x(6)) + 1/abs(x(1) - x(7)) + 1/abs(x(2) - x(3)) + 1/abs(x(2) - x(4)) + 1/abs(x(2) - x(5)) + 1/abs(x(2) - x(6)) + 1/abs(x(2) - x(7)) + 1/abs(x(3) - x(4)) + 1/abs(x(3) - x(5)) + 1/abs(x(3) - x(6)) + 1/abs(x(3) - x(7)) + 1/abs(x(4) - x(5)) + 1/abs(x(4) - x(6)) + 1/abs(x(4) - x(7)) + 1/abs(x(5) - x(6)) + 1/abs(x(5) - x(7)) + 1/abs(x(6) - (7)));
    % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
    options = optimset('TolFun', 1e-25, 'TolX', 1e-25, 'MaxFunEvals', 10^14, 'MaxIter', 10^14);
    x_start = [(-sqrt(a)-2) (-sqrt(a)) (-sqrt(a)+1) (-sqrt(a)+2) (sqrt(a)-2) (sqrt(a)-1) (sqrt(a)+1)];
    x_start = [(-sqrt(a)-2) (-sqrt(a)) (-sqrt(a)+1) 0 (sqrt(a)-2) (sqrt(a)-1) (sqrt(a)+1)];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    Eq_pos_5(i, :) = sort(x0);
    FuncVal = fval0;
end

disp('Done with 7 particle')

figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 5 particle')
xlabel('\alpha')
ylabel('\chi')
plot(alpha, Eq_pos_5(:, 1), '.-')
plot(alpha, Eq_pos_5(:, 2), '.-')
plot(alpha, Eq_pos_5(:, 3), '.-')
plot(alpha, Eq_pos_5(:, 4), '.-')
plot(alpha, Eq_pos_5(:, 5), '.-')
plot(alpha, Eq_pos_5(:, 6), '.-')
plot(alpha, Eq_pos_5(:, 7), '.-')
plot(alpha, sqrt(alpha), 'k')
plot(alpha, -sqrt(alpha), 'k')
yline(0)
hold off