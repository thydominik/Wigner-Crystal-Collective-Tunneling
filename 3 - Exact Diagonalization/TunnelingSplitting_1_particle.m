clc
clear all

disp('1 particle tunneling splitting calculation.')


%% 1 particle Equilibrium positions

Na = 800;
Nx = 5000;
alpha = 0:0.3:7.2;

Eq_Pos = [];

for i = 1:length(alpha)
    a = alpha(i);
    Potential = @(x) 0.25 * (x^2 - a)^2;
    % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
    options = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'MaxFunEvals', 10^9, 'MaxIter', 10^9);
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
ylabel('\chi_0')
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
% plot(alpha, (Eq_pos - sqrt(alpha)))
set(gca, 'Yscale', 'log')
hold off

%% Creating the potnential
for alphaInd = 1 : length(alpha)
    a = alpha(alphaInd);

    S   = linspace(-6, 6, Nx);
    dx  = S(2) - S(1);
    V   = 0.25 * (S.^2 - a).^2;

%     figure(3)
%     clf(figure(3))
%     hold on
%     title('Potentail')
%     xlabel('x')
%     ylabel('V(x)')
%     plot(S, V/max(V))
%     hold off
    
    PotentialMtx = sparse(diag(V));
    for i = 2:Nx
        K = 1/(2 * dx^2);
        KineticMtx(i - 1, i) = -K;
        KineticMtx(i, i - 1) = -K;
    end
    Hamiltonian = PotentialMtx + sparse(KineticMtx);

    [Psi, E]    = eig(full(Hamiltonian));
    E           = diag(E);

    figure(4)
    clf(figure(4))
    hold on
    title(a)
%     xlim([S(1) S(end)])
    ylim([0 1])
    plot(S, 10*V/max(V))
    plot(S, Psi(:, 1)/max(Psi(:, 1)))
    plot(S, Psi(:, 2)/max(Psi(:, 2)))
    yline(0)
    
    Fr = getframe(gcf);
    [Im(:, :, 1, alphaInd), Map] = rgb2ind(Fr.cdata, 8);

    hold off

    Split(alphaInd, 1) = E(2) - E(1);
    Split(alphaInd, 2) = a;

end


%%
figure(5)
clf(figure(5))
hold on
plot(Split(2:end, 2), Split(2:end, 1))
plot(-DMRG(:, 1), DMRG(:, 2))
xline([0 0.0725 0.6])
set(gca, 'Yscale', 'log')
hold off
data(:, 1) = alpha;
data(:, 2) = Split(:, 1);

name = ['EDSplitting_1_particle_Nx_' num2str(Nx)];
save(name, 'data');