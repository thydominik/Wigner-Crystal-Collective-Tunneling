    clc
clear all
N = 350;
x = linspace(-3, 3, N);
a = 0;
b = 0;
c = 10;

V = [100 zeros(1, 49) a*ones(1, 50) zeros(1, 50) c*ones(1, 50) zeros(1, 50) b*ones(1, 50) zeros(1, 49) 100];
%V = 7.447 * (x.^2 - 0.61).^2;

Hami        = zeros(N, N);
Kinetic     = zeros(N, N);
Potential   = zeros(N, N);

K = 1/(2 * (x(2) - x(1))^2);

for i = 2:N
    Kinetic(i - 1, i) = -K;
    Kinetic(i, i - 1) = -K;
end

for i = 1:N
    Potential(i, i) = V(i);
end

Hami = Potential + Kinetic;
[Psi, Spectra] = eig(Hami);
Spectra = diag(Spectra);

figure(1)
clf(figure(1))
hold on
fact = 1;
plot(x, V/max(V), '.-')
plot(x, abs(fact * Psi(:, 1)), 'o-')
plot(x, abs(fact * Psi(:, 2)))
% plot(x, fact * Psi(:, 4).^2)
% plot(x, fact * Psi(:, 5).^2)
%ylim([0 1000])
hold off

