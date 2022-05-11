function [Psi, Spectra] = Schrodinger_VSF(Sf, VSF)
%
L = length(Sf);

Hami = zeros(L, L);
Kinetic = Hami;
Potential = Kinetic;

m = 1;
hbar = 1;

dx = abs(Sf(1) - Sf(2));
K = hbar^2 / (2 * dx^2 * m);

for i = 2:L
    Kinetic(i - 1, i)   = -K;
    Kinetic(i, i - 1)   = -K;
end

for i = 1:L
    Potential(i, i) = VSF(i);
end

Hami = Potential + Kinetic;

[Psi, Spectra] = eig(Hami);
Spectra = diag(Spectra);
end

