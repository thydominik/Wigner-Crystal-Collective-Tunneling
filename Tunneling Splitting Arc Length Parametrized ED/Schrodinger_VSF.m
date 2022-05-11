function [Psi, Spectra] = Schrodinger_VSF(Sf, VSf)
%
L = length(Sf);
%VSf = smooth(VSf);
%Sf = smooth(Sf);
Hami        = sparse(L, L);
Kinetic     = Hami;
Potential   = Kinetic;

m       = 1;
hbar    = 1;

dx  = abs(Sf(2) - Sf(1));
K   = hbar^2 / (2 * dx^2 * m);

Kinetic(1, 1) = 0;
for i = 2:L
    Kinetic(i - 1, i)   = -K;
    Kinetic(i, i - 1)   = -K;
end

for i = 1:L
    Potential(i, i) = VSf(i);
end

Hami = full(Potential + Kinetic);

[Psi, Spectra] = eig(Hami);
Spectra = diag(Spectra);

end

