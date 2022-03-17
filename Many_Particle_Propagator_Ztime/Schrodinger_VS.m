function [Spectra] = Schrodinger_VS(VS, S, a, b, c)
%Schrodinger_VS: The goal is to calculate an effectively one particle
%Schrödinger numerical eq. solver for the potnetial that we get from the
%arc length parametrized potential.

% Parameters:
    % VS - The arc length parametrized effectively 1 dimensional potential
    % S - The arc length
    % a,b and c - fitting parameters for the 4 th order potential

% -------------------------------------------------------------------------

% Output:
    % Spectra -  Energy spectrum from the original VS
% -------------------------------------------------------------------------

% Calculation paramteres:
N           = 3000;             % Spatial Resolution
Hamiltonian = sparse(N, N);     % Hamiltonian
m           = 1;                % Particles mass
hbar        = 1;                % Hbar constant

%Potentail conditioning:
bound = max(S) + max(S)/2;      % The bound tells how overextended should the potential be above and below the minima

x   = linspace(-bound, bound, N); % New Spatial variable for the calc.
dx  = x(2) - x(1);

% Hamiltonian matrix elements:
Kinetic         = sparse(N, N);

Kin_factor = hbar^2 / (2 * dx^2 * m);
Kinetic(1, 1)   = 2 * Kin_factor;
for i = 2:N
    Kinetic(i - 1, i)   = -1 * Kin_factor;
    Kinetic(i, i - 1)   = -1 * Kin_factor;
    Kinetic(i, i)       = 2 * Kin_factor;
end

Potential = sparse(N, N);

for i = 1:N
    Potential(i, i) =a + b * (x(i)^2 - c)^2;
end

Hamiltonian = full(Kinetic + Potential);

% Exact diagonalization:
[Psi, Spectra] = eig(Hamiltonian);

plot(diag(Potential))
end

