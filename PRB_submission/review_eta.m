clc
clear all

hbar = 6.582 * 10^-16; % [eV]
w = 0.018 /hbar
eps_0 = 55.263 * 10^1;     % [e^2 eV^-1 m^-1]
eps_r = 11.4;
m_0 = 9.11 * 10^-31;    % [kg]
m = 0.191 * m_0;
m = 0.916 * m_0;
e = 1.602 * 10^-19;

a_0 = sqrt(hbar/(m * w))
E_orb   = hbar * w
E_ee    = 1/ ( 4 * pi * eps_0 * eps_r * a_0) 
R_w     = E_ee/E_orb * 1000
a_B = hbar^2 * eps_0 * eps_r / (m)


mu = a_0/a_B
