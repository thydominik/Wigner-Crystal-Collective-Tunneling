clear all
clc

N = 10000;      %number of points in the trajectory
N_p = 3;        %Number of particles

%Here we will need the alpha values and the associated eq. positions.
state = 41;     %this gets a particular alpha & eq_pos values
%later we can iterate through some of these
equilibrium_positions   = load('eq_pos');
alpha                   = equilibrium_positions.eqpos(4,state);
eq_position             = equilibrium_positions.eqpos(1:3,state);

%dimensionless interaction strengh
eta = 18.813 ;

%the energy shift due to the interaction
E_0 = f_initial_energy(N_p, alpha, eq_position, eta);
disp(E_0)

%time variables
eps         = 10^-3;
r           = 3;
z_reduced   = linspace(-1 + eps, 1 - eps, N);
z           = linspace(-1, 1, N);
y           = atanh(z_reduced) * r;
%y           = linspace(-10, 10, N);
%trajectory
%chi = zeros(N_p, N);

%the boundary conditions are the eq_positions

%the 0th step will be a guess at first, then I will try to find the
%solution using interval halving methods.

GuessValue = 10^-15;

%THIS WILL HAVE TO CHANGE------------------------------------------
%For now I wont bother with the other derivatives

chi_temp = zeros(N_p,1);
chi(:,1) = eq_position;

% for i = 1:N_p
%     chi(i, 1) = chi(i, 1) + GuessValue;
% end


%THIS WILL HAVE TO CHaNGE-------------------------------------------

%0th step:
dy = y(2) - y(1);
chi_temp = chi(:,1)' + dy * f_calculation(chi(:,1), N_p, alpha, eta, E_0) + GuessValue;

for i = 2:N-1
    
    dy = y(i) - y(i - 1);
    dy2 = dy/2;
    
    chi(:, i) = chi(:, i - 1)' + dy2 * (f_calculation(chi(:, i - 1), N_p, alpha, eta, E_0) + f_calculation(chi_temp, N_p, alpha, eta, E_0));

    chi_temp = chi(:, i)' + dy * f_calculation(chi(:, i), N_p, alpha, eta, E_0);
end

disp('Done')

figure(1)
clf(figure(1))
hold on
plot(y(1:end-1),real(chi(1,:)))
plot(y(1:end-1),chi(2,:))
plot(y(1:end-1),chi(3,:))
ylim([-5 5])
xlim([y(1) y(end)])
hold off










