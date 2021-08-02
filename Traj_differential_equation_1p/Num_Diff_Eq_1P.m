%Here I will try to obtain a fairly precise solution numerically for the 1
%particle instanton, then move to the 3,5, ...N particle case if this one
%works well :)
clear all
clc

%Heun method for the equation: 
%chi'[z] = (1/sqrt(2)) * r/(1-z^2) * (chi[z]^2 - alpha)

%number of points in the curve
N = 10000;

%constants:
alpha = 3;  eps = 10^-10;   r = 3;

%imaginary time
z_reduced = linspace(-1 + eps, 1 - eps, N);
z = linspace(0, 1, N);

%time step:
dz = z(2) - z(1);
dz2= dz/2;

%the trajectory:
chi = zeros(1, N);

%boundary conditions:
chi(1)          = 0;
%chi((N/2))      = 0;
chi(end)        = sqrt(alpha);

%0th step:
chi_temp = chi(1) + dz * (r/sqrt(2))/(1-z(1)^2) * (chi(1)^2 - alpha);  %bc the diff.eq. is divergent at Z = 1

for i = 2:N
    %the z(i) and z(i+1) prefactors:
    prefactor1 = (sqrt(2)*(1 - z(i - 1)^2)/r);
    prefactor2 = (sqrt(2)*(1 - z(i)^2)/r);
    
    %calculting the next term in the trajectory
    if i == N
        chi(i) = chi(i - 1) - (dz2 * ((prefactor1^-1 * (chi(i - 1)^2 - alpha)) ));
    else
        chi(i) = chi(i - 1) - (dz2 * ((prefactor1^-1 * (chi(i - 1)^2 - alpha)) + (prefactor2^-1 * (chi_temp^2 - alpha))));
    end
    
    %the next temporary chi value
    chi_temp = chi(i) + dz * (prefactor2^-1 * (chi(i)^2 - alpha));

end
disp('done')

%the full trajectory from Z in[0, 1] curve:
chi_full    = zeros(1, 2*N);
z_full      = linspace(-1, 1, 2*N);
for i = 1:2*N
    if i < N
        chi_full(i) = -chi(end - i + 1);
    elseif i > N
        chi_full(i) = chi(i - N);
    end    
end

%the analytical function:
chi_a = sqrt(alpha) * tanh(atanh(z_full)*sqrt(alpha/2)*r);

figure(1)
clf(figure(1))
hold on

plot(z_full, chi_full)
scatter(z_full(1:(N/50):end), chi_a(1:(N/50):end))
hold off



