clc
clear all

%Script for calculating the classical equilibriumm positions in a quartic potential 
particles   = 7;        % # of particles in the system 

alpha_start = -9;                                           % The initial value of alpha
alpha_fin   = -21;                                          % The Final value of alpha
positions   = abs(alpha_fin - alpha_start) * 10 + 1;        % # of alpha values ( right now it's 0.1 increments)
alpha       = linspace(alpha_start, alpha_fin, positions);  % All alpha values
eq_pos      = zeros(particles, positions);                 

disp(['Positions = ' num2str(positions)])
%disp(['Time ~ ' num2str(positions * 3.05)])


for i = 1:positions             % Loop: different potentials// ifferent alphas
    disp(num2str(i))
    
    iter        = 2 * 10^6;           % # of iterations for one particular alpha value
    eta         = 20; %18.813;      % Interaction strength
    E_0         = 0.478;            % Energy unit -- Irrelevant for now
    Ld          = 161.07;           % Length Unit
    well_shift  = zeros(1, iter);
    T_init      = 600;                                      % Initial Temperature
    T           = T_init * exp(-(linspace(0,80,iter))/2);   % T(iteration)
    sigma       = 0.1 * sqrt(T);                            % New position divergence metric
    
    position    = initpos(particles);
    E0          = energy(position, particles, alpha(i), well_shift(1), eta, E_0);
    discarded   = [];       %discarded iterations
    
    for j = 1:iter  %Loop: particular alpha iterations
        [new_position, sigma_new]   = new_pos(position, T(j), alpha(i), well_shift(j), eta, sigma(j));
        E                           = energy(new_position, particles, alpha(i), well_shift(j), eta, E_0);
        E_diff                      = E0 - E;
        
        %Simulated annealing part
        if E_diff > 0
            position = new_position;
            E0=E;
        elseif rand() < exp(E_diff/T(j))
            position = new_position;
            E0=E;
        else
            discarded(i)=i;
        end
     
    end
    %Saving stuff
    Energy(i)   = E0;
    eq_pos(:,i) = position;
end

eq_pos = organizer(eq_pos, particles);
disp('Done!')


figure(1)
clf(figure(1))
hold on
ylabel('\chi_i^0','FontSize',20)
xlabel('$$\tilde{a}$$', 'Interpreter', 'LaTeX', 'FontSize', 20)
plot(alpha, eq_pos(1,:), '.', 'DisplayName', '1st p.')
plot(alpha, eq_pos(2,:), '.', 'DisplayName', '2st p.')
plot(alpha, eq_pos(3,:), '.', 'DisplayName', '3st p.')
plot(alpha, eq_pos(4,:), '.', 'DisplayName', '4st p.')
plot(alpha, eq_pos(5,:), '.', 'DisplayName', '5st p.')
plot(alpha, eq_pos(6,:), '.', 'DisplayName', '6st p.')
plot(alpha, eq_pos(7,:), '.', 'DisplayName', '7st p.')
plot(alpha, -sqrt(-alpha), 'k-', 'DisplayName', '-sqrt(\alpha)')
plot(alpha, sqrt(-alpha), 'k-', 'DisplayName', 'sqrt(\alpha)')
legend
grid on
hold off


for i = 1: positions
    eqpos(:,i) = eq_pos(:,i); 
end
eqpos(particles + 1, :) = alpha;
FileName                = ['Eq_Pos_eta_20_particles_' num2str(particles)];
save(FileName, 'eqpos')