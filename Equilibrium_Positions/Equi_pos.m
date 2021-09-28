clc
clear all
close all

positions   = 12;
particles   = 3;

eq_pos      = zeros(particles, positions);

alpha_start = -4;
alpha_fin   = -15;

alpha       = linspace(alpha_start, alpha_fin,positions)

for i = 1:positions  
    disp(num2str(i))
    
    iter                        = 5*10^5;
    well_shift                  = zeros(iter,1);
%     well_shift(1:(end-iter/10)) = 10*exp(-(linspace(0,50,positions - iter/10))); 
    eta                         = 20; %18.813;
    E_0                         = 0.478;
    Ld                          = 161.07;
    
    T_init                      = 600;
    T                           = T_init * exp(-(linspace(0,80,iter))/2);
    sigma                       = 0.1 * sqrt(T);
    
    position    = initpos(particles);
    E0          = energy(position, particles, alpha(i), well_shift(1), eta, E_0);
    discarded   = [];
    
    for j = 1:iter
        [new_position, sigma_new]   = new_pos(position, T(j), alpha(i), well_shift(j), eta, sigma(j));
        E                       = energy(new_position, particles, alpha(i), well_shift(j), eta, E_0);
        E_diff                  = E0 - E;
        
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
    Energy(i)   = E0;
    eq_pos(:,i) = position;
    
end
%%
eq_pos = organizer(eq_pos);
%%

figure(1)
clf(figure(1))
hold on
ylabel('\chi_i^0','FontSize',20)
xlabel('$$\tilde{a}$$', 'Interpreter', 'LaTeX', 'FontSize', 20)
scatter(alpha, eq_pos(1,:),'k', 'LineWidth',2)
scatter(alpha, eq_pos(2,:),'r', 'LineWidth',2)
scatter(alpha, eq_pos(3,:), 'LineWidth',2)
hold off

%%
for i = 1: positions
    eqpos(:,i) = eq_pos(:,i); 
end
eqpos(4,:)      = alpha;