% monte carlo for determining the classical equilbirium positions
clc
clear all
close all

q_points = 1;                 %this gives the number of equilibrium positions
N        = 1;                   % # of particles
eq_pos   = zeros(q_points,N);  % ebben tárolom az összes pozit


a_init = -15;
a_fin = -7.4;
alpha = linspace(a_init,a_fin,q_points);

for Q = 1:q_points
    
    iteration = 1*10^5;
    %iteráció során exponenciálisan lecseng? potenciál döntés, hogy az
    %egyensúlyi pozíciók ne kerüljenek bármelyik oldalra.
    cc = linspace(0,50,iteration);
    c_init = 1;
    c = c_init * exp(-cc);
    c(iteration) = 0;
    
    %electron kcsh paraméter
    eta = 18.813;    %15.3149; %18.813;
    E_0 = 0.71738;    %0.71738; %0.478;
    
    %hosszskála
    Ld = 161;        %160; %131.07;
    
    %potenciált határozó paraméter
    a = alpha(Q);
    
    %h?mérséklet
    T_init = 3000;
    T = T_init * exp(-(linspace(0,100,iteration)));
    sig = 0.1 * sqrt(T);
    
    %initializing the problem
    position = initpos(N);
    E0 = energy(position,N,a,c(1),eta,E_0);
    discarded = [];
    
    for i = 1:iteration
        sigma = sig(i);
        
        [new_position, sigma] = new_pos(position,T(i),a,c(i),eta,sigma);
        E = energy(new_position,N,a,c(i),eta,E_0);
        E_diff = E0 - E;
        
        if E_diff > 0
            position = new_position;
            E0=E;
        elseif rand() < exp(E_diff/T(i))
            position = new_position;
            E0=E;
        else
            discarded(i)=i;
        end
        
        EE(i) = E0;
        xx(i,:) = position;
        
    end
    disp(num2str(Q))
    
    x_test = linspace(-6,6,1000);
    U = E_0 * (a/2 * x_test.^2 + (1/4) * x_test.^4 );



    eq_pos(Q,1:length(position)) = position;
    
end
%save('equilibrium_positions', eq_pos);

%%
    clf(figure(1))
    figure(1)
    
    hold on
    pos = position;
    %title('')
    xlabel('x [nm]')
    ylabel('V(x) [meV]')
    plot(Ld * x_test,(U - min(U)),'k','LineWidth',2)
    scatter(Ld * pos,E_0 * (a/2 *pos.^2 + (1/4) * pos.^4)-min(U) ,50,'filled','o','r')
    drawnow
    ylim([-0.1 24])
    xlim([-700 700])
    hold off
    
    %%
    figure(2)
    hold on
    plot(alpha,eq_pos(:,1))
    plot(alpha,eq_pos(:,2))
    plot(alpha,eq_pos(:,3))   
    plot(alpha,eq_pos(:,4))
    plot(alpha,eq_pos(:,5)) 
    hold off