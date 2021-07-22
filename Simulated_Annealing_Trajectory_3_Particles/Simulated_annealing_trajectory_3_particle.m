% 3 particle case
clc
clear all
eqpos = load('eq_pos');
eq_pos = eqpos.eqpos;

r = 0.9;              %time reparametrization free parameter
eps = 10^-6;        %z_time cutoff
N = 50*2;             %# of points in the curve
state = 30;         %
a = eq_pos(4,state);%a --> \alpha parameter of the potential
rs = 18.813;        %dimensionless interaction strength
iter = 2000000;     %# of iterations

z = linspace(-1+eps,1-eps,N);   %z time
dz = abs(z(1)-z(2));            %time difference for the integrations

p1_in = eq_pos(1,state);        %initial and final positions of the particles
p1_fi = -eq_pos(3,state);
p2_in = eq_pos(2,state);
p2_fi = -eq_pos(2,state);
p3_in = eq_pos(3,state);
p3_fi = -eq_pos(1,state);
    
T_init = 10;                     %starting temperature which is exponentially decreases
T = T_init * exp(-(linspace(0,40,iter))/2);

sigma = 0.1 * sqrt(T);          %new step deviance

[position, shift] = initpos(N,p1_in,p1_fi,p2_in,p2_fi,p3_in,p3_fi,z,r,rs,a);
E_0 = actioncalc(position,r,a,rs,N,z,dz,shift);
disp("E_0 = " + num2str(E_0))
disp("Shift = " + num2str(shift))

discarded = zeros(iter,1);
E = zeros(iter,1);

for i = 1:iter
    sig = sigma(i);
    pos_new = newstep(position,N,sig,p1_in,p1_fi,p2_in,p2_fi,p3_in,p3_fi);
    E_new = actioncalc(pos_new,r,a,rs,N,z,dz,shift);
    E_diff = E_0 - E_new;
    
    if E_diff > 0
        position = pos_new;
        E_0 = E_new;
    elseif rand() <= exp(E_diff/T(i))
        position = pos_new;
        E_0 = E_new;
    else
        discarded(i) = i;
    end
    E(i) = E_0;
    
    if rem(i,10000) == 0
        figure(1)
        clf(figure(1))
        hold on
        plot(z,position(1,:))
        plot(z,position(2,:))
        plot(z,position(3,:))
        hold off
        disp("iter= " + num2str(i) + "   "+ "E_0= " + num2str(E_0))
    end   

end

figure(2)
clf(figure(2))
hold on
z = linspace(-1,1,N);
plot(z,position(2,:),'r','LineWidth',2)
plot(z,position(1,:),'b','LineWidth',2)
plot(z,position(3,:),'k','LineWidth',2)
ylabel('\chi_i(z)','FontSize', 22)
xlabel('z','FontSize', 22)
ylim([-3.2 3.2])
hold off
disp('Done')

%%
figure(3)
clf(figure(3))
hold on
position = f_spline_fit(position, z);
plot3(position(1,:) - 1*min(position(1,:)),position(2,:) - 1*min(position(2,:)),position(3,:) - 1*min(position(3,:)))
%quiver3(0,0,0,0.4219,0.8025,0.4319)
%quiver3(0,0,0,-0.7071,0.0,0.7071)
%quiver3(0,0,0,-0.5675,-0.5966,0.5675)
hold off

figure(4)
clf(figure(4))
hold on
title('discarded steps')
hist(nonzeros(discarded),iter/10000)
hold off

figure(5)
clf(figure(5))
hold on
set(gca, 'YScale', 'log')
plot(E-min(E))
hold off

figure(6)
clf(figure(6))
hold on
set(gca, 'YScale', 'log')
plot(E)
hold off
