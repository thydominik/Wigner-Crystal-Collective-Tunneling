% 3 particle case
clc
clear all
eqpos = load('eq_pos');
eq_pos = eqpos.eqpos;

clf(figure(1))
clf(figure(2))
clf(figure(3))
clf(figure(4))
clf(figure(5))
clf(figure(6))
clf(figure(7))

r = 1;
eps = 10^-3;
N = 50;
khi0 = 1;
state = 30;
a = eq_pos(4,state);
rs = 18.813;
iter = 1000000;

z = linspace(-1+eps,1-eps,N);
dz = abs(z(1)-z(2));

p1_in = eq_pos(1,state);
p1_fi = -eq_pos(3,state);
p2_in = eq_pos(2,state);
p2_fi = -eq_pos(2,state);
p3_in = eq_pos(3,state);
p3_fi = -eq_pos(1,state);
    
T_init = 1;
T = T_init * exp(-(linspace(0,50,iter))/2);

sigma = 0.1 * sqrt(T);

[position, shift] = initpos(N,p1_in,p1_fi,p2_in,p2_fi,p3_in,p3_fi,z,r,rs,a);
E_0 = actioncalc(position,r,a,rs,khi0,N,z,dz,shift);
disp(E_0)

discarded = zeros(iter,1);
E = zeros(iter,1);

for i = 1:iter
    sig = sigma(i);
    pos_new = newstep(position,N,sig,p1_in,p1_fi,p2_in,p2_fi,p3_in,p3_fi);
    E_new = actioncalc(pos_new,r,a,rs,khi0,N,z,dz,shift);
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
        figure(4)
        clf(figure(4))
        hold on
%         position(1,:) = smooth(position(1,:),3);
%         position(2,:) = smooth(position(2,:),3);
%         position(3,:) = smooth(position(3,:),3);
        plot(z,position(2,:))
        plot(z,position(1,:))
        plot(z,position(3,:))
        hold off
        i 

        
    end   

end
positions = zeros(3,N+2);
positions(:,1) = [-3.637448784374272; -2.343826382887309; 3.216216171306582];
positions(:,N+2) = [-3.216216171306582; 2.343826382887309; 3.637448784374272];
positions(:,2:N+1) = position(:,1:N);
%%
figure(4)
clf(figure(4))
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
clf(figure(3))
figure(3)
hold on

plot3(position(1,:) - 1*min(position(1,:)),position(2,:) - 1*min(position(2,:)),position(3,:) - 1*min(position(3,:)))
quiver3(0,0,0,0.4219,0.8025,0.4319)
quiver3(0,0,0,-0.7071,0.0,0.7071)
quiver3(0,0,0,-0.5675,-0.5966,0.5675)
hold off



%%
clf(figure(7))
clf(figure(8))
clf(figure(9))
figure(7)
hold on
title('discarded steps')
hist(nonzeros(discarded),iter/10000)
hold off

figure(8)
hold on
set(gca, 'YScale', 'log')
plot(E-min(E))
hold off

figure(9)
hold on
set(gca, 'YScale', 'log')
plot(E)
hold off
