clc
clear all
clf(figure(1))

%1 particle trajectory calculation with simulated annealing

%constants:
iteration = 10^6;
N = 2^7;                %how much points are there in the curve
r = 2;                  %the time rescaling parameter set to some arbitrary number ordo 1
eps = 10^-10;            %due to some divergencies at the -infty and infty parts some cutoff value is a must.

%potential paramters:
E_0 = 0.71738;          %the energy scale of the hamiltonian
E_p = 0.71738;            %the energy scale of the potential
a = 1;             %the paramter that describes the potential.
l_d = 161.07;

%curve paramters:
z = linspace(-1 + eps, 1 -eps,N);
dz = abs(z(1) -z(2));

%temperature & deviation:
T_init = 100;
T = T_init * exp(-(linspace(0,40,iteration))/2);
sigma = 0.1 * sqrt(T);

%simulted annealing part:
%first the initial position:
khi = initpos(N,a);
%then the initial action calculation:
E0 = actioncalc(khi,r,N,z,dz,a);

figure(1)
hold on
%the initial curves plot:
plot(z,khi)
hold off


for i= 1:iteration
    khi_new = newstep(khi,N,sigma(i),a);
    E_new = actioncalc(khi_new,r,N,z,dz,a,E_p);
    E_diff = E0 - E_new;
    
    if E_diff > 0
        khi = khi_new;
        E0 = E_new;
    elseif rand() <= exp(E_diff/T(i))
        khi = khi_new;
        E0 = E_new;
    else
        discarded(i) = i;
    end
    E(i) = E0;
    if rem(i,10000) == 0
        %disp(i)
        disp(E(i))
        figure(2)
        clf(figure(2))
        hold on
        ylim([(-sqrt(a)-0.1) (sqrt(a)+0.1)])
        scatter(z,khi)
        plot(z,sqrt(a)*tanh(atanh(z)*r*sqrt(a)/sqrt(2)))
        xlabel(['i z(r); ', 'r = ',num2str(r),'; \alpha = ',num2str(a),'; E_p = ',num2str(E_p)],'FontSize',19)
        ylabel('\chi(z)','FontSize',19)
        hold off
        disp(['iter: ', num2str(i)])
        
    end
    
end

%%

figure(2)
        clf(figure(2))
        hold on
        ylim([(-sqrt(a)-0.1) (sqrt(a)+0.1)])
        scatter(z,khi,'k')
        plot(z,2*tanh(atanh(z)*6/sqrt(2)),'r','LineWidth',2)
        xlabel(['i z(r); ', 'r = ',num2str(r),'; \alpha = ',num2str(a),'; E_p = ',num2str(E_p)],'FontSize',19)
        xlabel('z','FontSize',25)
        ylabel('\chi(z)','FontSize',25)
        hold off
        disp(['iter: ', num2str(i)])
%%
figure(5)
hold on
plot(z(1:length(z)-1),diff(khi))
plot(z(1:length(z)-2),diff(diff(khi)))
hold off

%%
xx = linspace(-1,1,2^8);
func = tanh(atanh(xx)*r/sqrt(2));

figure(6)
clf(figure(6))
hold on
plot(xx(1:length(xx)-1),diff(func))
plot(xx(1:length(xx)-2),diff(diff(func)))
hold off