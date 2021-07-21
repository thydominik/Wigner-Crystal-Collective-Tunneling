clc
clear all
clf(figure(1))

%1 particle trajectory calculation with simulated annealing
tic
%constants:
iteration = 1.2*10^6;
N = 2^7;                %how much points are there in the curve
r = 3;                  %the time rescaling parameter set to some arbitrary number ordo 1
eps = 10^-15;            %due to some divergencies at the -infty and infty parts some cutoff value is a must.
Q = 10;

%potential paramters:
E_0 = 0.471738;          %the energy scale of the hamiltonian
E_p = 0.471738;            %the energy scale of the potential
aa = linspace(0.5,10,Q);             %the paramter that describes the potential.
l_d = 161.07;

%curve paramters:
z = linspace(-1 + eps, 1 -eps,N);
dz = abs(z(1) -z(2));

%temperature & deviation:
T_init = 300;
T = T_init * exp(-(linspace(0,40,iteration))/2);
sigma = 0.1 * sqrt(T);
for j=1:Q
     a = aa(j);
     j
%simulted annealing part:
%first the initial position:
khi = initpos(N,a);
%then the initial action calculation:
E0 = actioncalc(khi,r,N,z,dz,a,E_p);

% figure(1)
% hold on
% %the initial curves plot:
% plot(z,khi)
% hold off


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
    if rem(i,100000) == 0
        i
        figure(2)
        clf(figure(2))
        hold on
        ylim([(-sqrt(a)-0.1) (sqrt(a)+0.1)])
        scatter(z,khi)
        plot(z,tanh(atanh(z)*r/sqrt(2)))
        xlabel(['i z(r); ', 'r = ',num2str(r),'; \alpha = ',num2str(a),'; E_p = ',num2str(E_p)],'FontSize',19)
        ylabel('\chi(z)','FontSize',19)
        hold off
        disp(['iter: ', num2str(i)])
        
    end
    
end

khii(j,:) = khi;
Last_Energy(j) = E(iteration);
end
toc
figure(3)
hold on
plot(aa,exp(-Last_Energy))
hold off
%%
% x = linspace(1,10,100);
% func = Last_Energy(1)*x.^(3/2) + 0.001415;
% figure(4)
% clf(figure(4))
% hold on
% scatter(aa,Last_Energy)
% %plot(x,func)
% xlabel('\alpha')
% ylabel('S_E')
% hold off

figure(5)
hold on
plot(z,khii(1,:))
plot(z,khii(10,:))
plot(z,khii(20,:))
plot(z,khii(30,:))
plot(z,khii(40,:))
plot(z,khii(50,:))
xlabel('z')
ylabel('\chi(z)')
hold off