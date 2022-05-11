%% Simulated annealing
clc
clear all
%close all
eq_pos = zeros(100,3);
feloszt = 1;
alpha = linspace(-5,-4,feloszt);
for Q = 1:feloszt
     clearvars -except position Q eq_pos alpha feloszt
%important variables
N = 3;    % electron number
%m=10;

iter = 100000;                      %100K or more is recomended
%-3.171653056150049;%-4.119;          %the paramter that controls the two well depths
cc = linspace(0,50,iter);         % temporary shift of the wells
c = 10*exp(-cc);                %exponentially decaying shift
c(iter) = 0;                
el_stre = 18.813;%15.3149 ;%18.813;%               %parameter that controls the repuslion strength
E_0 = 0.478;% %1.752;                 %energy scaling paramter
Ld = 161.07; %131.07; %160;%133;                 %scaling paramter
a = alpha(Q);%alpha*(Ld*10^9)^2 /(2*(E_0 / 1.60217662 * 10^19))


%the exponentially decreasing temperature and possible movement size
%paramter
T_init = 600;
T = T_init*exp(-(linspace(0,80,iter))/2);
sigma2 = 0.1 * sqrt(T);

%nitializing the problem and calculating its energy
position = initpos(N);
E0=energy(position,N,a,c(1),el_stre,E_0);
discarded=[];       %keeping track of the discarded steps in the annealing

q = 0;
for i=1:iter
    sigma = sigma2(i);
    
    [new_position,sigma] = new_pos(position,T(i),a,c(i),el_stre,sigma);
    E = energy(new_position,N,a,c(i),el_stre,E_0);
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
    Sig(i) = sigma;
    
     if rem(i,10000) == 0 && i ~= 0
         i
%         i
%         q = q + 1;
%         EE_mean(q) = 0;
%         for k = i-99:i
%             EE_mean(q) = EE_mean(q) + EE(k);
%         end
%         EE_mean(q) = EE_mean(q)/100;
%         xx_std(q) = std(xx(i-99:i,1));
%     end   

end
Q
%disp('done');
position;
% Plots; the annealing part has to be executed first

%%
clf(figure(2))
x_test = linspace(-6,6,1000);
U = E_0 * (a/2 * x_test.^2 + (1/4) * x_test.^4 );

%clf(figure(2))
figure(2)
%load('potential')
hold on
pos = position;
title('pozíciók')
xlabel('x poz')
ylabel('Potenciál')
%scatter(x_dd/1000, V_dd - min(V_dd))
plot(Ld * x_test,(U - min(U)),'LineWidth',2)
scatter(Ld * pos,E_0 * (a/2 *pos.^2 + (1/4) * pos.^4)-min(U) ,'filled','o')
drawnow
ylim([-0.1 24])
xlim([-600 600])
hold off


eq_pos(Q,1:length(position)) = position;


end
%%
clf(figure(1))
clf(figure(3))
clf(figure(4))
clf(figure(5))
clf(figure(6))
clf(figure(8))
clf(figure(7))

figure(1)
hold on
set(gca, 'YScale', 'log')
title({'minimum energia' min(EE),'utolsó energia', EE(length(EE))})
ylabel('Energia')
xlabel('iter szám')
plot((EE-min(EE)))
%x=linspace(1,iter,length(EE_mean));
%plot(x,(EE_mean-min(EE_mean)))
hold off



figure(3)
hold on
title('részecskék poziciója')
for i=1:N
    plot((xx(:,i)))
end
hold off

dis = nonzeros(discarded);

figure(4)
hold on
title('Elutasított lépések hisztogrammja ~1000 lépésre')
ylabel('elutasított lépések')
xlabel('iter szám')
hist(dis,iter/10000)
hold off

figure(5)
hold on
set(gca, 'YScale', 'log')
xlabel('iter szám')
ylabel('T & simga ')
plot(Sig.^2)
plot(T)
legend('új lépés norm dist. szórása','T')
hold off

for i=1:(length(EE)-1)
   E_d(i) = EE(i) - EE(i+1); 
end

figure(6)
hold on

title({'Energia különbség:  E_i - E_{i+1} ', 'utolsó érték', E_d(length(E_d))})
xlabel('iter szám')
ylabel('dE')
plot((E_d))
hold off

figure(7)
hold on
title({'temperature', min(T)})
plot(T)
hold off

figure(8)
hold on
set(gca, 'YScale', 'log')
title({'elsõ x változó szórása 100 lépésre',min(xx_std)})
plot(xx_std)
hold off

format long
for i=1:N
   position(i)*Ld
end 
end

%% A test for the distribution of each new step

uzt=linspace(0,1,length(xx));
xxx=xx(:,1);
for i=2:(length(xxx))
    uzt(i)=xxx(i)-xxx(i-1);
end
figure(9)
uzt(uzt==0)=NaN;
hist(uzt,100)

%%
