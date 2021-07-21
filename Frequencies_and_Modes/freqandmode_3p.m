clc

eqpos = eq_pos';

len = 100;

x0 = 0;
y0 = 0;
z0 = 0;

ld  = 161.07;
eta = 18.813;                       

alpha   = eqpos(:,4);       
E_p     = 0.478;

figure(1)
clf(figure(1))
hold on
xlabel('\alpha','FontSize',20)
ylabel('\chi_{i}^0','FontSize',20)
plot(alpha,eqpos(:,1),'g')
plot(alpha,(eqpos(:,2)),'r')
plot(alpha,eqpos(:,3),'k')
hold off

%load('eq_pos_fitted_wider_200points');
eigenvalues = zeros(len(1),3);
eigenvectors = zeros(len(1),3,3);
for i=1:len

    x01 = eqpos(i,1);
    x02 = eqpos(i,2);
    x03 = eqpos(i,3);
    
    a = alpha(i);
%the potential derivative the first two letters denote the partial
%derivatives the third denotes the potential type.
xxV =  (a + 3*x01^2);
yyV =  (a + 3*x02^2);
zzV =  (a + 3*x03^2);

xxU = eta*2*(1/((x02 - x01)^3) + 1/((x03 - x01)^3));
yyU = eta*2*(1/((x03 - x02)^3) + 1/((x02 - x01)^3));
zzU = eta*2*(1/((x03 - x01)^3) + 1/((x03 - x02)^3));

xyU = -eta*(2/((x02 - x01)^3));
xzU = -eta*(2/((x03 - x01)^3));
yzU = -eta*(2/((x03 - x02)^3));


mtx =    [(xxV + xxU) xyU xzU;...
                    xyU (yyV + yyU) yzU;...
                    xzU yzU (zzV + zzU)];

[ModVec, EigVal] = eig(mtx);
eigenvalues(i,:) = diag(EigVal);
eigenvectors(i,:,:) = ModVec;

end

figure(2)
clf(figure(2))
hold on
scatter(alpha,sqrt((eigenvalues(:,1))),'r','Linewidth',2)
scatter(alpha,sqrt((eigenvalues(:,2))),'k','Linewidth',2)
scatter(alpha,sqrt((eigenvalues(:,3))),'b','Linewidth',2)
xlabel('$$\tilde{a}$$', 'Interpreter', 'LaTeX', 'FontSize', 16)
ylabel('$\hbar\omega [meV]$', 'Interpreter', 'latex','FontSize',16)
hold off
