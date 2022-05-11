clear all
clc

a = 3;
x = linspace(-5,5,1000);
func = a*(exp(4.2 * (x+5))-1) + 3;
logfunc = log(func);

figure(1)
clf(figure(1))
hold on
plot(x,func)
%set(gca,'Yscale','log')
hold off