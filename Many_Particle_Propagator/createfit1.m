function [cf] = createFit1(x_data, y_data, in, out)
% data
x = x_data';
y = y_data';
% Define Start points, fit-function and fit curve
x0 = [ 2 0];
a = (out + in)/2;
b = a - in;
fitfun = fittype( @(r,d,x) a + b*tanh(atanh(x)*r + d));
[fitted_curve,gof] = fit(x,y,fitfun,'StartPoint',x0);
% Save the coeffiecient values for a,b,c and d in a vector
coeffvals = coeffvalues(fitted_curve);
disp('illesztés adatai: ')
disp(num2str(coeffvals))
cf(1) = a;
cf(2) = b;
cf(3) = coeffvals(1);
cf(4) = coeffvals(2);





