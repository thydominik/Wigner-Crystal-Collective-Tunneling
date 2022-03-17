

%% 1 particle
clc
clear all
alpha = linspace(0, 10, 1000);

syms x real
assume(x, 'real')
for i = 1:length(alpha)
    f   = 0.25 * (x^2 - alpha(i))^2;
    df  = diff(f, x);
    f_0 = solve(df == 0, x);
    
    disp(f)
    disp(df)
    disp(f_0)

    figure(1)
    clf(figure(1))
    hold on
    fplot(f)
    fplot(df)
    xline([eval(f_0)])
    hold off

    pause
    clc
end

%% 2 particles coming soon....
clc
clear all
alpha = linspace(0, 10, 2);

syms x y
assume(x, 'real')
assume(y, 'real')

for i = 1:length(alpha)
    f = 0.25 * (x^2 - alpha(i))^2 + 0.25 * (y^2 - alpha(i))^2 + 20 * 1/(abs(x - y));
    dfx = diff(f, x);
    dfy = diff(f, y);
    
    figure(2)
    clf(figure(2))
    hold on
    fsurf(f)
    zlim([0 70])
%     ezplot(dfx)
%     ezplot(dfy)
    
    hold off
pause
end
%%
clc
clear all

%syms x y
%func = 3 * x^2 + 2*x*y + y^2 - 4*x + 5*y;
func        = @(x) 3*x(1)^2 + 2*x(1)*x(2) + x(2)^2 - 4*x(1) + 5*x(2) + 1/(x(1) - x(2));
Interval    = [-10, 10];
[xx, fval] = fminunc(func, Interval)
