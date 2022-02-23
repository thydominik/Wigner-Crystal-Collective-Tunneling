function [gof, fitted_curve] = f_fitting_VS(S, VS)
    %data:
    x = S;
    y = VS;

    starting_points = [0 0 0];

    fitfun              = fittype( @(a,b,c, x) a + b.*(x.^2 - c).^2);
    [fitted_curve, gof] = fit(x', y', fitfun, 'StartPoint', starting_points);

    coeffvals = coeffvalues(fitted_curve);
    disp('VS illesztés adatai: ')
    disp(num2str(coeffvals))
end