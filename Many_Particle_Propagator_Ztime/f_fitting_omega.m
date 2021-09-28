function [gof, fitted_curve] = f_fitting_omega(x_data, y_data)
    %data
    x = x_data;
    y = y_data;

    a = y_data(1);
    starting_points = [0 0 0 0];

    fitfun = fittype( @(b,c,d,e,x) a + b.*(x+1) + c.*(x+1).^2 + d.*(x+1).^3 + e.*(x+1).^4);
    [fitted_curve, gof] = fit(x, y, fitfun, 'StartPoint', starting_points);

    coeffvals = coeffvalues(fitted_curve);
    disp('illesztés adatai: ')
    disp(num2str(coeffvals))
end



