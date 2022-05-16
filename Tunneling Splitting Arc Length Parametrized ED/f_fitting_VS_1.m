function [gof, fitted_curve] = f_fitting_VS_1(S, VS)
    %data:
    x = S;
    y = VS;

    starting_points = [0 0 0];

    fitfun              = fittype( @(a,b,c,x) 0.5 * a * x.^2 + 0.25 * b * x.^4 + c);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.StartPoint = starting_points;
    [fitted_curve, gof] = fit(x', y', fitfun, opts);

    coeffvals = coeffvalues(fitted_curve);
    disp('VS illesztés adatai: ')
    disp(num2str(coeffvals))
end