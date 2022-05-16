function [gof, fitted_curve] = f_fitting_VS_3(S, VS)
%F_FITTING_VS_3
%data:
    x = S;
    y = VS;

    starting_points = [0];

    fitfun  = fittype( @(b, x) 0.25 * b * ((x - min(S)).^2) .* ((x + min(S)).^2) + max(y) - (0.25 * b * min(S)^4));

    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

    opts.StartPoint = starting_points;

    [fitted_curve, gof] = fit(x', y', fitfun, opts);

    coeffvals = coeffvalues(fitted_curve);
    disp('VS illesztés adatai: ')
    disp(num2str(coeffvals))
end

