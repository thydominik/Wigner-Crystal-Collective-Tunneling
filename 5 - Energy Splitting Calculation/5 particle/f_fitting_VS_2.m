function [gof, fitted_curve] = f_fitting_VS_2(S, VS)
    %data:
    x = S;
    y = VS;

    starting_points = [0];

    fitfun              = fittype( @(b, x) b.*(x - min(S)).^4);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.StartPoint = starting_points;
    [fitted_curve, gof] = fit(x', y', fitfun, opts);

    coeffvals = coeffvalues(fitted_curve);
    disp('VS illesztés adatai: ')
    disp(num2str(coeffvals))
end

