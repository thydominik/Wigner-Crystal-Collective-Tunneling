function [gof, fitted_curve] = f_fitting_VS_2(S, VS)
    %data:
    x = S;
    y = VS;

    %starting_points = [0 0 0];
    starting_points = [0];

    %fitfun              = fittype( @(a,b,c, x) a + b.*(x.^2 - c));
    fitfun              = fittype( @(b, x) b.*(x + min(S)).^2);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    %opts.Weights = w;
    opts.StartPoint = starting_points;
    %[fitted_curve, gof] = fit(x', y', fitfun, 'StartPoint', starting_points, opts);
    [fitted_curve, gof] = fit(x', y', fitfun, opts);

    coeffvals = coeffvalues(fitted_curve);
    disp('VS illesztés adatai: ')
    disp(num2str(coeffvals))
end

