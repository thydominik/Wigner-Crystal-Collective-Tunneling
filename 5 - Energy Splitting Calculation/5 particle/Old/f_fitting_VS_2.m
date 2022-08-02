function [gof, fitted_curve] = f_fitting_VS_2(S, VS)
    %data:
    x = S;
    y = VS;

    %starting_points = [0 0 0];
    starting_points = [0 0];
%     w = linspace(-3, 3, length(S));
%     sigma = 15;
%     w = exp(-w.^2 / sigma);
%     w = 1 - w;
%     w = w / max(w);
%     w(70:end-70) = 0;

%     w = (w + 0.5);
%     w = w / max(w);
    % fitfun              = fittype( @(a,b,c, x) a + b.*(x.^2 - c));
    fitfun              = fittype( @(a,b, x) a + b.*(x + min(S)).^2);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    %opts.Weights = w;
    opts.StartPoint = starting_points;
    %[fitted_curve, gof] = fit(x', y', fitfun, 'StartPoint', starting_points, opts);
    [fitted_curve, gof] = fit(x', y', fitfun, opts);

    coeffvals = coeffvalues(fitted_curve);
    disp('VS illesztés adatai: ')
    disp(num2str(coeffvals))
end

