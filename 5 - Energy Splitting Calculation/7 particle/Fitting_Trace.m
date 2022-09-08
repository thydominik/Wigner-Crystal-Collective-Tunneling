function [fitresult, gof] = Fitting_Trace(ztime, trace)
    %data:
    xData = ztime;
    yData = trace;

    % Exclude points:
    excludedPoints = (xData < -0.6) | (xData > 0);

    % Options:
    opts            = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.StartPoint = [3.96246703388466 -0.00628930817610063 0.24680401287015];
    opts.Exclude    = excludedPoints;

    % Fit type:
    ft = fittype( 'gauss1' );

    % Fitting:
    [fitresult, gof] = fit( xData, yData, ft, opts );
end

