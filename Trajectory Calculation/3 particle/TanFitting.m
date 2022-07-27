function [GOFMTX] = TanFitting(Pos, NoP, PN, z)
    for particleInd = 1:PN
        y = Pos(particleInd, 2:end-1).';
        x = z(2:end-1).';
        
        if particleInd == floor(PN/2) + 1
            starting_points = [1];
            b = abs(y(1));
            fitfun = fittype( @(c, x) b*tanh(atanh(x)* c));
            [fitted_curve, gof] = fit(x, y, fitfun, 'StartPoint', starting_points);
            coeffvals = coeffvalues(fitted_curve);
            disp('illesztés adatai: ')
            disp(num2str(coeffvals))
        else
            starting_points = [1 2 3 4];

            fitfun = fittype( @(a, b, c, d, x) a + b.*tanh(atanh(x).* c + d));
            [fitted_curve, gof] = fit(x, y, fitfun, 'StartPoint', starting_points);
            coeffvals = coeffvalues(fitted_curve);
            disp('illesztés adatai: ')
            disp(num2str(coeffvals))
        end
        GOFMTX(particleInd, :) = coeffvals;
    end
end

