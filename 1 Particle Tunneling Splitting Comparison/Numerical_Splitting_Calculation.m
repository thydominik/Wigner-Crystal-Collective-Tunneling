clc
clear all
%2022.02.23
%Replicating the anaylitical calcualtions using numerical data:
%Data read:
Chi     = load('Chi.mat');
Chi     = Chi.khii;
Alpha   = load('Alpha.mat');
Alpha   = Alpha.Alpha;
Action  = load('Action.mat');
Action  = Action.Action;

%rescaling param (from the code that produces the trajectories)
    %Mind that I did recalculate the trajectories for this so this r value
    %is unique to this batch of trajectories! --> Folder name
r = 2;

%Omega:
Omega = sqrt(2 * abs(Alpha)); 

for i = 1:length(Alpha)
    %Potential:
    V = 0.25 .* (Chi(i, :).^2 - abs(Alpha(i))).^2;
    %Classical Impulse:
    P       = sqrt(2 * V);
    P_0     = sqrt(2 * 0.25 * Alpha(i)^2);
    
    %Integral one by Chi:--------------------------------------------------
    Int1 = 0;
    func1 = Omega(i)./(P) - 1./(Chi(i, :) - Chi(i, 1));
    for k = 2:length(Chi(1,:))/2-1
        dz      = Chi(i, k + 1) - Chi(i, k);
        Int1 = Int1 + dz * 0.5 * (func1(k + 1) + func1(k));
    end
    %Milnikov prefactor with the First Integral:
    Split(1, i) = 2 * Omega(i) * sqrt(Alpha(i)) * sqrt(Omega(i)/pi)  * exp(Int1) * exp(-Action(i));

    %Integral two by tau:--------------------------------------------------
    %Time:
    eps = 10^-15;  
    z   = linspace(-1 + eps, 1 -eps, length(P));
    
    %Derivative of P:
    dP(1) = (P(2) - P(1))./(Chi(i, 2) - Chi(i, 1));
    for k = 2:length(P)-1
        dP(k) = (P(k+1) - P(k-1))./(Chi(i, k + 1) - Chi(i, k - 1));
    end
    dP(length(P)) = (P(end) - P(end-1))./(Chi(i, end) - Chi(i, end-1));
    
    plot(dP)
    plot(P)
    pause
    Int2 = 0;
    func2 = (r./(1 - z.^2)) .* (Omega(i) - dP);
    for k = 2:length(Chi(1, :))/2 - 1
        Int2 = Int2 + (z(k+1) - z(k)) * 0.5 * (func2(k + 1) + func2(k));
    end


    Split(2, i) =  sqrt(4 * Omega(i) / pi) * P(end/2) * exp(-Action(i)) * exp(Int2);

    %Analytical solution:--------------------------------------------------
    Action_a  = 2/3 * sqrt(2) * Alpha(i).^(3/2);   
    Split(3, i) = sqrt(Omega(i)/pi) * 4 * Omega(i) * sqrt(Alpha(i)) * exp(-Action_a);
end




figure(1)
clf(figure(1))
hold on
title('Numerical Action vs Analyitical Action')
xlabel('\alpha')
ylabel('S_0')
plot(Alpha, Action, '.-', 'DisplayName', 'Numerical S')
plot(Alpha, 2/3 * sqrt(2) * sqrt(Alpha) .* Alpha, 'o-', 'DisplayName', 'Analytical S')
legend
hold off

figure(2)
clf(figure(2))
hold on
plot(Alpha, Action - (2/3 * sqrt(2) * sqrt(Alpha) .* Alpha), '.-' )
title('Difference between the S')

hold off

figure(3)
clf(figure(3))
hold on
title('Numerical and Analytical splits')
ylabel('\Delta')
xlabel('\alpha')
plot(Alpha, Split(1, :),'.-', 'DisplayName', 'dz integral numerical')
plot(Alpha, Split(2, :),'x-', 'DisplayName', 'd\tau integral numerical')
plot(Alpha, Split(3, :), 'o-', 'DisplayName', 'Analytical Solution')
legend
hold off

figure(4)
clf(figure(4))
hold on
title('Numerical and Analytical splits on log scale')
ylabel('\Delta')
xlabel('\alpha')
plot(Alpha, Split(1, :),'.-', 'DisplayName', 'dz integral numerical')
plot(Alpha, Split(2, :),'x-', 'DisplayName', 'd\tau integral numerical')
plot(Alpha, Split(3, :), 'o-', 'DisplayName', 'Analytical Solution')
set(gca, 'Yscale', 'log')
legend
hold off

figure(5)
clf(figure(5))
hold on
title(' Differences between the splittings')
xlabel('\alpha')
ylabel('\Delta')
plot(Alpha, 100 * (Split(3,:) ./ Split(1, :) - 1), '.-', 'DisplayName', 'Analytical - dz integral ')
plot(Alpha, 100 * (Split(3,:) ./ Split(2, :) - 1), '.-', 'DisplayName', 'Analytical - d\tau integral ')
%set(gca, 'Yscale', 'log')
legend
hold off

figure(6)
hold on
title('Splittings Numerical & Analytical Calculation')
plot(Alpha, Split(1, :), 'o-', 'DisplayName', 'Numerical Milnikov')
xlabel('-\alpha')
hold off

