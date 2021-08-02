function [C] = f_Coefficients(S1,S2, N_Lp)
%this function is not yet working
for i = 1:N_Lp
    for j = 1:N_Lp
        if i == j
            S2(i,j) = S2(i,j) * 2;
        end
    end
end
C = linsolve(S2, -S1');
C(2:2:end) = 0;
C = C';
end

