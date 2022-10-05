function [C] = f_Coefficients(S1,S2, N_Lp)
%this function is not yet working
S2_new = zeros(N_Lp, N_Lp);

for i = 1:N_Lp
    for j = 1:N_Lp
        S2_new(i,j) = S2(i,j) + S2(j,i);
    end
end

C = linsolve(- S2_new, S1'); 
%C(2:2:end) = 0;
C = C';


end

