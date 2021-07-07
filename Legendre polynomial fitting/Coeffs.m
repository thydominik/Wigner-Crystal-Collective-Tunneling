function [C] = Coeffs(S1,S2, N)
%     for i = 1:N
%         C(i) = (-S1(i))/S2(i);
%     end
C = linsolve(S2,-S1);
C(2:2:end) = 0;
C = C';
end

