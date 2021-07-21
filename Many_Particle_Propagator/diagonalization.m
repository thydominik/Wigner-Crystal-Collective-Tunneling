function [P,D] = diagonalization(B)
[vec, val] = eig(B);
D       = val;

P       = vec;
P_inv   = inv(vec);
disp('B')
disp(num2str(sqrtm(B)))
disp('pDP')
disp( num2str(P*D*P_inv))


end

