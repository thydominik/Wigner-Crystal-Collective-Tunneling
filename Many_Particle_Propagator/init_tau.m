function [ModVec1, ModVec2, EigVal] = init_tau(trajectory, a, v)

x0  = trajectory(1,1);
y0  = trajectory(2,1);
z0  = trajectory(3,1);

eta = 18.813;

xxV =  (a + 3*x0^2);
yyV =  (a + 3*y0^2);
zzV =  (a + 3*z0^2);

xxU = eta*2*(1/((y0 - x0)^3) + 1/((z0 - x0)^3));
yyU = eta*2*(1/((z0 - y0)^3) + 1/((y0 - x0)^3));
zzU = eta*2*(1/((z0 - x0)^3) + 1/((z0 - y0)^3));

xyU = -eta*(2/((y0 - x0)^3));
xzU = -eta*(2/((z0 - x0)^3));
yzU = -eta*(2/((z0 - y0)^3));

mtx =  [(xxV + xxU) xyU xzU;...
        xyU (yyV + yyU) yzU;...
        xzU yzU (zzV + zzU)]

                          
[ModVec1, EigVal] = eig(mtx);

disp('Eigenvectors (columns)')
disp(num2str(ModVec1))
disp('Frequency Eigevalues (dimles)')
disp(num2str(((EigVal))))
%
% temp_vec = [v(1,1); v(2,1); v(3,1)];
% temp_norm = norm(temp_vec);
% % 
% % %a = temp_vec/temp_norm;
% % % v = cross(ModVec(:,1),a);
% % % c = ModVec(:,1)' * a;
% % % R = eye(3) + [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] + [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] * [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] * (1/(1+c));
% % %ModVec = R*ModVec;
% % 
% % 
% ModVec2(:,1) = temp_vec / temp_norm;
% ModVec2(:,2) = [ -ModVec2(2,1) ModVec2(1,1) 0];
% ModVec2(:,2) = ModVec2(:,2)/norm(ModVec2(:,2));
% ModVec2(:,3) = cross(ModVec2(:,1), ModVec2(:,2));
% disp('Eigenvectors (columns)')
% disp(num2str(ModVec2))
% % 
% % %now I will rotate the tau2 and tau3 vectors
% % 
%  B =  mtx;
%  
%  T = [ModVec2(:,2) ModVec2(:,3)];
% 
%  B_tilde = B * T;
% % 
%  B_tilde = T' * B_tilde;
% % 
%  [B_ModVec, B_EigVal] = eig(B_tilde);
%  D = B_EigVal
% % P = B_ModVec
% % 
% % T1 = T * P
% % T2 = inv(P)*T';
% % 
% % T2 * B * T1 
% % % pause
% % 
% % ModVec2(:,2) = T1(:,1)/norm(T1(:,1));
% % ModVec2(:,3) = T1(:,2)/norm(T1(:,2));
% % disp('Eigenvectors (columns)')
% % disp(num2str(ModVec2))
%  EigVal(2:3,2:3) = D ;
%  %pause
 ModVec2 = ModVec1;
 end

