clear all 
close all  


xMin = -6; % xmin in units of ld
xMax = 6; %xmax in units of ld;
beta = 0.6;

%parameters for the run. 

alpha =0.0; 
eta = 5.0; 
num=20; 

%for alpha = 0:-0.1:-10

    fname = sprintf('E_Schrodinger_3e_eta_%.2f_N_%d.dat', eta, num);

file=fopen(fname, 'a');

clear K_matrix U_matrix H_matrix E

x=linspace(xMin, xMax, num);
dx = (x(2)-x(1));
dx2 = dx^2;

%second derivative matrix corresponding to kinetic energy; 

off = ones(num-3,1);
K_matrix = sparse(-2*eye(num-2) + diag(off,1) + diag(off,-1))/dx2;

K1_matrix = kron(kron(K_matrix, speye(num-2)),speye(num-2));
K2_matrix = kron(kron(speye(num-2), K_matrix), speye(num-2));
K3_matrix = kron(kron(speye(num-2), speye(num-2)),K_matrix);
%potential energy matrix; 

K_matrix = K1_matrix+K2_matrix+K3_matrix;

clear K1_matrix K2_matrix K3_matrix;

U= 0.25*x.^4 + 0.5*alpha*x.^2;
U = U-min(U);
U_matrix =sparse(diag(U(2:end-1)));


U1_matrix = kron(kron(U_matrix, speye(num-2)), speye(num-2)); %potential matrix; 
U2_matrix = kron(kron(speye(num-2), U_matrix), speye(num-2));
U3_matrix = kron(kron(speye(num-2), speye(num-2)), U_matrix);

U_matrix = U1_matrix+U2_matrix+U3_matrix;

clear U1_matrix U2_matrix U3_matrix;

%Coulomb energy
x_matrix = sparse(zeros(num-2));

for j =1:(num-2)
    x_matrix(j,j) = x(j+1);
end

KU_matrix = -1/2*K_matrix + U_matrix;

clear K_matrix U_matrix;

x1_matrix = kron(kron(x_matrix, speye(num-2)), speye(num-2));
x2_matrix = kron(kron(speye(num-2), x_matrix), speye(num-2));
x3_matrix = kron(kron(speye(num-2), speye(num-2)), x_matrix);



fprintf('Generating UC_matrix...\n');

x12 = diag(x1_matrix)-diag(x2_matrix);
x13 = diag(x1_matrix)-diag(x3_matrix);
x23 = diag(x2_matrix)-diag(x3_matrix);

x12=sqrt(x12.^2+beta^2);
x13=sqrt(x13.^2+beta^2);
x23=sqrt(x23.^2+beta^2);


UC_matrix = diag(sparse(eta./x12) + sparse(eta./x13) + sparse(eta./x23));

%Hamiltonian matrix;
%clear x1_matrix x2_matrix x3_matrix;

H_matrix = KU_matrix + UC_matrix;

%clear KU_matrix UC_matrix;

H_matrix = (H_matrix +H_matrix')/2;
% Find Eignevalues E_n and Eigenfunctions psi_N
%[~, e_values] = eig(H_matrix);
%the first N states; 
N=10;
fprintf('Diagonalizing...\n')
e_values = eigs(H_matrix, N, 'sm');

e_values = sort(e_values);


for j=1:N
E(j)=e_values(j);
end

%compare energies and wave functions with the exact results; 


fprintf('E = %20.10f\n',E);
fprintf ('dE = %20.10f\n', E-E(1));

%fprintf(file, '%g\t%.10f\t%.10f\n',alpha, E(1), E(2));
fprintf(file, '%20.10f', alpha);
fprintf(file, '%20.10f', E);
fprintf(file, '\n');

%end
%keyboard;

fclose(file);
%end

