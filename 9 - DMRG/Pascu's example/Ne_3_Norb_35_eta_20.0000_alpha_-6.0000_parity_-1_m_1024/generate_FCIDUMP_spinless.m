Ne = 3;
Norb = 2;

N=Norb;%total number of orbitals.

xi = 0.3162;
%xi = 0.50;
%xi =0.1;  % this the effective screening. 

alpha = -6.0000;
p = 0.0000;
eta = 20.0000;
parity = -1;


X=80; 
X1 = 10; 

% lengthscale for the HO lenght
ld =0.65; 

nfft = 2^14;%number of points in the fft; 
dxn = X/nfft; 
Fs = 1/dxn; 
df = Fs/nfft;

xn = linspace(-X/2+0.5*dxn, X/2-0.5*dxn ,nfft);
fn = linspace( - Fs/2+0.5*df/2 , Fs/2 - 0.5*df , nfft );
omegan = 2*pi*fn;

%Fourier transform of the kernel ;
K_x= 1./sqrt(xn.^2+xi^2);
K_fft = 1/sqrt(2*pi)*fftshift(fft(ifftshift(K_x)))*dxn;


clasic_pos = linspace(-X1, X1, N);

error = 1e-10;
intervals = 8001; 
x = linspace(-X/2, X/2, intervals);
x1=x(1:end-1);
x2= x(2:end-1);

dx = (x(end)-x(1))/(intervals-1);
dx2=dx*dx;

V_pot = 1/4*x.^4 +1/2*(alpha)*x.^2-p*x;

V_min = min(V_pot);
V_pot = V_pot-V_min;

HO_WF=@(n, z, x0) interp1(x, 1./sqrt(2^(n-1)*factorial(n-1))*(1/(pi*ld^2))^0.25*exp(-(x-x0).^2/(2*ld^2)).*hermiteH(n-1,(x-x0)/ld), z, 'spline');
n=1;

WF = zeros(N, numel(x));

parity_list = [];

fprintf ('Construct wave functions\n');
  n=1;
  k=0;

fprintf('n=%d\n',n);

while n<=N-1
  k=k+1;
 
    WF(n,:)=0.5*(HO_WF(1, x, clasic_pos(k))+HO_WF(1, x, clasic_pos(N-k+1))) ;
    parity_list = [parity_list,1];
    n=n+1
     
 
   WF(n,:)=0.5*(HO_WF(1, x, clasic_pos(k))-HO_WF(1, x, clasic_pos(N-k+1))) ;
    parity_list = [parity_list,-1];
    n=n+1

end 

 WF(N,:)=HO_WF(1, x, clasic_pos((N+1)/2)) ;
 parity_list = [parity_list,1];

fprintf('Gram Schmidt Orthogonalization\n');

%Now we construct an orthogonal set using the Gram-Schmidt transformation;

WF = WF';
Q = zeros(numel(x), N);

[m,n] = size(WF);
% compute QR using Gram-Schmidt
for j = 1:n
    v = WF(:,j);
    for i=1:j-1
        R(i,j) = Q(:,i)'*WF(:,j);
        v = v - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end

WF=WF';
WF1=Q';

%proper normalization;

for j=1:N
norm1 = trapz(x, WF1(j,:).^2);
WF1(j,:)=WF1(j,:)/sqrt(norm1);
end

fprintf('Interpolate\n');

psi=@(n, z) interp1(x, WF1(n,:), z, 'spline');
%first derivative
psi1=@(n,z) interp1(x1, diff(WF1(n,:))./dx, z, 'spline');
%second derivative
psi2=@(n,z) interp1(x2, diff(diff(WF1(n,:)))./dx2, z, 'spline');


fprintf('Compute potentials\n');

%two body Coulomb integrals.
V = zeros(N, N, N, N);

for i=1:N
    tic;
    fprintf('i=%d\n', i);
    for j=1:N
        for k=1:N
            for l=1:N
               % tic
                psi_ij = psi(i, xn).*psi(j, xn);
                psi_ij_fft = 1/sqrt(2*pi)*fftshift(fft(ifftshift(psi_ij)))*dxn;
                psi_ij_fft = fliplr(psi_ij_fft);
                
                psi_kl = psi(k, xn).*psi(l, xn);
                psi_kl_fft = 1/sqrt(2*pi)*fftshift(fft(ifftshift(psi_kl)))*dxn;
               
             
                V(i, j, k, l) = V(i, j, k, l)+trapz(omegan, psi_ij_fft.*K_fft.*psi_kl_fft);
                V(i,j,k,l) = V(i,j,k,l)*sqrt(2*pi);
                
            end             
        end
    end
    toc;
end

%on site energy corresponding to diagonal term;
t=zeros(N);
t0=zeros(N);
t1 = zeros(N);
t2 = zeros(N);
t4 = zeros(N);

for i=1:N
    for j=1:N
       % integrand=@(z) psi(i, z).*psi2(j,z);
       % int = quad(integrand, -X/2, X/2);
        t0(i, j)= trapz(x,psi(i, x).*psi2(j,x) );
    end
end
%how compute the matrix elements of the z

for i=1:N
    for j=1:N
        %integrand=@(z) psi(i, z).*(z).*psi(j,z);
        %int = quad(integrand, -X/2, X/2);
        t1(i, j)= trapz(x, psi(i, x).*(x).*psi(j,x));
    end
end


for i=1:N
    for j=1:N
        %integrand=@(z) psi(i, z).*(z.^2).*psi(j,z);
        %int = quad(integrand, -X/2, X/2);
        t2(i, j)= trapz(x,psi(i, x).*(x.^2).*psi(j,x) );
    end
end


for i=1:N
    for j=1:N
        %integrand=@(z) psi(i, z).*(z.^4).*psi(j,z);
        %int = quad(integrand, -X/2, X/2);
        t4(i, j)= trapz(x,psi(i, x).*(x.^4).*psi(j,x));
    end
end


%t= -0.5*t0-p*t1+0.5*(kappa)*t2;
t= -0.5*t0-p*t1+0.5*(alpha)*t2+0.25*t4;

t=0.5*(t+t');

%here we save the data to the file.

filename = sprintf('FCIDUMP.dat');
file = fopen(filename, 'w');
fprintf(file, '&FCI NORB= %d,NELEC= %d,MS2= 0,\n',N, Ne);
fprintf(file, 'ORBSYM=%s\n', repmat('1,', 1, N));
fprintf(file, 'ISYM=1\n');
fprintf(file, '&END\n');

for i=1:N
    for j=1:N
        for k=1:N
            for l=1:N
                if abs(real(0.5*eta*V(i, j, k, l)))>error
                    fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n', 0.5*real(eta*V(i, j, k, l)), i, j,k,l );
                end
            end
        end
    end
end



for i=1:N
    for j=1:N
        if abs(t(i, j))>error
            fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n',t(i,j), i, j, 0, 0 );
        end
    end
end

fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n',0, 0, 0, 0, 0 );
fclose(file);


model_filename = 'model.input.wc_spinless'; 

%% Reading the template model file.
% Replace the data with the new variables.

fid = fopen(model_filename,'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;

while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);



fprintf('writing %s file\n', model_filename);
% Replace the necessary variables
fid = fopen(model_filename, 'w');



for i = 1:numel(A)
    if ischar(A{i})
        
      
        
        if startsWith(A{i}, 'CL')
            A{i} = sprintf('CL = %d;', N);
        end

       
        if startsWith(A{i}, 'par_i')
            A{i} = sprintf('par_i = [ %s ];', num2str(parity_list));
        end
        
        order = [1:N];
        if startsWith(A{i}, 'order')
            A{i} = sprintf('order = [ %s ];', num2str(order));
        end
        
    end
end

% rewrite the whole file

for i = 1:numel(A)-1
    fprintf(fid,'%s\n', A{i});
end
fclose(fid);


