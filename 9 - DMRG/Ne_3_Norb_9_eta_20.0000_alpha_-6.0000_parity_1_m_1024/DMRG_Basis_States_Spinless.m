clc
clear all

% Physicial constants:
Screening   = 10^-1;    % Effective coulomb screening
Alpha       = 1;        % Coeff of x^2
Beta        = 1;        % Coeff of x^4
Epsilon     = 0;        % Coeff of x^1
Eta         = 20;       % Coulomb dimensionless strength

Ls_d    = 160;      % Length Scale of the potential [nm]
L_HO    = 0.65;     % Length scale of the Hamornic oscillator

NumberOf_Electrons  = 3;
NumberOf_Orbitals   = 9;

NumberOf_Tot_Orbitals = NumberOf_Orbitals;  %total number of orbitals.

xi = 0.3162; 

Alpha       = -6;
Epsilon     = 0;
Eta         = 20;
Parity      = 1;

Length_X   = 80; 
Length_X1  = 10; 

% lengthscale for the HO lenght
L_d = 0.65; 

NumberOf_FFTPoints  = 2^14;%number of points in the fft; 
dx_n                = Length_X/NumberOf_FFTPoints; 
Length_Fourier      = 1/dx_n; 
df                  = Length_Fourier/NumberOf_FFTPoints;

Points_x        = linspace(-(Length_X/2) + (0.5 * dx_n),        (Length_X/2) - (0.5 * dx_n),    NumberOf_FFTPoints);
Points_fn       = linspace( -(Length_Fourier/2) + (0.5 * df)/2, (Length_Fourier/2) - (0.5*df),  NumberOf_FFTPoints);
Points_omegan   = 2 * pi * Points_fn;

%Fourier transform of the kernel ;
Kernel_x    = 1./sqrt(Points_x.^2 + xi^2);
Kernel_fft  = 1/sqrt(2 * pi) * fftshift(fft(ifftshift(Kernel_x))) * dx_n;

classic_pos = linspace(-Length_X1, Length_X1, NumberOf_Tot_Orbitals);

Error       = 1e-10;
Intervals   = 8001; 
x_0th_ord   = linspace(-Length_X/2, Length_X/2, Intervals);
x_1st_ord   = x_0th_ord(1:end-1);
x_2nd_ord   = x_0th_ord(2:end-1);

dx      = (x_0th_ord(end) - x_0th_ord(1))/(Intervals - 1);
dx2     = dx * dx;

V_potential = 1/4 * x_0th_ord.^4 + 1/2 * (Alpha) * x_0th_ord.^2 - Epsilon * x_0th_ord;
V_minimum   = min(V_potential);
V_potential = V_potential - V_minimum;

HarmonicOscillator_WF=@(n, z, x0) interp1(x_0th_ord, 1./sqrt(2^(n-1)*factorial(n-1))*(1/(pi*L_d^2))^0.25*exp(-(x_0th_ord-x0).^2/(2*L_d^2)).*hermiteH(n-1,(x_0th_ord-x0)/L_d), z, 'spline');
QuantumNumber_n = 1;

WaveFunction = zeros(NumberOf_Tot_Orbitals, numel(x_0th_ord));

parity_list = [];

fprintf ('Construct wave functions\n');
  QuantumNumber_n=1;
  k=0;

fprintf('n=%d\n',QuantumNumber_n);

while QuantumNumber_n <= NumberOf_Tot_Orbitals - 1
    k = k + 1;
 
    WaveFunction(QuantumNumber_n,:) = 0.5 * (HarmonicOscillator_WF(1, x_0th_ord, classic_pos(k)) + HarmonicOscillator_WF(1, x_0th_ord, classic_pos(NumberOf_Tot_Orbitals - k + 1))) ;
    parity_list                     = [parity_list,1];
    QuantumNumber_n                 = QuantumNumber_n + 1
     
 
    WaveFunction(QuantumNumber_n,:) = 0.5 * (HarmonicOscillator_WF(1, x_0th_ord, classic_pos(k)) - HarmonicOscillator_WF(1, x_0th_ord, classic_pos(NumberOf_Tot_Orbitals - k + 1))) ;
    parity_list                     = [parity_list,-1];
    QuantumNumber_n                 = QuantumNumber_n + 1

end 

 WaveFunction(NumberOf_Tot_Orbitals,:)  = HarmonicOscillator_WF(1, x_0th_ord, classic_pos((NumberOf_Tot_Orbitals + 1)/2)) ;
 parity_list                            = [parity_list,1];

fprintf('Gram Schmidt Orthogonalization\n');

%Now we construct an orthogonal set using the Gram-Schmidt transformation;

WaveFunction    = WaveFunction';
Q               = zeros(numel(x_0th_ord), NumberOf_Tot_Orbitals);

[m, QuantumNumber_n] = size(WaveFunction);
% compute QR using Gram-Schmidt
for j = 1:QuantumNumber_n
    v = WaveFunction(:,j);
    for index = 1:j-1
        R(index,j)  = Q(:,index)' * WaveFunction(:, j);
        v       = v - R(index, j)*Q(:, index);
    end
    R(j, j) = norm(v);
    Q(:, j) = v/R(j, j);
end

WaveFunction = WaveFunction';
WaveFunction_1 = Q';

%proper normalization;
for j = 1:NumberOf_Tot_Orbitals
    norm1               = trapz(x_0th_ord, WaveFunction_1(j,:).^2);
    WaveFunction_1(j,:) = WaveFunction_1(j,:)/sqrt(norm1);
end

fprintf('Interpolate\n');

psi_0th_ord = @(n, z) interp1(x_0th_ord, WaveFunction_1(n,:), z, 'spline');

% first derivative:
psi_1st_ord = @(n,z) interp1(x_1st_ord, diff(WaveFunction_1(n,:))./dx, z, 'spline');

%second derivative:
psi_2nd_ord = @(n,z) interp1(x_2nd_ord, diff(diff(WaveFunction_1(n,:)))./dx2, z, 'spline');

fprintf('Compute potentials\n');

%two body Coulomb integrals:
V = zeros(NumberOf_Tot_Orbitals, NumberOf_Tot_Orbitals, NumberOf_Tot_Orbitals, NumberOf_Tot_Orbitals);

for index = 1:NumberOf_Tot_Orbitals
    tic;
    fprintf('i=%d\n', index);
    for j = 1:NumberOf_Tot_Orbitals
        for k = 1:NumberOf_Tot_Orbitals
            for l = 1:NumberOf_Tot_Orbitals
               % tic
                psi_ij      = psi_0th_ord(index, Points_x) .* psi_0th_ord(j, Points_x);
                psi_ij_fft  = 1/sqrt(2 * pi) * fftshift(fft(ifftshift(psi_ij))) * dx_n;
                psi_ij_fft  = fliplr(psi_ij_fft);
                
                psi_kl      = psi_0th_ord(k, Points_x) .* psi_0th_ord(l, Points_x);
                psi_kl_fft  = 1/sqrt(2 * pi) * fftshift(fft(ifftshift(psi_kl))) * dx_n;
               
                V(index, j, k, l)   = V(index, j, k, l) + trapz(Points_omegan, psi_ij_fft .* Kernel_fft .* psi_kl_fft);
                V(index, j, k, l)   = V(index, j, k, l) * sqrt(2 * pi);
            end             
        end
    end
    toc;
end

%on site energy corresponding to diagonal term;
t   = zeros(NumberOf_Tot_Orbitals);
t0  = zeros(NumberOf_Tot_Orbitals);
t1  = zeros(NumberOf_Tot_Orbitals);
t2  = zeros(NumberOf_Tot_Orbitals);
t4  = zeros(NumberOf_Tot_Orbitals);

for index = 1:NumberOf_Tot_Orbitals
    for j = 1:NumberOf_Tot_Orbitals
        t0(index, j) = trapz(x_0th_ord, psi_0th_ord(index, x_0th_ord) .* psi_2nd_ord(j, x_0th_ord));
    end
end

%how compute the matrix elements of the z
for index = 1:NumberOf_Tot_Orbitals
    for j = 1:NumberOf_Tot_Orbitals
        t1(index, j) = trapz(x_0th_ord, psi_0th_ord(index, x_0th_ord) .* (x_0th_ord) .* psi_0th_ord(j, x_0th_ord));
    end
end

for index = 1:NumberOf_Tot_Orbitals
    for j = 1:NumberOf_Tot_Orbitals
        t2(index, j) = trapz(x_0th_ord, psi_0th_ord(index, x_0th_ord) .* (x_0th_ord.^2) .* psi_0th_ord(j, x_0th_ord));
    end
end

for index = 1:NumberOf_Tot_Orbitals
    for j = 1:NumberOf_Tot_Orbitals
        t4(index, j) = trapz(x_0th_ord, psi_0th_ord(index, x_0th_ord) .* (x_0th_ord.^4) .* psi_0th_ord(j, x_0th_ord));
    end
end

t = -0.5 * t0 - Epsilon * t1 + 0.5 * (Alpha) * t2 + 0.25 * t4;
t = 0.5 * (t + t');

% Here we save the data to the file.
filename    = sprintf('FCIDUMP.dat');
file        = fopen(filename, 'w');
fprintf(file, '&FCI NORB= %d,NELEC= %d,MS2= 0,\n',NumberOf_Tot_Orbitals, NumberOf_Electrons);
fprintf(file, 'ORBSYM=%s\n', repmat('1,', 1, NumberOf_Tot_Orbitals));
fprintf(file, 'ISYM=1\n');
fprintf(file, '&END\n');

for index = 1:NumberOf_Tot_Orbitals
    for j = 1:NumberOf_Tot_Orbitals
        for k = 1:NumberOf_Tot_Orbitals
            for l = 1:NumberOf_Tot_Orbitals
                if abs(real(0.5 * Eta * V(index, j, k, l)))>Error
                    fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n', 0.5 * real(Eta * V(index, j, k, l)), index, j, k, l);
                end
            end
        end
    end
end



for index = 1:NumberOf_Tot_Orbitals
    for j = 1:NumberOf_Tot_Orbitals
        if abs(t(index, j))>Error
            fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n',t(index,j), index, j, 0, 0 );
        end
    end
end

fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n',0, 0, 0, 0, 0 );
fclose(file);


model_filename = 'model.input.wc_spinless'; 

%% Reading the template model file.

% Replace the data with the new variables.
fid         = fopen(model_filename,'r');
index       = 1;
tline       = fgetl(fid);
A{index}    = tline;

while ischar(tline)
    index       = index + 1;
    tline       = fgetl(fid);
    A{index}    = tline;
end
fclose(fid);



fprintf('writing %s file\n', model_filename);
% Replace the necessary variables
fid = fopen(model_filename, 'w');



for index = 1:numel(A)
    if ischar(A{index})
        if startsWith(A{index}, 'CL')
            A{index} = sprintf('CL = %d;', NumberOf_Tot_Orbitals);
        end    

        if startsWith(A{index}, 'par_i')
            A{index} = sprintf('par_i = [ %s ];', num2str(parity_list));
        end

        order = [1:NumberOf_Tot_Orbitals];

        if startsWith(A{index}, 'order')
            A{index} = sprintf('order = [ %s ];', num2str(order));
        end
    end
end

% rewrite the whole file
for index = 1:numel(A)-1
    fprintf(fid, '%s\n', A{index});
end
fclose(fid);


