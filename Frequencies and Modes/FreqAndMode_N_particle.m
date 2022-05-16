clc
clear all

disp('Frequencies and Modes for N particle systems')

%Loading the equilbirium positions and Alpha values:
%Data = load('Eq_Pos_eta_20_particles_1');
Data = load('Eq_Pos_eta_20_particles_1.mat');

% Distributing the alpha and the position values:
SubData = Data.eqpos;

Alpha   = SubData(end, :); 
Eq_Pos  = SubData(1:end-1, :);      %Equilibrium positions

disp(['\alpha ranges from ' num2str(Alpha(1)) ' to ' num2str(Alpha(end)) ' with ' num2str(length(Alpha)) ' values.'])

% Constant of the Build:
N       = length(Eq_Pos(:, 1));     %Particle Number
eta     = 20;                       % Dimensionless Coulomb int. strength [-]
ld      = 161.07;                   % Length unit [nm]
E_p     = 0.478;                    % Energy unit [meV]

% Building the Potential and Interaction Matrix
    % V is the potential matrix and U is the interaction

% Declaring the frequency matrix and the EigMode matrix:
Frequencies = zeros(length(Alpha), N);
EigModes    = zeros(length(Alpha), N, N);

for alpha_i = 1:length(Alpha)               % loop: for alpha values

    V = [];             %This will hold the second derivatives of the quartic potential
    U_diag      = zeros(N, 1);
    U_offdiag   = zeros(N, N);
    for part_i = 1:N     % loop: for particles

        V(part_i) = (3 * Eq_Pos(part_i, alpha_i)^2 + Alpha(alpha_i));
        
        for part_j = 1:N
            if part_i == part_j

            else
                U_diag(part_i)              = U_diag(part_i) + eta * (2/abs((Eq_Pos(part_i, alpha_i) - Eq_Pos(part_j, alpha_i))^3));
                U_offdiag(part_i, part_j)   = - eta * (2/abs((Eq_Pos(part_i, alpha_i) - Eq_Pos(part_j, alpha_i))^3));
            end
        end
    end

    % matrix elements are calculated
    % Constructing the Matrix

    for i = 1:N
        for j = 1:N
            if i == j
                Matrix(i, j) = V(i) + U_diag(i);
            else
                Matrix(i, j) = U_offdiag(i, j);
            end
        end
    end

    % Egzakt diagonalization:
    [ModVec, EigVal]        = eig(Matrix);
    Frequencies(alpha_i, :) = diag(EigVal);
    EigModes(alpha_i, :, :) = ModVec;

end
disp('Done!')


figure(1)
clf(figure(1))
hold on
for n = 1:N
    CurveName = [num2str(n) ' \omega'];
    plot(Alpha, sqrt(Frequencies(:, n)), '.', 'DisplayName', CurveName)
end
minimum = min(min(Frequencies));
[x, y] = find(Frequencies == minimum);
xline(Alpha(x), 'DisplayName', ['\alpha_c \sim ' num2str(Alpha(x))])
grid on
legend
xlabel('\alpha', 'FontSize', 22)
ylabel('\omega_i', 'FontSize', 22)

hold off

FileName = ['Omega_sq_Freqs_N_' num2str(N)];
save(FileName, 'EigModes')

